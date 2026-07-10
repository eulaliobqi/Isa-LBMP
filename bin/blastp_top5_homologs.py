"""Roda BLASTP remoto (NCBI nr, restrito ao genero Urochloa) para cada proteina em
ptn-arroz-milho-arabidopsis.txt e baixa as 5 sequencias distintas mais proximas
(excluindo a propria query) via efetch.

Uso: python bin/blastp_top5_homologs.py
Saida: results/blastp_homologs/<query_id>_top5.fasta (+ combined.fasta, summary.tsv)
"""
import re
import sys
import time
from pathlib import Path

from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "eulalio.santos@ufv.br"

ROOT = Path(__file__).resolve().parent.parent
INPUT_FASTA = ROOT / "ptn-arroz-milho-arabidopsis.txt"
OUTDIR = ROOT / "blastp_homologs"  # rastreado no git (dado de referencia curado, nao output de pipeline)
OUTDIR.mkdir(parents=True, exist_ok=True)

N_HITS_WANTED = 5
HITLIST_SIZE = 50  # margem para filtrar self-hits e duplicatas identicas
ENTREZ_QUERY = "Urochloa[Organism]"  # restringe BLASTP ao genero Urochloa
EXPECT = 10  # relaxado: banco filtrado por taxon e pequeno, poucos hits esperados


def parse_fasta(path):
    entries = []
    header, seq_lines = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        elif line.strip():
            seq_lines.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def accession_of(header):
    return header.split()[0]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_one(header, seq, idx, total):
    acc = accession_of(header)
    log(f"({idx}/{total}) BLASTP {acc} — enviando para NCBI nr (entrez_query={ENTREZ_QUERY})...")
    result_handle = NCBIWWW.qblast(
        "blastp",
        "nr",
        seq,
        hitlist_size=HITLIST_SIZE,
        expect=EXPECT,
        entrez_query=ENTREZ_QUERY,
    )
    log(f"({idx}/{total}) {acc} — resultado recebido, parseando XML...")
    record = NCBIXML.read(result_handle)
    result_handle.close()

    seen_seqs = set()
    picked = []
    for alignment in record.alignments:
        hit_acc = alignment.accession
        if hit_acc == acc or acc.split(".")[0] == hit_acc.split(".")[0]:
            continue  # pula a propria query
        hsp = alignment.hsps[0]
        pct_id = 100.0 * hsp.identities / hsp.align_length
        picked.append(
            {
                "hit_acc": hit_acc,
                "title": alignment.title,
                "pct_id": pct_id,
                "evalue": hsp.expect,
                "bitscore": hsp.bits,
            }
        )

    # dedup por sequencia real (baixando fasta) mantendo ordem de melhor score
    final = []
    for hit in picked:
        if len(final) >= N_HITS_WANTED:
            break
        log(f"  fetching {hit['hit_acc']} ...")
        try:
            handle = Entrez.efetch(
                db="protein", id=hit["hit_acc"], rettype="fasta", retmode="text"
            )
            fasta_text = handle.read()
            handle.close()
        except Exception as e:
            log(f"  falhou {hit['hit_acc']}: {e}")
            continue
        seq_body = "".join(
            l.strip() for l in fasta_text.splitlines() if not l.startswith(">")
        )
        if seq_body in seen_seqs or not seq_body:
            continue
        seen_seqs.add(seq_body)
        hit["fasta"] = fasta_text.strip() + "\n"
        final.append(hit)
        time.sleep(0.4)  # cortesia NCBI

    return acc, final


def main():
    entries = parse_fasta(INPUT_FASTA)
    log(f"{len(entries)} proteinas de query encontradas em {INPUT_FASTA.name}")

    combined_path = OUTDIR / "combined_top5.fasta"
    summary_path = OUTDIR / "summary.tsv"
    combined_lines = []
    summary_rows = [
        "query_acc\tquery_desc\thit_acc\thit_desc\tpct_identity\tevalue\tbitscore"
    ]

    for idx, (header, seq) in enumerate(entries, 1):
        acc, hits = run_one(header, seq, idx, len(entries))
        out_path = OUTDIR / f"{acc.replace('.', '_')}_top5.fasta"
        with out_path.open("w", encoding="utf-8") as fh:
            for h in hits:
                fh.write(h["fasta"])
                combined_lines.append(h["fasta"])
                summary_rows.append(
                    "\t".join(
                        [
                            acc,
                            header,
                            h["hit_acc"],
                            h["title"][:120],
                            f"{h['pct_id']:.1f}",
                            f"{h['evalue']:.2e}",
                            f"{h['bitscore']:.1f}",
                        ]
                    )
                )
        log(f"({idx}/{len(entries)}) {acc}: {len(hits)} homologos salvos em {out_path.name}")

    combined_path.write_text("".join(combined_lines), encoding="utf-8")
    summary_path.write_text("\n".join(summary_rows) + "\n", encoding="utf-8")
    log(f"Concluido. Combined: {combined_path}")
    log(f"Sumario: {summary_path}")


if __name__ == "__main__":
    main()
