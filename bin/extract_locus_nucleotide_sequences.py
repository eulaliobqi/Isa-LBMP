# -*- coding: utf-8 -*-
"""Extrai a sequencia de nucleotideo real das 6 loci candidatas confirmadas.

Para locus_001/005/009/015 (que tem modelo genico do miniprot): concatena
os CDS na ordem genomica e traduz de uma vez so (MESMA logica de
bin/translate_miniprot_gff.py -- NAO traduzir exon a exon, miniprot tolera
frameshift dentro de um exon isolado). O resultado e validado comparando a
traducao com a proteina ja publicada em candidate_proteins_confirmed.fasta;
se qualquer uma das 4 falhar, o script aborta sem escrever nada.

Para locus_002/012 (sem modelo de exon, so evidencia TBLASTN): extrai o
span genomico bruto inteiro. Isso NAO e um CDS limpo -- pode conter introns
e/ou UTR -- e e rotulado como tal no header e no log.

Roda de dentro do repo:
    python bin/extract_locus_nucleotide_sequences.py
"""
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
GENOME = ROOT / "genoma.fa"
GFF3 = ROOT / "results" / "results" / "results" / "miniprot_hits.gff3"
FINAL_TABLE = ROOT / "results" / "results" / "results" / "final_candidates_table.tsv"
PROTEIN_FASTA = ROOT / "results" / "results" / "results" / "candidate_proteins_confirmed.fasta"
OUT_DIR = ROOT / "candidate_sequences"

# locus_id -> miniprot mRNA ID (so as 4 loci com modelo genico)
LOCUS_TO_MRNA = {
    "locus_001": "MP000001",
    "locus_005": "MP000003",
    "locus_009": "MP000004",
    "locus_015": "MP000005",
}
# loci sem modelo de exon (so span genomico bruto)
RAW_SPAN_LOCI = {"locus_002", "locus_012"}


def parse_gff3_cds(path, wanted_mrna_ids):
    mrnas = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, feat, start, end, _score, strand, _phase, attrs = fields
            attr = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
            if feat == "mRNA":
                mrna_id = attr.get("ID")
                if mrna_id in wanted_mrna_ids:
                    mrnas.setdefault(mrna_id, {"chrom": chrom, "strand": strand, "cds": []})
            elif feat == "CDS":
                parent = attr.get("Parent")
                if parent in wanted_mrna_ids:
                    mrnas.setdefault(parent, {"chrom": chrom, "strand": strand, "cds": []})
                    mrnas[parent]["cds"].append((int(start), int(end)))
    return mrnas


def parse_final_table(path, wanted_loci):
    loci = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            locus_id = fields[idx["locus_id"]]
            if locus_id in wanted_loci:
                loci[locus_id] = {
                    "chrom": fields[idx["chrom"]],
                    "start": int(fields[idx["start"]]),
                    "end": int(fields[idx["end"]]),
                    "strand": fields[idx["strand"]],
                }
    return loci


def load_known_proteins(path):
    proteins = {}
    for rec in SeqIO.parse(path, "fasta"):
        locus_id = rec.id.split("|")[0]
        proteins[locus_id] = str(rec.seq)
    return proteins


def main():
    if not GENOME.exists():
        print(f"ERRO: genoma nao encontrado em {GENOME}", file=sys.stderr)
        sys.exit(1)

    print("Indexando genoma...", file=sys.stderr)
    genome = SeqIO.index(str(GENOME), "fasta")

    known_proteins = load_known_proteins(PROTEIN_FASTA)

    # ---- 4 loci com modelo genico (miniprot) -------------------------------
    mrna_ids = set(LOCUS_TO_MRNA.values())
    mrnas = parse_gff3_cds(GFF3, mrna_ids)

    validated = {}
    log_lines = []
    all_ok = True

    for locus_id, mrna_id in LOCUS_TO_MRNA.items():
        info = mrnas.get(mrna_id)
        if not info or not info["cds"]:
            print(f"ERRO: nao achei CDS de {mrna_id} ({locus_id}) no GFF3", file=sys.stderr)
            all_ok = False
            continue
        chrom = info["chrom"]
        if chrom not in genome:
            print(f"ERRO: cromossomo {chrom} nao esta no genoma baixado", file=sys.stderr)
            all_ok = False
            continue
        cds_sorted = sorted(info["cds"], key=lambda iv: iv[0])
        nt_seq = "".join(str(genome[chrom].seq[s - 1:e]) for s, e in cds_sorted)
        seq_obj = Seq(nt_seq)
        if info["strand"] == "-":
            seq_obj = seq_obj.reverse_complement()
            nt_seq = str(seq_obj)
        translated = str(seq_obj.translate(to_stop=False)).rstrip("*")

        known = known_proteins.get(locus_id, "")
        status = "PASS" if translated == known else "FAIL"
        if status == "FAIL":
            all_ok = False
        log_lines.append(
            f"{locus_id} ({mrna_id}): {status} | CDS={len(nt_seq)}bp "
            f"| traduzido={len(translated)}aa | esperado={len(known)}aa"
        )
        if status == "FAIL":
            log_lines.append(f"  traduzido : {translated}")
            log_lines.append(f"  esperado  : {known}")

        span = f"{cds_sorted[0][0]}-{cds_sorted[-1][1]}"
        header = (
            f"{locus_id} CDS_reconstructed miniprot_model={len(cds_sorted)}_exons "
            f"{chrom}:{span}({info['strand']})"
        )
        validated[locus_id] = (header, nt_seq)

    if not all_ok:
        print("\nVALIDACAO FALHOU -- abortando sem escrever candidate_sequences/", file=sys.stderr)
        print("\n".join(log_lines), file=sys.stderr)
        sys.exit(1)

    print("Validacao das 4 loci com modelo genico: todas PASS.", file=sys.stderr)

    # ---- 2 loci sem modelo de exon (span bruto) ----------------------------
    raw_loci = parse_final_table(FINAL_TABLE, RAW_SPAN_LOCI)
    raw_entries = {}
    for locus_id, info in raw_loci.items():
        chrom = info["chrom"]
        if chrom not in genome:
            print(f"ERRO: cromossomo {chrom} nao esta no genoma baixado", file=sys.stderr)
            sys.exit(1)
        nt_seq = str(genome[chrom].seq[info["start"] - 1:info["end"]])
        if info["strand"] == "-":
            nt_seq = str(Seq(nt_seq).reverse_complement())
        header = (
            f"{locus_id} RAW_GENOMIC_SPAN exon_structure=UNKNOWN "
            f"pode_incluir_introns_ou_UTR {chrom}:{info['start']}-{info['end']}"
            f"({info['strand']}) length={len(nt_seq)}bp"
        )
        raw_entries[locus_id] = (header, nt_seq)
        log_lines.append(
            f"{locus_id}: INFORMATIVO (sem validacao possivel, sem modelo de exon) "
            f"| span_bruto={len(nt_seq)}bp"
        )

    # ---- escreve saidas -----------------------------------------------------
    OUT_DIR.mkdir(exist_ok=True)

    with open(OUT_DIR / "nucleotide_sequences.fasta", "w", encoding="utf-8") as out:
        for locus_id in ["locus_001", "locus_002", "locus_005", "locus_009", "locus_012", "locus_015"]:
            if locus_id in validated:
                header, seq = validated[locus_id]
            else:
                header, seq = raw_entries[locus_id]
            out.write(f">{header}\n{seq}\n")

    with open(OUT_DIR / "validation_log.txt", "w", encoding="utf-8") as out:
        out.write("Log de validacao -- extracao de nucleotideo (candidate_sequences/)\n")
        out.write("=" * 70 + "\n")
        out.write("\n".join(log_lines) + "\n")

    # copia protein fasta (fonte em results/ e git-ignored)
    with open(PROTEIN_FASTA, encoding="utf-8") as src, \
         open(OUT_DIR / "protein_sequences.fasta", "w", encoding="utf-8") as dst:
        dst.write(src.read())

    print(f"OK -- escrito em {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
