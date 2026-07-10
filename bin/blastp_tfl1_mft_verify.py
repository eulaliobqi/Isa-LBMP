"""Verifica por BLASTP (nao por nome de gene, que colidiu para RCN1 em arroz) os
homologos mais proximos de TFL1 e MFT de Arabidopsis em arroz e milho.

Uso: python bin/blastp_tfl1_mft_verify.py
"""
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "eulalio.santos@ufv.br"

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "tfl1_mft_verify.log"
OUT.parent.mkdir(exist_ok=True)

refs = {
    "TFL1_arabidopsis": "NP_196004.1",
    "MFT_arabidopsis": "NP_173250.1",
}


def log(fh, msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def main():
    with OUT.open("w", encoding="utf-8") as fh:
        seqs = {}
        for label, acc in refs.items():
            h = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
            rec = SeqIO.read(h, "fasta")
            h.close()
            seqs[label] = str(rec.seq)
            log(fh, f"{label} {acc} {len(rec.seq)} aa")

        for label, seq in seqs.items():
            for org in ["Oryza sativa", "Zea mays"]:
                log(fh, f"=== BLASTP {label} vs {org} — enviando ===")
                rh = NCBIWWW.qblast(
                    "blastp",
                    "nr",
                    seq,
                    hitlist_size=10,
                    expect=1e-10,
                    entrez_query=f"{org}[Organism]",
                )
                rec = NCBIXML.read(rh)
                rh.close()
                log(fh, f"=== BLASTP {label} vs {org} — resultado ===")
                for aln in rec.alignments[:8]:
                    hsp = aln.hsps[0]
                    pct = 100.0 * hsp.identities / hsp.align_length
                    log(
                        fh,
                        f"  {aln.accession}\t{pct:.1f}%\tevalue={hsp.expect:.1e}\t"
                        f"bits={hsp.bits:.1f}\t{aln.title[:90]}",
                    )
                time.sleep(1)
        log(fh, "CONCLUIDO")


if __name__ == "__main__":
    main()
