#!/usr/bin/env python3
"""Traduz os modelos gênicos (mRNA/CDS) do GFF3 do miniprot para proteína.

NÃO usa gffread: miniprot tolera pequenos frameshifts dentro de um éxon
individual (o comprimento de um CDS isolado pode não ser múltiplo de 3),
mas a soma de todos os CDS de um mRNA É múltipla de 3. Traduzir exon a
exon respeitando a coluna "phase" isoladamente (abordagem do gffread)
produz proteína fora de frame; concatenar todos os CDS em ordem genômica
e traduzir de uma vez a partir da posição 0 reproduz corretamente o
alinhamento pretendido pelo miniprot.
"""
import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out-fasta", required=True)
    return ap.parse_args()


def parse_gff3_cds(path):
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
                if mrna_id:
                    mrnas.setdefault(mrna_id, {"chrom": chrom, "strand": strand, "cds": []})
            elif feat == "CDS":
                parent = attr.get("Parent")
                if not parent:
                    continue
                mrnas.setdefault(parent, {"chrom": chrom, "strand": strand, "cds": []})
                mrnas[parent]["cds"].append((int(start), int(end)))
    return mrnas


def main():
    args = parse_args()
    genome = SeqIO.index(args.genome, "fasta")
    mrnas = parse_gff3_cds(args.gff3)

    n_written = 0
    with open(args.out_fasta, "w", encoding="utf-8") as out:
        for mrna_id, info in mrnas.items():
            chrom = info["chrom"]
            if chrom not in genome or not info["cds"]:
                continue
            cds_sorted = sorted(info["cds"], key=lambda iv: iv[0])
            nt_seq = "".join(str(genome[chrom].seq[s - 1:e]) for s, e in cds_sorted)
            seq_obj = Seq(nt_seq)
            if info["strand"] == "-":
                seq_obj = seq_obj.reverse_complement()
            protein = str(seq_obj.translate(to_stop=False)).rstrip("*")
            if not protein:
                continue
            out.write(f">{mrna_id}\n{protein}\n")
            n_written += 1

    print(f"Proteínas traduzidas (concatenação de CDS, frame único): {n_written}", file=sys.stderr)


if __name__ == "__main__":
    main()
