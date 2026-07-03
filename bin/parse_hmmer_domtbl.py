#!/usr/bin/env python3
"""Converte hmmsearch --domtblout em TSV limpo + FASTA só com domínio PEBP confirmado.

Ver Passo 5 do plano científico. Reaproveita o padrão de parsing de domtblout já
usado em Kerson-paper/analyses/06_domain_architecture.
"""
import argparse
import sys

import pandas as pd
from Bio import SeqIO

DOMTBL_COLS = [
    "target_name", "target_accession", "tlen", "query_name", "query_accession",
    "qlen", "evalue_full", "score_full", "bias_full", "dom_num", "dom_of",
    "c_evalue", "i_evalue", "score_dom", "bias_dom",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fasta", required=True)
    return ap.parse_args()


def parse_domtblout(path):
    rows = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split(None, len(DOMTBL_COLS))
                row = dict(zip(DOMTBL_COLS, fields[: len(DOMTBL_COLS)]))
                row["description"] = fields[len(DOMTBL_COLS)].strip() if len(fields) > len(DOMTBL_COLS) else ""
                rows.append(row)
    except FileNotFoundError:
        pass
    df = pd.DataFrame(rows)
    if not df.empty:
        for col in ["tlen", "qlen", "dom_num", "dom_of", "hmm_from", "hmm_to",
                    "ali_from", "ali_to", "env_from", "env_to"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        for col in ["evalue_full", "score_full", "bias_full", "c_evalue",
                    "i_evalue", "score_dom", "bias_dom", "acc"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def main():
    args = parse_args()
    df = parse_domtblout(args.domtbl)

    out_cols = [
        "target_name", "query_name", "i_evalue", "score_dom",
        "ali_from", "ali_to", "acc",
    ]
    df.reindex(columns=out_cols).to_csv(args.out_tsv, sep="\t", index=False)

    confirmed_ids = set(df["target_name"]) if not df.empty else set()
    print(f"Alvos confirmados com domínio PEBP: {len(confirmed_ids)}", file=sys.stderr)

    written = 0
    with open(args.out_fasta, "w", encoding="utf-8") as out_fh:
        for rec in SeqIO.parse(args.candidates, "fasta"):
            if rec.id in confirmed_ids:
                out_fh.write(f">{rec.id}\n{rec.seq}\n")
                written += 1
    print(f"Proteínas escritas em candidate_proteins_confirmed.fasta: {written}", file=sys.stderr)


if __name__ == "__main__":
    main()
