#!/usr/bin/env python3
"""Determina Reciprocal Best Hit (RBH) a partir dos TSVs fwd_<sp>.tsv / rev_<sp>.tsv
gerados pelo módulo RECIPROCAL_BEST_HIT.

RBH confirmado quando o melhor hit de volta (rev) é o próprio candidato original,
ou um dos genes de referência originais (query_proteins.fasta) — ver Passo 9 do
plano científico. IDs de candidatos seguem o padrão `locus_NNN|...` gerado por
bin/consolidate_loci.py; qualquer outro ID no BLAST db combinado é um dos genes
de referência originais.
"""
import argparse
import sys

import pandas as pd

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--out-tsv", required=True)
    return ap.parse_args()


def read_blast_tsv(path):
    try:
        df = pd.read_csv(path, sep="\t", names=BLAST_COLS, header=None)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=BLAST_COLS)
    return df


def best_hit_per_query(df):
    if df.empty:
        return {}
    idx = df.groupby("qseqid")["bitscore"].idxmax()
    best = df.loc[idx]
    return dict(zip(best["qseqid"], best["sseqid"]))


def is_candidate_id(seq_id):
    return seq_id.startswith("locus_")


def candidate_locus(seq_id):
    return seq_id.split("|", 1)[0]


def strip_id_wrapper(seq_id):
    """makeblastdb -parse_seqids faz o blastp reportar sseqid como
    "ref|ACCESSION|" (formato NCBI), mas blastdbcmd -entry_batch extrai (e o
    blastp reciproco lê de volta) o header ORIGINAL, accession puro -- sem
    normalizar os dois lados pro mesmo formato antes de cruzar, o lookup
    nunca bate mesmo em auto-hits triviais."""
    parts = str(seq_id).split("|")
    if len(parts) >= 3 and parts[1]:
        return parts[1]
    return seq_id


def main():
    args = parse_args()
    records = []

    for sp in args.species:
        fwd = read_blast_tsv(f"fwd_{sp}.tsv")
        rev = read_blast_tsv(f"rev_{sp}.tsv")
        fwd_best = best_hit_per_query(fwd)   # candidate_id -> species_protein_id
        rev_best_raw = best_hit_per_query(rev)   # species_protein_id -> combined_ref_id
        rev_best = {strip_id_wrapper(k): v for k, v in rev_best_raw.items()}

        for candidate_id, species_protein_id in fwd_best.items():
            rev_hit = rev_best.get(strip_id_wrapper(species_protein_id))
            reciprocal_confirmed = False
            if rev_hit is not None:
                if is_candidate_id(rev_hit):
                    reciprocal_confirmed = candidate_locus(rev_hit) == candidate_locus(candidate_id)
                else:
                    # hit de volta é um dos genes de referência originais (FT-like esperado,
                    # ou TFL1/MFT se o conjunto de referência tiver sido expandido — a
                    # interpretação fica em best_hit_gene para revisão manual)
                    reciprocal_confirmed = True

            fwd_row = fwd[(fwd.qseqid == candidate_id) & (fwd.sseqid == species_protein_id)].iloc[0]
            records.append({
                "locus_id": candidate_locus(candidate_id),
                "species": sp,
                "forward_hit_protein": species_protein_id,
                "forward_hit_evalue": fwd_row.evalue,
                "forward_hit_bitscore": fwd_row.bitscore,
                "best_hit_gene": rev_hit,
                "rev_hit_is_query": (rev_hit is not None and not is_candidate_id(rev_hit)),
                "reciprocal_confirmed": reciprocal_confirmed,
            })

    out_df = pd.DataFrame(records)
    out_df.to_csv(args.out_tsv, sep="\t", index=False)
    n_confirmed = int(out_df["reciprocal_confirmed"].sum()) if not out_df.empty else 0
    print(f"RBH confirmados: {n_confirmed} / {len(out_df)}", file=sys.stderr)


if __name__ == "__main__":
    main()
