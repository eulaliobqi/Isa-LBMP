#!/usr/bin/env python3
"""Consolida hits do TBLASTN e modelos gênicos do miniprot em loci candidatos.

Ver Passo 4 do plano científico:
C:\\Users\\eulal\\.claude\\plans\\objetivo-principal-usar-o-synthetic-peach.md
"""
import argparse
import re
import sys

import pandas as pd
from Bio import SeqIO

TBLASTN_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
    "qlen", "slen", "sseq",
]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tblastn", required=True)
    ap.add_argument("--miniprot-gff3", required=True)
    ap.add_argument("--miniprot-proteins", required=True)
    ap.add_argument("--cluster-distance-bp", type=int, default=20000)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-fasta", required=True)
    return ap.parse_args()


def load_tblastn(path):
    try:
        df = pd.read_csv(path, sep="\t", names=TBLASTN_COLS, header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=TBLASTN_COLS)
    if df.empty:
        return df
    df["strand"] = df.apply(lambda r: "+" if r.sstart <= r.send else "-", axis=1)
    df["gstart"] = df[["sstart", "send"]].min(axis=1)
    df["gend"] = df[["sstart", "send"]].max(axis=1)
    return df


def cluster_intervals(df, cluster_distance_bp):
    """Agrupa hits do TBLASTN por sseqid+strand, mesclando intervalos próximos."""
    loci = []
    if df.empty:
        return loci
    for (chrom, strand), g in df.sort_values("gstart").groupby(["sseqid", "strand"]):
        cur_start = cur_end = None
        cur_hits = []
        for _, row in g.iterrows():
            if cur_start is None or row.gstart - cur_end > cluster_distance_bp:
                if cur_hits:
                    loci.append((chrom, strand, cur_start, cur_end, cur_hits))
                cur_start, cur_end, cur_hits = row.gstart, row.gend, [row]
            else:
                cur_end = max(cur_end, row.gend)
                cur_hits.append(row)
        if cur_hits:
            loci.append((chrom, strand, cur_start, cur_end, cur_hits))
    return loci


def parse_miniprot_gff3(path):
    """Extrai features mRNA do GFF3 do miniprot: (chrom, strand, start, end, id, target)."""
    models = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "mRNA":
                    continue
                chrom, _src, _feat, start, end, _score, strand, _frame, attrs = fields
                attr_dict = dict(
                    kv.split("=", 1) for kv in attrs.split(";") if "=" in kv
                )
                mrna_id = attr_dict.get("ID", "")
                target = attr_dict.get("Target", "")
                query_name = target.split(" ")[0] if target else mrna_id
                models.append({
                    "chrom": chrom,
                    "strand": strand,
                    "start": int(start),
                    "end": int(end),
                    "id": mrna_id,
                    "query": query_name,
                })
    except FileNotFoundError:
        pass
    return models


def overlaps_or_near(a_start, a_end, b_start, b_end, max_gap):
    return not (a_end + max_gap < b_start or b_end + max_gap < a_start)


def main():
    args = parse_args()

    tblastn_df = load_tblastn(args.tblastn)
    tblastn_loci = cluster_intervals(tblastn_df, args.cluster_distance_bp)
    miniprot_models = parse_miniprot_gff3(args.miniprot_gff3)

    miniprot_proteins = {}
    try:
        for rec in SeqIO.parse(args.miniprot_proteins, "fasta"):
            miniprot_proteins[rec.id] = str(rec.seq)
    except FileNotFoundError:
        pass

    records = []

    # Loci originados do TBLASTN, tentando casar com um modelo do miniprot
    for chrom, strand, start, end, hits in tblastn_loci:
        queries = sorted({h.qseqid for h in hits})
        best_hit = max(hits, key=lambda h: h.bitscore)
        matched_model = None
        for m in miniprot_models:
            if m["chrom"] == chrom and overlaps_or_near(
                start, end, m["start"], m["end"], args.cluster_distance_bp
            ):
                matched_model = m
                break
        source = "tblastn+miniprot" if matched_model else "tblastn_only"
        predicted_id = matched_model["id"] if matched_model else None
        predicted_seq = miniprot_proteins.get(predicted_id) if predicted_id else None
        if predicted_seq is None:
            predicted_seq = best_hit.sseq.replace("-", "")
        records.append({
            "chrom": chrom,
            "strand": strand,
            "start": min(start, matched_model["start"]) if matched_model else start,
            "end": max(end, matched_model["end"]) if matched_model else end,
            "queries": ";".join(sorted(set(queries) | ({matched_model["query"]} if matched_model else set()))),
            "source": source,
            "n_queries_hit": len(queries),
            "best_bitscore": best_hit.bitscore,
            "predicted_protein_id": predicted_id or f"{chrom}_{start}-{end}_tblastn",
            "predicted_protein_seq": predicted_seq,
        })

    # Modelos do miniprot sem nenhum hit correspondente do TBLASTN (splicing capturou
    # algo que a triagem ampla perdeu, ex. exon curto/divergente)
    matched_chrom_ranges = [
        (r["chrom"], r["start"], r["end"]) for r in records
    ]
    for m in miniprot_models:
        already_matched = any(
            m["chrom"] == c and overlaps_or_near(m["start"], m["end"], s, e, args.cluster_distance_bp)
            for c, s, e in matched_chrom_ranges
        )
        if already_matched:
            continue
        predicted_seq = miniprot_proteins.get(m["id"], "")
        records.append({
            "chrom": m["chrom"],
            "strand": m["strand"],
            "start": m["start"],
            "end": m["end"],
            "queries": m["query"],
            "source": "miniprot_only",
            "n_queries_hit": 1,
            "best_bitscore": None,
            "predicted_protein_id": m["id"],
            "predicted_protein_seq": predicted_seq,
        })

    out_df = pd.DataFrame(records)
    if not out_df.empty:
        out_df = out_df.sort_values("n_queries_hit", ascending=False).reset_index(drop=True)
        out_df.insert(0, "locus_id", [f"locus_{i+1:03d}" for i in range(len(out_df))])
        out_df["predicted_protein_len"] = out_df["predicted_protein_seq"].str.len()

    cols = [
        "locus_id", "chrom", "start", "end", "strand", "queries", "source",
        "n_queries_hit", "best_bitscore", "predicted_protein_id", "predicted_protein_len",
    ]
    out_df.reindex(columns=cols).to_csv(args.out_summary, sep="\t", index=False)

    with open(args.out_fasta, "w", encoding="utf-8") as fh:
        for _, row in out_df.iterrows():
            seq = row.get("predicted_protein_seq") or ""
            if not seq:
                continue
            fh.write(f">{row.locus_id}|{row.predicted_protein_id}|source={row.source}\n{seq}\n")

    print(f"Loci consolidados: {len(out_df)}", file=sys.stderr)


if __name__ == "__main__":
    main()
