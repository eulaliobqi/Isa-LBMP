#!/usr/bin/env python3
"""Consolida loci_summary + pebp_confirmed + rbh_summary + árvore filogenética
em uma tabela final de candidatos e duas figuras resumo.

Ver Passo 8 do plano científico. n_evidencias (0-3) = domínio PEBP confirmado
+ RBH confirmado (qualquer espécie) + presença na filogenia.
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

try:
    from Bio import Phylo
except ImportError:
    Phylo = None


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--loci-summary", required=True)
    ap.add_argument("--pebp-confirmed", required=True)
    ap.add_argument("--rbh-summary", required=True)
    ap.add_argument("--treefile", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-heatmap", required=True)
    ap.add_argument("--out-diagram", required=True)
    return ap.parse_args()


def safe_read_tsv(path, **kwargs):
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def locus_ids_with_pebp(pebp_df):
    if pebp_df.empty or "target_name" not in pebp_df.columns:
        return set()
    return {str(t).split("|", 1)[0] for t in pebp_df["target_name"].dropna()}


def locus_ids_with_rbh(rbh_df):
    if rbh_df.empty or "reciprocal_confirmed" not in rbh_df.columns:
        return set()
    confirmed = rbh_df[rbh_df["reciprocal_confirmed"] == True]  # noqa: E712
    return set(confirmed["locus_id"].dropna())


def locus_ids_in_tree(treefile):
    if Phylo is None:
        print("AVISO: Biopython Bio.Phylo indisponível, pulando checagem de árvore.", file=sys.stderr)
        return set()
    try:
        tree = Phylo.read(treefile, "newick")
    except Exception as exc:  # noqa: BLE001
        print(f"AVISO: não foi possível ler {treefile}: {exc}", file=sys.stderr)
        return set()
    leaf_names = {leaf.name for leaf in tree.get_terminals() if leaf.name}
    return {name.split("|", 1)[0] for name in leaf_names if name.startswith("locus_")}


def build_heatmap(loci_df, out_path):
    if loci_df.empty:
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "Nenhum locus candidato", ha="center", va="center")
        plt.axis("off")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        return

    rows = []
    for _, row in loci_df.iterrows():
        genes = str(row.get("queries", "")).split(";")
        for gene in genes:
            gene = gene.strip()
            if gene:
                rows.append({"locus_id": row.locus_id, "gene": gene, "bitscore": row.get("best_bitscore") or 0})
    long_df = pd.DataFrame(rows)
    if long_df.empty:
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "Sem convergência de genes de referência", ha="center", va="center")
        plt.axis("off")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        return

    pivot = long_df.pivot_table(index="locus_id", columns="gene", values="bitscore", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(pivot.columns) + 2), max(3, 0.4 * len(pivot.index) + 1)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Bitscore por locus candidato x gene de referência convergente")
    fig.colorbar(im, ax=ax, label="bitscore (TBLASTN)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def build_diagram(loci_df, out_path):
    if loci_df.empty or "chrom" not in loci_df.columns:
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "Nenhum locus candidato", ha="center", va="center")
        plt.axis("off")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        return

    chroms = sorted(loci_df["chrom"].dropna().unique())
    fig, ax = plt.subplots(figsize=(10, max(2, 0.6 * len(chroms) + 1)))
    for i, chrom in enumerate(chroms):
        sub = loci_df[loci_df["chrom"] == chrom]
        ax.hlines(y=i, xmin=0, xmax=sub[["start", "end"]].max().max() * 1.02, color="lightgray", linewidth=1)
        for _, row in sub.iterrows():
            ax.plot([row.start, row.end], [i, i], linewidth=6, solid_capstyle="butt")
            ax.annotate(row.locus_id, (row.start, i), textcoords="offset points", xytext=(0, 6), fontsize=7)
    ax.set_yticks(range(len(chroms)))
    ax.set_yticklabels(chroms)
    ax.set_xlabel("Posição no cromossomo (pb)")
    ax.set_title("Loci candidatos por cromossomo")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()

    loci_df = safe_read_tsv(args.loci_summary)
    pebp_df = safe_read_tsv(args.pebp_confirmed)
    rbh_df = safe_read_tsv(args.rbh_summary)

    pebp_loci = locus_ids_with_pebp(pebp_df)
    rbh_loci = locus_ids_with_rbh(rbh_df)
    tree_loci = locus_ids_in_tree(args.treefile)

    if not loci_df.empty:
        loci_df["pebp_confirmed"] = loci_df["locus_id"].isin(pebp_loci)
        loci_df["rbh_confirmed"] = loci_df["locus_id"].isin(rbh_loci)
        loci_df["in_phylogeny"] = loci_df["locus_id"].isin(tree_loci)
        loci_df["n_evidencias"] = (
            loci_df["pebp_confirmed"].astype(int)
            + loci_df["rbh_confirmed"].astype(int)
            + loci_df["in_phylogeny"].astype(int)
        )
        loci_df = loci_df.sort_values(["n_evidencias", "best_bitscore"], ascending=[False, False])

    loci_df.to_csv(args.out_table, sep="\t", index=False)

    build_heatmap(loci_df, args.out_heatmap)
    build_diagram(loci_df, args.out_diagram)

    if not loci_df.empty:
        top = loci_df.iloc[0]
        print(
            f"Locus mais forte: {top.locus_id} ({top.chrom}:{top.start}-{top.end}), "
            f"n_evidencias={top.n_evidencias}/3, genes convergentes={top.queries}",
            file=sys.stderr,
        )
    else:
        print("Nenhum locus candidato encontrado.", file=sys.stderr)


if __name__ == "__main__":
    main()
