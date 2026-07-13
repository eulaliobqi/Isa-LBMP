#!/usr/bin/env python3
"""Visualiza árvore filogenética colorindo clados por função proteica.

Estratégia: desenho simples em estilo cladograma vertical.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent
TREEFILE = ROOT / "results" / "results" / "results" / "phylogeny" / "urochloa_ft_phylo.treefile"
OUT_PNG = ROOT / "results" / "results" / "results" / "figures" / "phylogeny_ete3_colored.png"
OUT_SVG = ROOT / "results" / "results" / "results" / "figures" / "phylogeny_ete3_colored.svg"

COLORS = {
    "FT-like": "#2e8b57",
    "TFL1-like": "#4682b4",
    "MFT-like": "#cd853f",
}

TAXA_CLASSIFICATION = {
    "NP_173250.1": "MFT-like",
    "NP_001410922.1": "MFT-like",
    "NP_001408118.1": "TFL1-like",
    "NP_001408117.1": "TFL1-like",
    "NP_001320342.1": "TFL1-like",
    "NP_001106249.1": "TFL1-like",
    "NP_001106241.1": "TFL1-like",
    "NP_001106247.1": "FT-like",
}


def classify_taxa(name):
    if not name:
        return "FT-like"
    if name in TAXA_CLASSIFICATION:
        return TAXA_CLASSIFICATION[name]
    if "locus_" in name:
        return "FT-like"
    return "FT-like"


def count_leaves(clade):
    """Conta folhas recursivamente."""
    if clade.is_terminal():
        return 1
    return sum(count_leaves(c) for c in clade.clades)


def draw_clade(clade, y, x_min, x_max, ax, depth=0, max_depth=1):
    """Desenha um clado recursivamente (layout HORIZONTAL). Retorna (y_usado, x_min, x_max)."""

    if clade.is_terminal():
        # Folha: ponto + label (nomes saem para baixo)
        x_center = (x_min + x_max) / 2
        func = classify_taxa(clade.name)
        color = COLORS[func]

        ax.scatter([x_center], [y], s=50, color=color, zorder=3, edgecolors="black", linewidth=0.5)

        # Simplificar nome: pegar só parte numérica
        label = clade.name if clade.name else "?"
        if "_" in label:
            # Extrair ID numérico (ex: "locus_001" → "001", "NP_123.1" → "123.1")
            parts = label.split("_")
            label = parts[-1]
        if len(label) > 20:
            label = label[:17] + "..."

        ax.text(x_center, y - 0.3, label, fontsize=8, ha="center", va="top", zorder=3,
                clip_on=False, family="monospace", rotation=90, weight="normal")

        return y - 0.5, x_center

    # Nó interno: desenha filhos e conecta
    if not clade.clades:
        return y, (x_min + x_max) / 2

    n_children = len(clade.clades)
    x_per_child = (x_max - x_min) / n_children

    child_ys = []
    child_xs = []

    for i, child in enumerate(clade.clades):
        child_x_min = x_min + i * x_per_child
        child_x_max = child_x_min + x_per_child

        child_y, child_x = draw_clade(child, y + 1, child_x_min, child_x_max, ax, depth + 1, max_depth)
        child_ys.append(child_y)
        child_xs.append(child_x)

    # Conectar filhos
    if child_xs:
        x_center = (min(child_xs) + max(child_xs)) / 2

        # Linha horizontal conectando todos os filhos
        ax.plot([min(child_xs), max(child_xs)], [y + 1, y + 1],
               color="#333333", linewidth=1.2, zorder=1)

        # Linhas verticais para cada filho
        for cx, cy in zip(child_xs, child_ys):
            ax.plot([cx, cx], [y + 1, cy],
                   color="#333333", linewidth=0.8, zorder=1)

        return max(child_ys), x_center

    return y, (x_min + x_max) / 2


def main():
    if not TREEFILE.exists():
        print(f"ERRO: {TREEFILE} nao existe", file=sys.stderr)
        sys.exit(1)

    print(f"Lendo arvore de {TREEFILE}")
    tree = Phylo.read(str(TREEFILE), "newick")

    # Contar folhas para dimensionar figura
    n_leaves = count_leaves(tree.root)
    print(f"Arvore tem {n_leaves} folhas")

    fig, ax = plt.subplots(figsize=(max(16, n_leaves * 0.4), 12), facecolor="white")

    print("Desenhando arvore filogenética...")
    draw_clade(tree.root, 0, 0, n_leaves, ax, max_depth=1)

    # Ajustes visuais
    ax.set_xlim(-1, n_leaves + 1)
    ax.set_ylim(-2, 10)
    ax.axis("off")

    # Títulos
    fig.suptitle("Urochloa FT Phylogeny", fontsize=18, fontweight="bold", y=0.98)
    ax.text(0.5, 0.95, "Colored by protein function (FT-like, TFL1-like, MFT-like)",
            fontsize=12, ha="center", transform=ax.transAxes, style="italic")

    # Legenda
    legend_handles = [Rectangle((0, 0), 1, 1, fc=color, ec="black", linewidth=1)
                     for color in COLORS.values()]
    ax.legend(legend_handles, COLORS.keys(), loc="upper left",
             fontsize=11, title="Protein Function", framealpha=0.95)

    print(f"Exportando PNG para {OUT_PNG}")
    fig.savefig(OUT_PNG, dpi=300, format="png", bbox_inches="tight", facecolor="white")

    print(f"Exportando SVG para {OUT_SVG}")
    fig.savefig(OUT_SVG, dpi=300, format="svg", bbox_inches="tight", facecolor="white")

    plt.close(fig)

    print("[OK] Arvore colorida salva com sucesso!")
    print(f"  PNG: {OUT_PNG}")
    print(f"  SVG: {OUT_SVG}")
    print("\nLegenda de cores:")
    for func, color in COLORS.items():
        print(f"  {func}: {color}")


if __name__ == "__main__":
    main()
