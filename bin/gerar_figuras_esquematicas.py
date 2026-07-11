# -*- coding: utf-8 -*-
"""Gera 4 figuras esquematicas adicionais para o relatorio Word:
  A) fluxograma do pipeline (metodologia)
  B) estrutura genica: predicao (U. ruziziensis) vs. real com UTR (U. brizantha)
  C) diagrama de convergencia de evidencias para locus_001
  D) mapa dos 36 cromossomos do assembly de U. decumbens (subgenoma b/rd)

Paleta: skill dataviz (references/palette.md), ordem categorica fixa.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import matplotlib.font_manager as fm

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "results" / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# ---- paleta (dataviz skill, references/palette.md) ------------------------
BLUE = "#2a78d6"
AQUA = "#1baf7a"
YELLOW = "#eda100"
GREEN = "#008300"
VIOLET = "#4a3aa7"
RED = "#e34948"
MAGENTA = "#e87ba4"
ORANGE = "#eb6834"

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"
GOOD = "#0ca30c"

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["text.color"] = INK
plt.rcParams["axes.edgecolor"] = MUTED


# =============================================================================
# A) Fluxograma do pipeline
# =============================================================================
def fig_pipeline():
    stages = [
        ("1. BUILD_BLASTDB", "Indexa o genoma alvo", BLUE),
        ("2. TBLASTN", "Triagem ampla proteína × genoma\n(E-value 1e-5)", BLUE),
        ("2b. BLASTN_CDS", "Evidência cruzada CDS × genoma\n(E-value 1e-3)", BLUE),
        ("3. MINIPROT", "Estrutura gênica spliced-aware\n+ tradução", BLUE),
        ("4. CONSOLIDATE_LOCI", "Funde as 3 fontes em loci\ncandidatos (cluster <20 kb)", BLUE),
        ("5. HMMSEARCH_PEBP", "Confirma domínio PEBP\n(PF01161, HMMER)", AQUA),
        ("6. DOWNLOAD_REF_PROTEOMES", "Baixa proteomas de\narroz/milho/Arabidopsis", AQUA),
        ("7. RECIPROCAL_BEST_HIT", "BLASTP recíproco — confirma\nortologia (E-value 1e-10)", AQUA),
        ("8. PHYLOGENY", "MAFFT + trimAl + IQ-TREE\n(outgroup MFT-like)", VIOLET),
        ("9. REPORT", "Tabela final + figuras", VIOLET),
    ]

    fig, ax = plt.subplots(figsize=(9, 13))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(stages) * 2 + 1)
    ax.axis("off")

    box_w, box_h = 8.6, 1.5
    y = len(stages) * 2 - 0.5
    centers = []
    for title, desc, color in stages:
        cx, cy = 5, y
        centers.append((cx, cy))
        box = FancyBboxPatch(
            (cx - box_w / 2, cy - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.08,rounding_size=0.12",
            linewidth=1.4, edgecolor=color, facecolor=SURFACE, zorder=2,
        )
        ax.add_patch(box)
        # colored left tab (identity marker)
        tab = FancyBboxPatch(
            (cx - box_w / 2, cy - box_h / 2), 0.18, box_h,
            boxstyle="round,pad=0,rounding_size=0.06",
            linewidth=0, facecolor=color, zorder=3,
        )
        ax.add_patch(tab)
        ax.text(cx - box_w / 2 + 0.45, cy + 0.28, title, fontsize=11.5, fontweight="bold",
                 color=INK, ha="left", va="center", zorder=4)
        ax.text(cx - box_w / 2 + 0.45, cy - 0.32, desc, fontsize=9.5, color=INK2,
                 ha="left", va="center", zorder=4, linespacing=1.4)
        y -= 2

    for i in range(len(centers) - 1):
        x0, y0 = centers[i]
        x1, y1 = centers[i + 1]
        arrow = FancyArrowPatch(
            (x0, y0 - box_h / 2), (x1, y1 + box_h / 2),
            arrowstyle="-|>", mutation_scale=16, linewidth=1.6,
            color=MUTED, zorder=1,
        )
        ax.add_patch(arrow)

    # legenda de categorias
    handles = [
        mpatches.Patch(color=BLUE, label="Descoberta de loci"),
        mpatches.Patch(color=AQUA, label="Confirmação de ortologia"),
        mpatches.Patch(color=VIOLET, label="Filogenia e relatório"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
              ncol=1, frameon=False, fontsize=10)

    plt.tight_layout()
    out = FIG / "esquema_pipeline.png"
    plt.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close()
    print("salvo", out)


# =============================================================================
# B) Estrutura gênica: predição vs. real com UTR
# =============================================================================
def fig_gene_structure():
    # Coordenadas reais (URODEC1_LOCUS91967, brizantha subgenome, cromossomo 36b)
    # mRNA: join(7157787..7158312,7158518..7158579,7158740..7158780,7159117..7159709)
    # CDS:  join(7158121..7158312,7158518..7158579,7158740..7158780,7159117..7159346)
    gene_start = 7157787
    exons_mrna = [(7157787, 7158312), (7158518, 7158579), (7158740, 7158780), (7159117, 7159709)]
    cds_start, cds_end = 7158121, 7159346

    def rel(x):
        return x - gene_start

    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.set_xlim(-50, rel(7159709) + 50)
    ax.set_ylim(0, 4.5)
    ax.axis("off")

    thin_h, exon_h = 0.10, 0.55

    def draw_model(y, label, has_utr):
        # linha do intron/base
        ax.plot([0, rel(7159709) if has_utr else rel(cds_end)], [y, y],
                 color=MUTED, linewidth=1.3, zorder=1)
        for (s, e) in exons_mrna:
            if not has_utr:
                # so a parte CDS (predicao original, sem UTR)
                s2, e2 = max(s, cds_start), min(e, cds_end)
                if s2 >= e2:
                    continue
                ax.add_patch(mpatches.FancyBboxPatch(
                    (rel(s2), y - exon_h / 2), rel(e2) - rel(s2), exon_h,
                    boxstyle="round,pad=0,rounding_size=0.0", linewidth=0,
                    facecolor=BLUE, zorder=2))
            else:
                # UTR (parte fora do CDS) em tom claro sequencial; CDS em azul solido
                if s < cds_start:
                    utr_e = min(e, cds_start)
                    ax.add_patch(mpatches.FancyBboxPatch(
                        (rel(s), y - exon_h / 2), rel(utr_e) - rel(s), exon_h,
                        boxstyle="round,pad=0,rounding_size=0.0", linewidth=0.8,
                        edgecolor=BLUE, facecolor="#cde2fb", zorder=2))
                cs, ce = max(s, cds_start), min(e, cds_end)
                if cs < ce:
                    ax.add_patch(mpatches.FancyBboxPatch(
                        (rel(cs), y - exon_h / 2), rel(ce) - rel(cs), exon_h,
                        boxstyle="round,pad=0,rounding_size=0.0", linewidth=0,
                        facecolor=BLUE, zorder=2))
                if e > cds_end:
                    utr_s = max(s, cds_end)
                    ax.add_patch(mpatches.FancyBboxPatch(
                        (rel(utr_s), y - exon_h / 2), rel(e) - rel(utr_s), exon_h,
                        boxstyle="round,pad=0,rounding_size=0.0", linewidth=0.8,
                        edgecolor=BLUE, facecolor="#cde2fb", zorder=2))
        ax.text(-80, y, label, fontsize=10.5, ha="right", va="center", color=INK, fontweight="bold")

    draw_model(3.2, "locus_001\n(predição via\nU. ruziziensis,\nsó miniprot)", has_utr=False)
    draw_model(1.4, "URODEC1_LOCUS91967\n(subgenoma real de\nU. brizantha,\nassembly anotado)", has_utr=True)

    # escala
    scale_y = 0.3
    ax.plot([0, 500], [scale_y, scale_y], color=INK, linewidth=1.5)
    ax.text(250, scale_y + 0.2, "500 pb", fontsize=9, ha="center", color=INK2)

    # legenda
    handles = [
        mpatches.Patch(facecolor=BLUE, edgecolor="none", label="CDS (região codificante)"),
        mpatches.Patch(facecolor="#cde2fb", edgecolor=BLUE, label="UTR (5'/3', não-codificante)"),
        mpatches.Patch(facecolor="none", edgecolor=MUTED, label="íntron"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.32),
              ncol=3, frameon=False, fontsize=9.5)

    ax.set_title(
        "Estrutura gênica: predição original vs. anotação real no subgenoma brizantha",
        fontsize=12, fontweight="bold", color=INK, pad=48,
    )
    plt.tight_layout()
    out = FIG / "esquema_estrutura_genica.png"
    plt.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close()
    print("salvo", out)


# =============================================================================
# C) Convergência de evidências para locus_001
# =============================================================================
def fig_evidence_hub():
    evidences = [
        ("TBLASTN", "bitscore 155,0\n(maior do conjunto)", BLUE),
        ("BLASTN_CDS", "confirmação cruzada\nnucleotídeo × genoma", BLUE),
        ("miniprot", "4 éxons / 3 íntrons\nsem stop prematuro", AQUA),
        ("RBH", "3/3 espécies\n(arroz, milho, Arabidopsis)", AQUA),
        ("Filogenia", "clado FT-like\nbootstrap 99", VIOLET),
        ("Subgenoma\nbrizantha real", "99,4% identidade\n(1 aa de diferença)", GREEN),
    ]
    import numpy as np

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.axis("off")
    ax.set_aspect("equal")

    # centro
    center = Circle((0, 0), 1.35, facecolor=SURFACE, edgecolor=INK, linewidth=2.2, zorder=3)
    ax.add_patch(center)
    ax.text(0, 0.15, "locus_001", fontsize=14, fontweight="bold", ha="center", va="center", zorder=4)
    ax.text(0, -0.25, "candidato\nconfirmado", fontsize=9.5, color=INK2, ha="center", va="center", zorder=4)

    n = len(evidences)
    radius = 4.6
    for i, (title, desc, color) in enumerate(evidences):
        angle = np.pi / 2 - i * (2 * np.pi / n)
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        # linha
        edge_x = 1.35 * np.cos(angle)
        edge_y = 1.35 * np.sin(angle)
        ax.plot([edge_x, x * 0.72], [edge_y, y * 0.72], color=color, linewidth=1.8, zorder=1)
        # marcador check
        ax.scatter([x * 0.72], [y * 0.72], s=110, color=GOOD, zorder=2, marker="o")
        ax.text(x * 0.72, y * 0.72, "✓", fontsize=10, color="white", ha="center", va="center",
                 zorder=3, fontweight="bold")
        # caixa de texto
        box = FancyBboxPatch(
            (x - 1.5, y - 0.55), 3.0, 1.1,
            boxstyle="round,pad=0.05,rounding_size=0.12",
            linewidth=1.4, edgecolor=color, facecolor=SURFACE, zorder=2,
        )
        ax.add_patch(box)
        ax.text(x, y + 0.18, title, fontsize=10.5, fontweight="bold", ha="center", va="center", color=INK)
        ax.text(x, y - 0.22, desc, fontsize=8.5, ha="center", va="center", color=INK2, linespacing=1.3)

    plt.tight_layout()
    out = FIG / "esquema_convergencia_evidencias.png"
    plt.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close()
    print("salvo", out)


# =============================================================================
# D) Mapa de subgenomas do assembly de U. decumbens (36 cromossomos, b vs rd)
# =============================================================================
def fig_subgenome_map():
    chrom_labels = [
        "1b", "2b", "3rd", "4rd", "5rd", "6rd", "7b", "8b", "9rd", "10rd",
        "11b", "12b", "13rd", "14rd", "15b", "16b", "17b", "18b", "19rd", "20rd",
        "21rd", "22rd", "23rd", "24b", "25rd", "26rd", "27b", "28b", "29rd", "30rd",
        "31b", "32b", "33rd", "34rd", "35b", "36b",
    ]
    highlight = {"33rd": "cópia rd\n(100% id.)", "34rd": "cópia rd\n(99,4% id.)", "36b": "cópia real\nbrizantha\n(99,4% id.)"}

    fig, ax = plt.subplots(figsize=(13, 5.5))
    ncols = 9
    nrows = 4
    ax.set_xlim(0, ncols)
    ax.set_ylim(0, nrows + 1.6)
    ax.axis("off")

    for i, label in enumerate(chrom_labels):
        col = i % ncols
        row = nrows - 1 - (i // ncols)
        is_b = label.endswith("b")
        color = BLUE if is_b else AQUA
        is_hit = label in highlight
        box = FancyBboxPatch(
            (col + 0.08, row + 0.08), 0.84, 0.84,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=2.6 if is_hit else 0.8,
            edgecolor=RED if is_hit else "white",
            facecolor=color, zorder=2,
        )
        ax.add_patch(box)
        ax.text(col + 0.5, row + 0.5, label, fontsize=10, color="white", fontweight="bold",
                 ha="center", va="center", zorder=3)
        if is_hit:
            ax.annotate(
                highlight[label], xy=(col + 0.5, row + 0.08),
                xytext=(col + 0.5, row - 0.55),
                ha="center", va="top", fontsize=7.8, color=INK,
                arrowprops=dict(arrowstyle="-", color=RED, linewidth=1.2),
            )

    handles = [
        mpatches.Patch(color=BLUE, label="Subgenoma U. brizantha diploide (\"b\", 18 cromossomos)"),
        mpatches.Patch(color=AQUA, label="Subgenoma U. ruziziensis/decumbens diploide (\"rd\", 18 cromossomos)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.18),
              ncol=2, frameon=False, fontsize=10.5)
    ax.set_title(
        "Composição de subgenomas do assembly U. decumbens cv. Basilisk (GCA_964030465.3)\n"
        "e posição das cópias homeólogas de locus_001",
        fontsize=11.5, fontweight="bold", color=INK, pad=44,
    )

    plt.tight_layout()
    out = FIG / "esquema_subgenomas_decumbens.png"
    plt.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close()
    print("salvo", out)


if __name__ == "__main__":
    fig_pipeline()
    fig_gene_structure()
    fig_evidence_hub()
    fig_subgenome_map()
