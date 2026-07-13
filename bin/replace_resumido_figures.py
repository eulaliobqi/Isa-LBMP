#!/usr/bin/env python3
"""Substitui a arvore filogenetica e o heatmap de bitscore embutidos em
relatorio_isa_lbmp-resumido.docx pelas versoes corrigidas, sem regenerar
o documento inteiro (preserva todo o resto do docx).

image3.png (dentro do docx) = bitscore heatmap antigo (6 colunas) -> figures/bitscore_heatmap.png (4 colunas, corrigido)
image4.png (dentro do docx) = arvore filogenetica antiga -> figures/arvore-filo-correta.png (nova, legivel)

A arvore nova tem aspect ratio diferente da antiga, entao o wp:extent/a:ext
do bloco XML referente a rId9 (image4) e recalculado para manter a largura
original e evitar distorcao.
"""
import re
import shutil
import zipfile
from pathlib import Path
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
DOCX = ROOT / "relatorio_isa_lbmp-resumido.docx"
FIG = ROOT / "results" / "results" / "results" / "figures"

NEW_BITSCORE = FIG / "bitscore_heatmap.png"
NEW_PHYLO = FIG / "arvore-filo-correta.png"

# mapeamento fixo, confirmado por inspecao manual do docx extraido
TARGET_MEDIA_BITSCORE = "word/media/image3.png"
TARGET_MEDIA_PHYLO = "word/media/image4.png"
PHYLO_RID = "rId9"


def get_drawing_block(xml, rid):
    """Extrai o bloco <w:drawing>...</w:drawing> que contem o r:embed=rid."""
    idx = xml.find(f'r:embed="{rid}"')
    if idx == -1:
        raise ValueError(f"rid {rid} nao encontrado no document.xml")
    start = xml.rfind("<w:drawing>", 0, idx)
    end = xml.find("</w:drawing>", idx) + len("</w:drawing>")
    return start, end, xml[start:end]


def recompute_extent(block, new_w, new_h):
    """Atualiza todos os cx/cy dentro do bloco (wp:extent e a:ext) preservando
    a largura original e recalculando a altura pela nova proporcao."""
    m = re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"', block)
    if not m:
        raise ValueError("wp:extent nao encontrado no bloco")
    old_cx = int(m.group(1))
    new_cy = round(old_cx * new_h / new_w)

    old_cx_str, old_cy_str = m.group(1), m.group(2)
    block = block.replace(f'cx="{old_cx_str}" cy="{old_cy_str}"', f'cx="{old_cx_str}" cy="{new_cy}"')
    return block, old_cx, new_cy


def main():
    if not DOCX.exists():
        raise FileNotFoundError(DOCX)
    if not NEW_BITSCORE.exists() or not NEW_PHYLO.exists():
        raise FileNotFoundError("figuras novas nao encontradas")

    backup = DOCX.with_suffix(".docx.bak")
    shutil.copy2(DOCX, backup)
    print(f"Backup salvo em {backup}")

    with zipfile.ZipFile(DOCX, "r") as zin:
        names = zin.namelist()
        data = {name: zin.read(name) for name in names}

    document_xml = data["word/document.xml"].decode("utf-8")

    # Recalcula o extent do bloco da arvore filogenetica (aspect ratio mudou)
    with Image.open(NEW_PHYLO) as im:
        new_w, new_h = im.size
    start, end, block = get_drawing_block(document_xml, PHYLO_RID)
    new_block, old_cx, new_cy = recompute_extent(block, new_w, new_h)
    document_xml = document_xml[:start] + new_block + document_xml[end:]
    print(f"Extent da arvore atualizado: cx={old_cx} cy={new_cy} (era proporcional a imagem antiga)")

    data["word/document.xml"] = document_xml.encode("utf-8")

    # Substitui os bytes das imagens (bitscore mantem mesmas dimensoes, sem precisar mexer no XML)
    data[TARGET_MEDIA_BITSCORE] = NEW_BITSCORE.read_bytes()
    data[TARGET_MEDIA_PHYLO] = NEW_PHYLO.read_bytes()

    tmp_path = DOCX.with_suffix(".docx.tmp")
    with zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_DEFLATED) as zout:
        for name in names:
            zout.writestr(name, data[name])

    tmp_path.replace(DOCX)
    print(f"[OK] {DOCX} atualizado com sucesso.")
    print(f"  - {TARGET_MEDIA_BITSCORE} <- {NEW_BITSCORE.name}")
    print(f"  - {TARGET_MEDIA_PHYLO} <- {NEW_PHYLO.name}")


if __name__ == "__main__":
    main()
