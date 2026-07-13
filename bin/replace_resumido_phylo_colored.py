#!/usr/bin/env python3
"""Substitui apenas a arvore filogenetica em relatorio_isa_lbmp-resumido.docx
pela versao colorida por clado + legenda (filogenia-certa-colorida-legenda.png).

Mesma tecnica do replace_resumido_figures.py: edita o zip do docx diretamente,
recalculando o wp:extent/a:ext do bloco rId9 para a nova proporcao da imagem
(evita esticar/distorcer).
"""
import re
import shutil
import zipfile
from pathlib import Path
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
DOCX = ROOT / "relatorio_isa_lbmp-resumido.docx"
FIG = ROOT / "results" / "results" / "results" / "figures"

NEW_PHYLO = FIG / "filogenia-certa-colorida-legenda.png"
TARGET_MEDIA_PHYLO = "word/media/image4.png"
PHYLO_RID = "rId9"


def get_drawing_block(xml, rid):
    idx = xml.find(f'r:embed="{rid}"')
    if idx == -1:
        raise ValueError(f"rid {rid} nao encontrado no document.xml")
    start = xml.rfind("<w:drawing>", 0, idx)
    end = xml.find("</w:drawing>", idx) + len("</w:drawing>")
    return start, end, xml[start:end]


def recompute_extent(block, new_w, new_h):
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
    if not NEW_PHYLO.exists():
        raise FileNotFoundError(NEW_PHYLO)

    backup = DOCX.with_suffix(".docx.bak2")
    shutil.copy2(DOCX, backup)
    print(f"Backup salvo em {backup}")

    with zipfile.ZipFile(DOCX, "r") as zin:
        names = zin.namelist()
        data = {name: zin.read(name) for name in names}

    document_xml = data["word/document.xml"].decode("utf-8")

    with Image.open(NEW_PHYLO) as im:
        new_w, new_h = im.size
    start, end, block = get_drawing_block(document_xml, PHYLO_RID)
    new_block, old_cx, new_cy = recompute_extent(block, new_w, new_h)
    document_xml = document_xml[:start] + new_block + document_xml[end:]
    print(f"Extent atualizado: cx={old_cx} cy={new_cy}")

    data["word/document.xml"] = document_xml.encode("utf-8")
    data[TARGET_MEDIA_PHYLO] = NEW_PHYLO.read_bytes()

    tmp_path = DOCX.with_suffix(".docx.tmp")
    with zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_DEFLATED) as zout:
        for name in names:
            zout.writestr(name, data[name])

    tmp_path.replace(DOCX)
    print(f"[OK] {DOCX} atualizado.")
    print(f"  - {TARGET_MEDIA_PHYLO} <- {NEW_PHYLO.name}")


if __name__ == "__main__":
    main()
