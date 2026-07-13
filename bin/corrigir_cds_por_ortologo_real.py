# -*- coding: utf-8 -*-
"""Corrige artefatos de N-terminal nas sequencias de locus_001/005/009, detectados
por alinhamento par-a-par contra os ortologos reais confirmados no subgenoma
brizantha do assembly de U. decumbens (brizantha_search/real_orthologs_decumbens.fasta).

Achado: o miniprot, ao forcar o alinhamento completo (1..N) contra a proteina de
referencia usada como query (Hd3a/ZCN8/Arabidopsis-FT), por vezes estende a predicao
de CDS alguns codons para dentro da UTR 5' quando a extremidade N-terminal diverge da
referencia -- produzindo residuos extras no inicio da proteina predita que NAO
correspondem a nada no gene real. So detectavel comparando contra um ortologo
independente com anotacao real (nao predita).

Cada correcao abaixo foi verificada por:
1. Alinhamento par-a-par (Bio.Align.PairwiseAligner) predito x real, confirmando
   o deslocamento exato.
2. Confirmacao de que o CDS truncado no ponto do deslocamento comeca com ATG.
3. Traducao do CDS truncado sem stop interno, terminando em stop codon.
4. Traducao == proteina esperada pelo alinhamento.

Nao se aplica a locus_002/012 (sem modelo de exon, so fragmento de TBLASTN -- nao ha
CDS para truncar) nem a locus_015 (nao tem residuo extra; ao contrario, FALTA
sequencia real no inicio -- nada a corrigir, so documentar).

Roda de dentro do repo:
    python bin/corrigir_cds_por_ortologo_real.py
"""
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
SEQ_DIR = ROOT / "candidate_sequences"
NT_FASTA = SEQ_DIR / "nucleotide_sequences.fasta"
PROT_FASTA = SEQ_DIR / "protein_sequences.fasta"

# Deslocamento verificado (nt) para cada locus com artefato confirmado.
# Fonte: alinhamento par-a-par contra o ortologo real no subgenoma brizantha.
TRIMS_NT = {
    "locus_001": 15,   # 5 aa extras (LVAIG) vs CAL5051155.1
    "locus_005": 24,   # 8 aa extras (YNSTHPLV) vs CAL4902238.1
    "locus_009": 42,   # 14 aa extras (EDPSTRHRSASYRR) vs CAL4952393.1
}


def load_fasta_dict(path):
    return {r.id.split()[0].split("|")[0]: str(r.seq) for r in SeqIO.parse(path, "fasta")}


def main():
    nt = load_fasta_dict(NT_FASTA)
    prot = load_fasta_dict(PROT_FASTA)

    nt_out_order = ["locus_001", "locus_002", "locus_005", "locus_009", "locus_012", "locus_015"]

    orig_headers = {r.id.split()[0]: r.description for r in SeqIO.parse(NT_FASTA, "fasta")}

    log = []
    nt_records = []
    prot_records = []

    for locus in nt_out_order:
        full_nt = nt[locus]
        full_prot = prot[locus]

        if locus in TRIMS_NT:
            n = TRIMS_NT[locus]
            trimmed_nt = full_nt[n:]
            translated = str(Seq(trimmed_nt).translate(to_stop=False))
            assert trimmed_nt[:3] == "ATG", f"{locus}: CDS corrigido nao comeca com ATG"
            assert "*" not in translated[:-1], f"{locus}: stop interno apos correcao"
            assert translated.endswith("*"), f"{locus}: CDS corrigido nao termina em stop"
            corrected_prot = translated.rstrip("*")

            log.append(
                f"{locus}: CORRIGIDO. Predicao bruta do miniprot tinha {n // 3} aminoacido(s) "
                f"extra(s) no N-terminal ({full_prot[:n // 3]}) que nao correspondem ao gene "
                f"real (confirmado por alinhamento contra ortologo no subgenoma brizantha, "
                f"brizantha_search/real_orthologs_decumbens.fasta). CDS bruto {len(full_nt)}bp "
                f"-> corrigido {len(trimmed_nt)}bp. Proteina bruta {len(full_prot)}aa -> "
                f"corrigida {len(corrected_prot)}aa."
            )
            nt_records.append((
                locus,
                f"{locus} CDS_reconstructed_CORRIGIDO miniprot_model=4_exons "
                f"N-terminal_verificado_contra_ortologo_real removidos={n // 3}aa_espurios",
                trimmed_nt,
            ))
            nt_records.append((
                f"{locus}_bruto_miniprot",
                f"{locus}_bruto_miniprot CDS_reconstructed_SEM_CORRECAO (predicao original do "
                f"pipeline, mantida para rastreabilidade -- contem {n // 3} aa espurios no N-terminal)",
                full_nt,
            ))
            prot_records.append((locus, f"{locus} (proteina, CORRIGIDA)", corrected_prot))
            prot_records.append((
                f"{locus}_bruto_miniprot",
                f"{locus}_bruto_miniprot (proteina, predicao original SEM correcao)",
                full_prot,
            ))
        elif locus == "locus_015":
            log.append(
                f"{locus}: SEM CORRECAO POSSIVEL (nao ha residuo espurio a remover). Alinhamento "
                f"contra o ortologo real (CAL4952393.1, mesma copia de locus_009) mostra que "
                f"locus_015 esta FALTANDO os 4 primeiros aminoacidos do gene real (MRRG) -- "
                f"o miniprot simplesmente nao encontrou essa extensao N-terminal durante o "
                f"alinhamento spliced contra a proteina de Arabidopsis FT (Target=NP_001320342.1 "
                f"51 219, ou seja, comeca na posicao 51 de 219 da referencia). Predicao mantida "
                f"como esta -- e parcial mas nao contem erro, so incompletude."
            )
            nt_records.append((locus, f"{locus} CDS_reconstructed miniprot_model=4_exons "
                                       f"INCOMPLETO_no_N-terminal_ver_README", full_nt))
            prot_records.append((locus, f"{locus} (proteina, incompleta no N-terminal, ver README)", full_prot))
        else:  # locus_002, locus_012 -- sem modelo de exon
            log.append(
                f"{locus}: SEM MODELO DE EXON -- nao ha CDS reconstruido para corrigir (so "
                f"fragmento de proteina via TBLASTN). Mantido como span genomico bruto, ja "
                f"rotulado RAW_GENOMIC_SPAN. Contexto novo desta verificacao: alinhamento contra "
                f"o ortologo real mostra que o fragmento corresponde a uma janela interna da "
                f"proteina completa (ver README para as coordenadas exatas)."
            )
            nt_records.append((locus, orig_headers[locus], full_nt))
            prot_records.append((locus, f"{locus} (proteina, fragmento)", full_prot))

    with open(NT_FASTA, "w", encoding="utf-8") as out:
        for locus, header, seq in nt_records:
            out.write(f">{header}\n{seq}\n")

    with open(PROT_FASTA, "w", encoding="utf-8") as out:
        for locus, header, seq in prot_records:
            out.write(f">{header}\n{seq}\n")

    with open(SEQ_DIR / "validation_log.txt", "a", encoding="utf-8") as out:
        out.write("\n" + "=" * 70 + "\n")
        out.write("Verificacao adicional 2026-07-13: correcao de artefato de N-terminal\n")
        out.write("(alinhamento par-a-par contra ortologos reais no subgenoma brizantha)\n")
        out.write("=" * 70 + "\n")
        out.write("\n".join(log) + "\n")

    print("Concluido. Log:")
    print("\n".join(log))


if __name__ == "__main__":
    main()
