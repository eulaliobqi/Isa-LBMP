# Busca de locus_001 no subgenoma real de *U. brizantha*

Gerado 2026-07-11. Não existe assembly próprio de *U. brizantha* (segue "órfã genômica"),
mas o assembly haplótipo-resolvido de *Urochloa decumbens* cv. Basilisk
(**GCA_964030465.3**, de Vega Group/Earlham Institute, depositado no ENA,
BioProject `PRJEB73762`, paper: de Jesus et al. via bioRxiv
10.1101/2024.09.25.614935, PMC12005165) é um alotetraploide cujos
cromossomos foram rotulados por origem de subgenoma após análise de
cobertura de reads: sufixo **"b" = subgenoma derivado de *U. brizantha*
diploide**, sufixo **"rd" = subgenoma derivado de *U. ruziziensis*/
*U. decumbens* diploide** (confirmado no texto do artigo, não assumido).

## Metodologia

1. Baixado `protein.faa.gz` (proteoma anotado, 150.710 proteínas,
   ~126 mil genes preditos pelo pipeline REAT) e `feature_table.txt.gz`
   (mapeamento proteína→cromossomo) via FTP do NCBI.
2. BLASTP local (BLAST+ 2.17.0, binário Windows baixado do NCBI) da
   proteína de `locus_001` (predita em *U. ruziziensis*, script
   `bin/gerar_relatorio_docx.py`/pipeline principal) contra o proteoma
   de *U. decumbens*.
3. Top hits cruzados contra `feature_table.txt` pra identificar em qual
   subgenoma (sufixo "b" vs "rd") cada hit está.

## Resultado

5 cópias quase idênticas do gene (99,4–100% identidade), refletindo os
múltiplos haplótipos do assembly:

| Accession | Identidade | Cromossomo | Subgenoma |
|---|---|---|---|
| CAL5036879.1 | 100,0% | 33rd | ruziziensis/decumbens |
| CAM0151951.1 | 99,4% | scaffold não localizado | — |
| CAM0146821.1 | 99,4% | scaffold não localizado | — |
| **CAL5051155.1** | **99,4%** | **36b** | ***brizantha*** |
| CAL5044486.1 | 99,4% | 34rd | ruziziensis/decumbens |

**`CAL5051155.1` (gene `URODEC1_LOCUS91967`, cromossomo `36b`,
OZ075146.1:7.157.787–7.159.709, fita +) é a cópia real do subgenoma
*brizantha*** — a que interessa para o projeto. Difere da predição de
`locus_001` (via *U. ruziziensis*) por **exatamente 1 aminoácido**
(posição ~53: V→I), confirmando altíssima conservação e validando a
identificação original.

**Ganho adicional importante**: diferente da predição em *U. ruziziensis*
(só CDS, via miniprot, sem UTR), esse gene tem **anotação completa com
UTRs reais** — mRNA de 1.222 nt vs. CDS de 525 nt (174 aa), ou seja,
~697 nt de UTR modelada a partir de evidência real do pipeline de
anotação REAT (provavelmente ancorada em RNA-seq/evidência transcricional
usada pelos autores do assembly, não confirmado se foi confirmação direta
por RNA-seq desta espécie ou modelo). Estrutura gênica: 4 éxons/3 íntrons,
mesma contagem predita para `locus_001` em *U. ruziziensis* — validação
cruzada da estrutura.

## Arquivos neste diretório

- `locus_001_query.fasta` — proteína de locus_001 usada como query
- `locus_001_vs_decumbens.tsv` — resultado bruto do BLASTP (top 20 hits)
- `brizantha_CDS.fasta` — CDS (525 bp) do gene `URODEC1_LOCUS91967` (cópia *brizantha*)
- `brizantha_mRNA_with_UTR.fasta` — mRNA completo com UTRs (1.222 nt)

Arquivos grandes de trabalho (proteoma completo, banco BLAST, binários) não
foram mantidos neste diretório após a extração — refazer via `datasets`/FTP
se precisar rodar de novo (accession `GCA_964030465.3`).
