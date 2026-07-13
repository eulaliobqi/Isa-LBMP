# Correlação locus_XXX (pipeline) → gene ID real no subgenoma *U. brizantha*

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

1. Baixado `protein.faa` (proteoma anotado, 150.710 proteínas,
   ~126 mil genes preditos pelo pipeline REAT) e `genomic.gff` via
   NCBI Datasets API (`GCA_964030465.3`), mais o `sequence_reports`
   (mapeamento accession→nome de cromossomo, ex. `OZ075146.1`→`36b`).
2. BLASTP local (BLAST+ 2.17.0, binário Windows baixado do NCBI) de
   cada uma das **6 proteínas candidatas confirmadas** do pipeline
   principal (`locus_001`, `002`, `005`, `009`, `012`, `015` — as únicas
   com `pebp_confirmed=True` e `rbh_confirmed=True` em
   `results/final_candidates_table.tsv`) contra o proteoma de
   *U. decumbens*.
3. Hits cruzados com `genomic.gff` (`protein_id=`/`locus_tag=`) e o
   mapeamento de cromossomo pra identificar em qual subgenoma
   (sufixo "b" = *brizantha* real, "rd" = *ruziziensis*/*decumbens*)
   cada hit está.

## Resultado — tabela de correlação (2026-07-13)

Tabela completa em `loci_brizantha_correlation.tsv`. Resumo:

| locus (pipeline) | classificação | melhor cópia subgenoma **b** | identidade | gene ID (locus_tag) |
|---|---|---|---|---|
| `locus_001` | FT-like (Hd3a/RFT1-type) | CAL5051155.1 | **99,4%** | URODEC1_LOCUS91967 (chr 36b) |
| `locus_002` | FT-like (fragmento fraco, 76 aa) | CAL5082316.1 / CAL5087382.1 (empate) | 100,0% | URODEC1_LOCUS109225 / URODEC1_LOCUS112169 (chr 7b/8b) |
| `locus_005` | FT-like (ZCN8-type) | CAL4902238.1 | 91,6% | URODEC1_LOCUS9836 (chr 11b) |
| `locus_009` | FT-like (par 009/015) | CAL4952393.1 | 98,8% | URODEC1_LOCUS39484 (chr 17b) |
| `locus_012` | **MFT-like** (reclassificado, não é FT-like) | CAL4950975.1 | 100,0% | URODEC1_LOCUS38694 (chr 17b) |
| `locus_015` | FT-like (par 009/015) | CAL4952393.1 | 98,8% | URODEC1_LOCUS39484 (chr 17b) |

**Achados notáveis**:
- `locus_001` é o único caso em que a cópia do subgenoma **b** já é, ao mesmo tempo, o melhor hit *geral* do proteoma (empatado ou superior às cópias rd) — nos outros 5 loci, o melhor hit geral cai numa cópia **rd**, e a cópia **b** correspondente aparece com identidade igual ou (no caso de `locus_005`) sensivelmente menor (91,6% vs. 100%).
- `locus_009` e `locus_015` (parálogos quase idênticos já identificados na filogenia, branch length ≈0) convergem exatamente para a mesma cópia b (`CAL4952393.1`/`URODEC1_LOCUS39484`) — replicação independente do sinal de parologia.
- `locus_012` está na tabela só por rastreabilidade — já foi reclassificado como **MFT-like**, não é candidato ao gene FT-like da proposta (ver `artigo.md` §3.4/§4).
- `locus_002` é o caso mais frágil: fragmento curto (76 aa, só evidência TBLASTN) com 3 cópias empatadas em 100% no proteoma decumbens (1 rd + 2 b) — não dá pra apontar uma única cópia b "correta" com esse nível de evidência.

**Detalhe do achado original de `locus_001`**: `CAL5051155.1` (gene
`URODEC1_LOCUS91967`, cromossomo `36b`, OZ075146.1:7.157.787–7.159.709,
fita +) difere da predição de `locus_001` (via *U. ruziziensis*) por
**exatamente 1 aminoácido** (posição ~53: V→I). Também ganha **anotação
completa com UTRs reais** — mRNA de 1.222 nt vs. CDS de 525 nt (174 aa),
~697 nt de UTR modelada pelo pipeline de anotação REAT do assembly
(não confirmado se ancorada em RNA-seq desta espécie especificamente ou
em modelo). Estrutura gênica: 4 éxons/3 íntrons, mesma contagem predita
para `locus_001` em *U. ruziziensis* — validação cruzada da estrutura.

## Arquivos neste diretório

- `locus_{001,002,005,009,012,015}_query.fasta` — proteína de cada locus usada como query
- `locus_{001,002,005,009,012,015}_vs_decumbens.tsv` — resultado bruto do BLASTP (outfmt 6, top 10-20 hits)
- `loci_brizantha_correlation.tsv` — tabela consolidada de correlação (gerada 2026-07-13)
- `brizantha_CDS.fasta` — CDS (525 bp) do gene `URODEC1_LOCUS91967` (cópia *brizantha* de locus_001)
- `brizantha_mRNA_with_UTR.fasta` — mRNA completo com UTRs (1.222 nt, locus_001)

Arquivos grandes de trabalho (proteoma completo, GFF completo ~880 MB,
banco BLAST, binários) não foram mantidos neste diretório após a
extração — refazer via NCBI Datasets API se precisar rodar de novo
(accession `GCA_964030465.3`).
