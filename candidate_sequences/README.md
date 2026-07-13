# Sequências candidatas — nucleotídeo e proteína (2026-07-13)

Sequências completas das 6 loci confirmadas (`pebp_confirmed=True` e `rbh_confirmed=True`
em `results/final_candidates_table.tsv`), geradas por `bin/extract_locus_nucleotide_sequences.py`.

## Fonte do genoma

Genoma alvo: *Urochloa ruziziensis* assembly **Embrapa_Uruz_1.0** (`GCA_015476505.1`),
baixado via NCBI Datasets API v2 (`GENOME_FASTA`). O arquivo (`genoma.fa`, ~600 MB) **não é
mantido no repositório** (está no `.gitignore`, mesmo padrão já usado para `brizantha_search/`)
— refazer o download se precisar rodar de novo.

## `protein_sequences.fasta`

Cópia direta de `results/results/results/candidate_proteins_confirmed.fasta` (a fonte fica em
`results/`, que é git-ignored — esta cópia torna o arquivo rastreável no repo). Nenhuma
sequência foi alterada.

## `nucleotide_sequences.fasta`

Três categorias, diferenciadas pelo header:

1. **CDS reconstruído e corrigido** (`locus_001`, `005`, `009`) — essas 3 loci têm modelo gênico
   do miniprot (spliced-aware, com éxons/íntrons definidos). O CDS foi reconstruído concatenando
   todos os éxons na ordem genômica (reverse-complement se fita `-`), validado por tradução, e
   **corrigido** removendo um artefato de N-terminal descoberto em 2026-07-13 (ver seção
   "Achado: artefato de N-terminal" abaixo). Header: `CDS_reconstructed_CORRIGIDO`. A predição
   original sem correção é mantida logo em seguida, com header
   `<locus>_bruto_miniprot CDS_reconstructed_SEM_CORRECAO`, para rastreabilidade completa.

2. **CDS reconstruído, incompleto no N-terminal** (`locus_015`) — tem modelo gênico do miniprot
   e foi validado normalmente, mas a comparação contra o ortólogo real mostrou que a predição
   **não tem** resíduos espúrios — ao contrário, **falta** um pedaço real do início do gene
   (ver abaixo). Não há correção possível sem reprocessar o alinhamento; mantido como está,
   com aviso no header (`INCOMPLETO_no_N-terminal_ver_README`).

3. **Span genômico bruto** (`locus_002`, `012`) — essas 2 loci têm evidência só de TBLASTN
   (proteína fragmentada, sem modelo de éxon/íntron do miniprot). O header diz explicitamente
   `RAW_GENOMIC_SPAN exon_structure=UNKNOWN pode_incluir_introns_ou_UTR`. **Isso não é um CDS
   limpo** — é o trecho genômico inteiro entre as coordenadas de início e fim do locus
   (`locus_002`: 4.674 pb; `locus_012`: 20.344 pb), que quase certamente inclui íntron(s) e/ou
   UTR não removidos. Não há como reconstruir um CDS confiável para essas 2 loci sem uma
   predição gênica adicional (miniprot não gerou modelo pra elas na rodagem original).

**Nota técnica sobre a reconstrução original**: concatena todos os CDS antes de traduzir, nunca
éxon a éxon — miniprot tolera pequenos frameshifts dentro de um éxon isolado, e traduzir exon a
exon (como faria `gffread -y`) produziria proteína fora de frame. Essa foi uma lição real do
pipeline (ver `artigo.md`), e o script de extração (`bin/extract_locus_nucleotide_sequences.py`)
reaplica a mesma lógica correta de `bin/translate_miniprot_gff.py`.

## Achado: artefato de N-terminal na predição do miniprot (2026-07-13)

**Pergunta que originou a investigação**: por que a proteína predita de `locus_001` não começa
com Met (M)/ATG, se M é sempre o primeiro resíduo de uma proteína real?

**Descoberta**: comparando `locus_001` (proteína predita, 179 aa) com o ortólogo real e
independentemente anotado no subgenoma *brizantha* (`CAL5051155.1`/`URODEC1_LOCUS91967`, já
confirmado em `brizantha_search/`), a proteína real bate exatamente com `locus_001` a partir do
6º resíduo — os 5 primeiros resíduos preditos (`LVAIG`) não existem no gene real.

**Causa provável**: o miniprot alinha a proteína inteira usada como referência (ex.: Hd3a, 179
aa, posições 1–179) contra o genoma. Quando a extremidade N-terminal do gene real diverge da
referência, o algoritmo pode estender o alinhamento — e portanto a predição de CDS — alguns
códons para dentro da UTR 5' (região não-codificante) para tentar completar o alinhamento
inteiro, "traduzindo" um pedaço de UTR como se fosse região codificante real.

**Extensão da verificação a todas as 6 loci** (a pedido do usuário): baixei as proteínas reais
dos ortólogos já identificados no subgenoma *brizantha* para as outras 5 loci
(`brizantha_search/real_orthologs_decumbens.fasta`, via NCBI Datasets API,
`GCA_964030465.3`) e alinhei cada predição contra seu ortólogo real
(`Bio.Align.PairwiseAligner`, modo local). Resultado:

| Locus | Situação | Ação |
|---|---|---|
| `locus_001` | 5 aa espúrios no N-terminal (`LVAIG`), confirmado vs. `CAL5051155.1` | **Corrigido**: CDS 540→525 pb, proteína 179→174 aa |
| `locus_005` | 8 aa espúrios no N-terminal (`YNSTHPLV`), confirmado vs. `CAL4902238.1` | **Corrigido**: CDS 528→504 pb, proteína 175→167 aa |
| `locus_009` | 14 aa espúrios no N-terminal (`EDPSTRHRSASYRR`), confirmado vs. `CAL4952393.1` | **Corrigido**: CDS 564→522 pb, proteína 187→173 aa |
| `locus_015` | Nenhum resíduo espúrio — ao contrário, **faltam** 4 aa reais no início (`MRRG`) vs. `CAL4952393.1`; miniprot alinhou só a partir da posição 51 de 219 da referência (*Arabidopsis* FT) | Sem correção possível; mantido incompleto, documentado |
| `locus_002` | Fragmento (76 aa) corresponde às posições 101–173 de uma proteína real de 178 aa (`CAL5082316.1`/`CAL5087382.1`), com 3 aa iniciais (`TFS`) também não batendo exatamente | Sem modelo de éxon — não há CDS pra corrigir; contexto documentado |
| `locus_012` | Fragmento (72 aa) corresponde às posições 100–171 de uma proteína real de 171 aa (`CAL4950975.1`), sem resíduo espúrio | Sem modelo de éxon — não há CDS pra corrigir; contexto documentado |

**Cada correção foi verificada, não só aplicada**: para as 3 loci corrigidas, o script
`bin/corrigir_cds_por_ortologo_real.py` confirma automaticamente que (1) o CDS truncado no ponto
do deslocamento começa com `ATG`, (2) a tradução não tem stop códon interno, (3) termina em stop
códon, e (4) a proteína resultante bate exatamente com o esperado pelo alinhamento — se qualquer
checagem falhar, o script aborta (`assert`) sem escrever o arquivo.

**Implicação prática**: os números "179 aa" (`locus_001`), "175 aa" (`locus_005`) e "187 aa"
(`locus_009`) que aparecem nas tabelas de resultado do pipeline original (`final_candidates_table.tsv`,
`artigo.md`, `relatorio_isa_lbmp.docx`) refletem a **predição bruta do miniprot**, correta como
registro do que o pipeline de fato produziu — mas as sequências corrigidas aqui em
`candidate_sequences/` (174, 167 e 173 aa) são as que devem ser usadas para qualquer desenho de
primer/dsRNA, por representarem o gene real com mais precisão.

## `validation_log.txt`

PASS/FAIL da extração original (as 4 loci com modelo gênico) + tamanhos informativos das 2 sem
modelo, seguido do log de correção do N-terminal (seção acima) apensado em 2026-07-13.
Se qualquer uma das 4 loci com modelo tivesse falhado na validação original, o script abortaria
sem escrever nenhum arquivo — nenhuma sequência não-validada chega a este diretório por engano.
