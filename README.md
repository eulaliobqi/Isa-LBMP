# Isa-LBMP — Homólogos FT-like (PEBP) em *Urochloa*

Pipeline Nextflow DSL2 para identificar genes homólogos a FLOWERING LOCUS T
(FT, família PEBP) no genoma de *Urochloa ruziziensis* (assembly
Embrapa_Uruz_1.0, proxy diploide para *U. brizantha*, que não tem assembly
cromossômico próprio), usando como referência proteínas ortólogas de Hd3a e
RFT1 (*Oryza sativa*), ZCN8 (*Zea mays*) e FT (*Arabidopsis thaliana*).

Contexto científico completo e justificativa de cada decisão metodológica:
`C:\Users\eulal\.claude\plans\objetivo-principal-usar-o-synthetic-peach.md`.

Este é o passo de identificação/validação de homólogos dentro do projeto
maior de silenciamento gênico (RNAi) de SBP5 + FT em *U. brizantha*, descrito
em `SBP5 + FT.docx`.

## Pré-requisitos

- Servidor `eulalio@200.235.143.10` (Debian, 32 cores) — o pipeline é
  100% CPU-bound (BLAST/HMMER/MAFFT/IQ-TREE não usam GPU).
- Acesso SSH configurado (chave própria, sem senha) — ver Passo 0 do plano
  científico.
- Ambiente `isa-lbmp` criado via mamba:
  ```bash
  mamba env create -f environment.yml
  ```
- **Antes de transferir o genoma (583 MB), confira espaço em disco:**
  ```bash
  df -h ~
  ```
  A partição `/home` do servidor já chegou a 96% de uso em outro projeto do
  mesmo usuário — confirme espaço livre suficiente antes do `scp`.

## Como transferir o genoma para o servidor

O genoma **não é versionado** neste repositório (`.gitignore` exclui os
FASTA grandes). Transferir separadamente via `scp` a partir do notebook
Windows:

```bash
ssh isalbmp-server "mkdir -p /home/eulalio/databases/urochloa_ruziziensis"
scp -C "sequencias cromossomos brizantha unidas.txt" \
    isalbmp-server:/home/eulalio/databases/urochloa_ruziziensis/urochloa_ruziziensis_Embrapa_Uruz_1.0.fasta
```

## Como rodar

Sempre dentro de `screen` — sem isso o processo morre por SIGTTOU ao
desconectar o SSH:

```bash
ssh isalbmp-server
screen -S isa-lbmp
mamba activate isa-lbmp   # ou deixe o profile -profile mamba cuidar disso
cd ~/isa-lbmp
git pull                  # sempre sincronizar com o último commit antes de rodar

nextflow run main.nf \
    -profile mamba,local \
    -params-file params.yaml \
    --genome /home/eulalio/databases/urochloa_ruziziensis/urochloa_ruziziensis_Embrapa_Uruz_1.0.fasta
```

`Ctrl+A, D` para desanexar da screen com segurança; `screen -r isa-lbmp` para
reconectar depois.

## Processos do pipeline

| # | Processo | O que faz | Output principal |
|---|---|---|---|
| 1 | `BUILD_BLASTDB` | `makeblastdb` no genoma Urochloa | `urochloa_db.*` |
| 2 | `TBLASTN` | Triagem ampla proteína×genoma (`-evalue 1e-5`) | `tblastn_raw.tsv` |
| 3 | `MINIPROT_INDEX` + `MINIPROT_ALIGN` | Estrutura gênica spliced-aware (introns) + extração da proteína predita via `gffread` | `miniprot_hits.gff3`, `miniprot_predicted_proteins.fasta` |
| 4 | `CONSOLIDATE_LOCI` | Funde hits do TBLASTN (clusterizados) com modelos do miniprot em loci candidatos | `loci_summary.tsv`, `candidate_proteins.fasta` |
| 5 | `HMMSEARCH_PEBP` | Confirma domínio PF01161 (PEBP) via HMMER | `pebp_confirmed.tsv`, `candidate_proteins_confirmed.fasta` |
| 6 | `DOWNLOAD_REF_PROTEOMES` | Baixa proteomas completos de Oryza/Zea/Arabidopsis | `<especie>.faa` + blastdb |
| 7 | `RECIPROCAL_BEST_HIT` | BLASTP recíproco — confirma ortologia real (não paráloga TFL1/MFT) | `rbh_summary.tsv` |
| 8 | `PHYLOGENY` | MAFFT + trimAl + IQ-TREE — candidatos no clado FT-like? | `urochloa_ft_phylo.treefile` |
| 9 | `REPORT` | Tabela final + heatmap + diagrama de loci | `final_candidates_table.tsv`, `figures/*.png` |

## Bloqueios de rede conhecidos

O servidor `200.235.143.10` só libera acesso de rede a `github.com` — todo
o resto (NCBI, EBI, etc.) é bloqueado ou não-garantido. Dois pontos do
pipeline dependem de download externo e têm fallback manual:

1. **`pfam_a_hmm`** (parâmetro em `nextflow.config`) aponta por padrão para
   `/home/eulalio/databases/pfam/Pfam-A.hmm`, já baixado pelo projeto
   irmão Kerson-paper no mesmo servidor — **reuso de dado de referência
   comum, não de código/execução** (compatível com a regra de isolamento
   entre projetos: repositórios e ambientes de execução nunca se misturam,
   mas bancos de dados de referência padrão podem ser compartilhados por
   path). Se esse arquivo não existir, o processo `HMMSEARCH_PEBP` tenta
   baixar o HMM de PF01161 via API do InterPro (best-effort, provavelmente
   bloqueado) e, se falhar, instrui a baixar manualmente no Windows e
   transferir via `scp`.
2. **`DOWNLOAD_REF_PROTEOMES`** tenta `datasets download genome accession ...`
   (NCBI Datasets CLI). Se falhar, baixe manualmente no Windows
   (NCBI Datasets → accession → Download → Protein FASTA), salve como
   `<nome_especie>.faa` e transfira via `scp` para
   `${params.local_ref_proteomes_dir}` (default
   `~/isa-lbmp/ref_proteomes_cache/`) antes de rodar — o processo usa esse
   cache local automaticamente se o arquivo já existir ali.

## Como estender a filogenia (Passo 10 do plano científico)

O domínio PEBP (PF01161) cobre FT, TFL1, MFT e CEN — genes com funções
opostas (florígeno vs. antiflorígeno). Para diferenciar corretamente FT-like
de TFL1-like/MFT-like na árvore, expanda o conjunto de referência
adicionando sequências dessas famílias das mesmas espécies e aponte
`--phylo_reference` para o novo FASTA, sem tocar em nenhum código:

```bash
nextflow run main.nf -profile mamba,local -params-file params.yaml \
    --phylo_reference reference_queries_expanded_tfl1_mft.fasta
```

## Isolamento do projeto

Este repositório é **totalmente independente** — `.git` próprio, sem
dependência de runtime de nenhum outro projeto do mesmo usuário. Convenções
(estrutura `modules/local/`, `conf/base.config`, uso de
`${System.getenv('HOME')}` em vez de `$HOME`) foram copiadas e adaptadas de
projetos irmãos (`MD-gromacs`, `soja-iac`), nunca importadas em runtime.
