# Identificação e validação in silico de ortólogos de FLOWERING LOCUS T (FT) em *Urochloa* spp.

> Documento vivo — atualizado a cada rodada de análise. Ver [Log de atualizações](#log-de-atualizações) no final.

## 1. Contexto e objetivo

Este documento registra a etapa de bioinformática (identificação e validação funcional *in silico* de ortólogos gênicos) do projeto "Modulação da transição vegetativo-reprodutiva em *Urochloa brizantha* por RNA de interferência exógeno: validação funcional de SBP5 e FT-like e impacto na produção de biomassa" (proposta completa em `projeto-isa.docx`, coordenação Drª Isabela Malaquias Dalto de Souza, LBMP/BIOAGRO-UFV).

O projeto maior visa silenciar transitoriamente, via aplicação exógena de dsRNA por pulverização foliar (não RNAi transgênico), os genes **SBP5** (família SPL, competência reprodutiva) e **FT-like** (florígeno, integrador final da via de indução floral) em *U. brizantha* cv. Marandu, testando a hipótese de que o silenciamento combinado prolonga a fase vegetativa e aumenta o acúmulo de biomassa.

**Este documento cobre exclusivamente o Objetivo Específico 1 da proposta**: *"Identificar regiões conservadas dos genes SBP5 e FT-like para o desenho de moléculas de dsRNA"* — restrito, por ora, ao gene **FT-like** (o gene SBP5 segue sem sequência de referência disponível; ver [§5.3](#53-lacunas-abertas)).

Como não existe assembly cromossômico próprio de *U. brizantha* (espécie agâmica/poliploide complexa), o genoma usado como alvo é o de **_Urochloa ruziziensis_** (assembly Embrapa_Uruz_1.0, 9 cromossomos), espécie diploide filogeneticamente próxima, adotada como *proxy*.

## 2. Metodologia

### 2.1 Pipeline de descoberta de loci (Nextflow DSL2, `main.nf`)

Pipeline reprodutível com 11 processos (repositório `github.com/eulaliobqi/Isa-LBMP`), executado em servidor Debian (32 cores) via `nextflow run main.nf -profile mamba,local -params-file params.yaml`, com ambientes conda/mamba criados automaticamente por processo:

| # | Processo | Função | Parâmetros-chave |
|---|---|---|---|
| 1 | `BUILD_BLASTDB` | Indexa o genoma alvo (`makeblastdb`) | — |
| 2 | `TBLASTN` | Triagem ampla proteína × genoma | e-value 1e-5 |
| 2b | `BLASTN_CDS` | Evidência cruzada CDS nucleotídeo × genoma | e-value 1e-3 |
| 3 | `MINIPROT_INDEX`/`ALIGN` | Estrutura gênica *spliced-aware* (introns) + tradução da proteína predita | — |
| 4 | `CONSOLIDATE_LOCI` | Funde hits das 3 fontes em loci candidatos | cluster a <20 kb |
| 5 | `HMMSEARCH_PEBP` | Confirma domínio PEBP (PF01161) via HMMER | e-value 1e-5 |
| 6 | `DOWNLOAD_REF_PROTEOMES` | Baixa proteomas completos de arroz/milho/Arabidopsis | assemblies RefSeq |
| 7 | `RECIPROCAL_BEST_HIT` | BLASTP recíproco — confirma ortologia real (não paralogia) | e-value 1e-10 |
| 8 | `PHYLOGENY` | MAFFT + trimAl + IQ-TREE | ver §2.3 |
| 9 | `REPORT` | Tabela final + figuras | — |

**Sequências de query FT-like (`reference_queries.fasta`)**: 4 proteínas de referência baixadas do NCBI — Hd3a (NP_001408118.1) e RFT1 (NP_001408117.1) de *Oryza sativa*, ZCN8 (NP_001106247.1) de *Zea mays*, FT (NP_001320342.1) de *Arabidopsis thaliana*.

### 2.2 Expansão da referência com homólogos reais de *Urochloa*

Para dar contexto filogenético mais próximo ao alvo (as 4 referências originais são de espécies distantes — arroz, milho, dicotiledônea), rodou-se BLASTP remoto (NCBI, banco `nr`) para cada uma das 4 proteínas de referência, restrito a `Urochloa[Organism]`, baixando as 5 sequências distintas mais próximas de cada uma (script `bin/blastp_top5_homologs.py`). Resultado: 20 hits, reduzidos a **15 sequências únicas** após remoção de duplicatas — todas de *U. humidicola* e *U. decumbens* (73–88% identidade); nenhuma de *U. ruziziensis*/*U. brizantha* especificamente, refletindo a escassez de proteínas anotadas dessas duas espécies em bancos públicos.

### 2.3 Filogenia com discriminação FT vs. TFL1/MFT

**Problema identificado**: o domínio PEBP (PF01161) é compartilhado por FT (florígeno), TFL1 (antiflorígeno) e MFT — genes de função oposta. Uma filogenia construída só com referências FT-like não consegue, por si, confirmar que um candidato é FT e não um parálogo antagonista.

**Solução**: busca por BLASTP (mesma técnica do §2.2) de sequências TFL1-like e MFT-like verificadas nas mesmas 3 espécies de referência (arroz, milho, Arabidopsis), usando como query TFL1 (NP_196004.1) e MFT (NP_173250.1) de Arabidopsis — ambos confirmados pelo campo `[Gene Name]` do NCBI, não por busca textual livre (que colidiu: o símbolo histórico `RCN1` de arroz foi reatribuído no NCBI atual a um gene não relacionado, um transportador ABC; script `bin/blastp_tfl1_mft_verify.py`).

Sequências TFL1-like e MFT-like identificadas (script `bin/blastp_tfl1_mft_verify.py`):

| Família | Arabidopsis | Arroz | Milho |
|---|---|---|---|
| TFL1-like | NP_196004.1 | XP_015624118.1 ("SELF-PRUNING"/CEN-like, 72,8% id.) | NP_001106241.1 (ZCN2, 74,0% id.) |
| MFT-like | NP_173250.1 | NP_001410922.1 (66,9% id.) | NP_001106249.1 (ZCN10, 63,0% id.) |

Conjunto de referência final para a filogenia (`reference_queries_expanded_tfl1_mft.fasta`, 25 sequências) = 4 FT-like originais + 15 homólogos de *Urochloa* (§2.2) + 3 TFL1-like + 3 MFT-like.

**Enraizamento**: MFT é basal ao *split* evolutivo FT/TFL1 (Wickland & Hanzawa, 2015). As 3 sequências MFT-like foram usadas como *outgroup* (`iqtree_outgroup` em `nextflow.config`/`params.yaml`, passado como `-o` ao IQ-TREE) — sem isso o Newick produzido é tecnicamente não-enraizado e sua topologia próxima à raiz não tem interpretação biológica confiável.

Alinhamento MAFFT (`--auto`) → aparo trimAl (`-automated1`) → IQ-TREE (`-m MFP`, `-bb 1000`, `-T AUTO`).

## 3. Resultados

### 3.1 Loci candidatos e confirmação por domínio PEBP

24 loci candidatos na triagem ampla (TBLASTN); **6 confirmam o domínio PEBP** via HMMER: `locus_001, 002, 005, 009, 012, 015`.

### 3.2 Reciprocal Best Hit (RBH) — confirmação de ortologia

Após correção de um bug de normalização de ID (`ref|...|` vs. accession puro — commit `9f96ac6`), rodagem limpa (sem `-resume`) confirmou:

| Locus | Evidências | Bitscore | Proteína | PEBP | RBH (3 espécies) |
|---|---|---|---|---|---|
| **locus_001** | tblastn+miniprot+blastn_cds | **155,0** | 179 aa | ✅ | ✅ **3/3** |
| locus_009 | tblastn+miniprot | 129,0 | 187 aa | ✅ | ✅ **3/3** |
| locus_015 | tblastn+miniprot | 129,0 | 169 aa | ✅ | ⚠️ 1/3 (só Arabidopsis) |
| locus_002 | tblastn | 123,0 | 76 aa (fragmento) | ✅ | ✅ (evidência fraca) |
| locus_005 | tblastn+miniprot | 122,0 | 175 aa | ✅ | ✅ |
| locus_012 | tblastn | 79,0 | 72 aa (fragmento) | ✅ | ✅ (2/3) |

`locus_015` só confirma reciprocidade com Arabidopsis — em arroz/milho, o melhor hit reverso cai em `locus_009` (parálogo quase idêntico, branch length ≈0 na árvore), sinal biológico de duplicação em tandem recente, não artefato de pipeline.

### 3.3 Filogenia (sem outgroup TFL1/MFT — versão anterior)

Análise da primeira árvore (`reference_queries_expanded_urochloa.fasta`, sem outgroup) revelou **3 clados FT-like distintos** em *Urochloa*:

1. **`locus_001`** + homólogos de Urochloa recuperados via BLASTP-Hd3a/RFT1 + o par Hd3a/RFT1 de arroz (bootstrap 69) — candidato mais próximo do florígeno "clássico".
2. **`locus_009`/`locus_015`** (bootstrap 90) + homólogos recuperados via BLASTP-Arabidopsis, próximo (bootstrap fraco 58) de `NP_001320342.1` e dos fragmentos curtos `locus_002`/`012`.
3. **`locus_005`** + `NP_001106247.1` (ZCN8) + homólogos recuperados via BLASTP-ZCN8 — bootstrap **100**, o clado mais bem suportado.

Sem outgroup, os nós profundos que conectam os 3 clados tinham suporte fraco (bootstrap 47–70), e a posição de `locus_001` "na base" da árvore era um artefato do Newick não-enraizado do IQ-TREE, sem interpretação biológica confiável.

### 3.4 Filogenia com outgroup MFT-like — CONCLUÍDA

Rodagem com `reference_queries_expanded_tfl1_mft.fasta` + `iqtree_outgroup` (as 3 MFT-like) concluída em 2026-07-10 no servidor (`-resume`, reaproveitando os processos anteriores ao `PHYLOGENY`). Árvore renderizada em `results/figures/phylogeny_tree_rooted_tfl1_mft.png`.

**A árvore agora é enraizada com significado biológico real** (raiz = `NP_173250.1`, MFT de Arabidopsis) e resolve de forma limpa os 3 grupos esperados:

1. **Clado MFT** (raiz): `NP_173250.1` (Arabidopsis) → (`NP_001410922.1` [arroz], bootstrap 62) → (**`locus_012`** + `NP_001106249.1` [ZCN10, milho], bootstrap **100**, branch length ≈0). **Achado novo**: `locus_012`, que tinha evidência fraca de RBH e era tratado como candidato FT-like fragmentado, na verdade **não é FT-like — é um ortólogo de MFT** (agrupa quase identicamente com ZCN10 de milho). Reclassificado.
2. **Clado TFL1-like** (bootstrap 98, sister ao clado FT): `NP_196004.1` (Arabidopsis TFL1) + (`XP_015624118.1` [arroz] + `NP_001106241.1` [milho, ZCN2], bootstrap 93). **Nenhum locus de Urochloa cai aqui.**
3. **Clado FT-like** (bootstrap 99, ingroup principal): contém `locus_001`, `locus_002`, `locus_005`, `locus_009`, `locus_015` — todos claramente do lado FT do split, separados do TFL1-like e do outgroup MFT. Dentro dele:
   - subclado ZCN8-type (bootstrap 100): `locus_005` + homólogos Urochloa + `NP_001106247.1` (ZCN8)
   - subclado Arabidopsis-FT/paralogos (bootstrap 56): `locus_002`+`NP_001320342.1` (84) e `locus_009`/`015`+homólogos Urochloa (100)
   - subclado Hd3a/RFT1-type (bootstrap 82): **`locus_001`** + `NP_001408118.1`(Hd3a)/`NP_001408117.1`(RFT1) + homólogos Urochloa densamente agrupados (92–95) — posição fina de `locus_001` entre os homólogos mais próximos tem suporte fraco (30–50, esperado dado a altíssima similaridade entre eles), mas sua alocação no clado FT-like como um todo é robusta.

**Critério de sucesso do §3.4 anterior confirmado**: `locus_001` cai dentro do clado FT, separado do TFL1-like e do outgroup MFT-like.

## 4. Discussão

- **`locus_001`** (CM027126.1:~30.720.251–30.721.504, fita +) é o candidato confirmado a ortólogo de FT em *U. ruziziensis*: única confirmação com as 3 fontes de evidência simultâneas (TBLASTN+BLASTN_CDS+miniprot), maior bitscore (155), proteína completa (179 aa), RBH limpo (3/3 espécies, sem ambiguidade de parálogo) **e agora também confirmação filogenética direta** — cai dentro do clado FT-like (bootstrap 99), claramente separado do clado TFL1-like (bootstrap 98) e do outgroup MFT-like, no subclado Hd3a/RFT1-type (bootstrap 82) junto com os homólogos mais próximos de *Urochloa*. **Está pronto para ser reportado como candidato publicável.**
- A existência de **3 clados FT-like distintos dentro do ingroup confirmado** (`locus_001`-type/Hd3a-RFT1, `locus_009/015`-type, `locus_005`≈ZCN8-type) é biologicamente esperada — gramíneas costumam ter pequenas famílias de genes FT-like (arroz: Hd3a+RFT1; milho: múltiplos ZCN) — e reforça que a triagem ampla capturou a diversidade real do locus, não ruído.
- **Reclassificação de `locus_012`**: a suspeita levantada pela coincidência do RBH (§3.2, `NP_001106249.1` como forward-hit) foi **confirmada pela filogenia enraizada** — `locus_012` cai dentro do clado MFT (bootstrap 100, quase idêntico a ZCN10 de milho), não do clado FT-like. Não deve mais ser tratado como candidato FT-like (nem como "fragmento fraco" da família FT) — é um ortólogo de MFT, com implicação biológica distinta (MFT regula germinação/dormência de semente e integra sinalização floral de forma indireta, não é o florígeno em si).
- **Fragmento curto restante**: `locus_002` (76 aa) segue com evidência mais fraca (só TBLASTN) dentro do clado FT-like — interpretar com cautela; pode ser um gene real truncado na montagem, pseudogene, ou artefato de triagem.

### 5.1 Próximo passo imediato

Com `locus_001` confirmado (evidência de sequência + RBH + filogenia enraizada), o próximo passo é **consolidar esse resultado como candidato final** para a etapa de desenho de dsRNA (extrair a região codificante completa, identificar segmentos conservados vs. específicos para o Objetivo Específico 1 da proposta).

### 5.2 Consideração metodológica pendente

Avaliar se `locus_009`/`locus_015` (segundo clado FT-like, bem suportado mas com parálogo quase idêntico) merece investigação adicional — pode ser uma segunda cópia funcional relevante para o desenho de dsRNA (maior cobertura) ou uma duplicação recente sem relevância funcional distinta.

### 5.3 Lacunas abertas

- **SBP5**: nenhuma sequência de referência identificada ainda — bloqueia completamente a etapa equivalente de descoberta de loci para esse gene (Objetivo Específico 1 da proposta cobre os dois genes).
- **Especificidade off-target**: a proposta científica (`projeto-isa.docx`) pede explicitamente análise de especificidade do dsRNA contra o transcriptoma/genoma de *U. brizantha* antes do desenho final dos primers — etapa ainda não implementada no pipeline.
- **Synteny/colinearidade**: não avaliada; seria evidência adicional relevante se o projeto avançar para publicação (ex.: GENESPACE/MCScanX).

## Log de atualizações

| Data | O que mudou |
|---|---|
| 2026-07-10 | Criação do documento. Consolida: descoberta de loci, RBH confirmado (`locus_001` 3/3), expansão com homólogos de *Urochloa*, identificação de TFL1-like/MFT-like reais e configuração do outgroup. Filogenia com outgroup rodando no servidor, resultado pendente. |
| 2026-07-10 | **Filogenia com outgroup MFT concluída.** `locus_001` confirmado dentro do clado FT-like (bootstrap 99), separado de TFL1-like (98) e MFT (raiz) — candidato pronto para reportar. Achado novo: `locus_012` reclassificado de "FT-like fraco" para **MFT-like** (bootstrap 100 com ZCN10 de milho). §3.4, §4 e §5.1 atualizados. |
