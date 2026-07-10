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

**Enraizamento**: MFT é basal ao *split* evolutivo FT/TFL1 — genes MFT-like estão presentes em angiospermas, gimnospermas, licófitas e briófitas, enquanto os clados FT-like e TFL1-like não incluem nenhum gene PEBP de musgo ou licófita (Karlgren et al., 2011), e múltiplos estudos convergem em tratar MFT como ancestral a FT e TFL1 (de Jesus et al., 2022; ver [§6.9](#6-justificativa-metodológica) para a citação completa e uma ressalva sobre o que cada uma dessas referências sustenta exatamente). As 3 sequências MFT-like foram usadas como *outgroup* (`iqtree_outgroup` em `nextflow.config`/`params.yaml`, passado como `-o` ao IQ-TREE) — sem isso o Newick produzido é tecnicamente não-enraizado e sua topologia próxima à raiz não tem interpretação biológica confiável.

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

## 5. Próximos passos e lacunas

### 5.1 Próximo passo imediato

Com `locus_001` confirmado (evidência de sequência + RBH + filogenia enraizada), o próximo passo é **consolidar esse resultado como candidato final** para a etapa de desenho de dsRNA (extrair a região codificante completa, identificar segmentos conservados vs. específicos para o Objetivo Específico 1 da proposta).

### 5.2 Consideração metodológica pendente

Avaliar se `locus_009`/`locus_015` (segundo clado FT-like, bem suportado mas com parálogo quase idêntico) merece investigação adicional — pode ser uma segunda cópia funcional relevante para o desenho de dsRNA (maior cobertura) ou uma duplicação recente sem relevância funcional distinta.

### 5.3 Lacunas abertas

- **SBP5**: nenhuma sequência de referência identificada ainda — bloqueia completamente a etapa equivalente de descoberta de loci para esse gene (Objetivo Específico 1 da proposta cobre os dois genes).
- **Especificidade off-target**: a proposta científica (`projeto-isa.docx`) pede explicitamente análise de especificidade do dsRNA contra o transcriptoma/genoma de *U. brizantha* antes do desenho final dos primers — etapa ainda não implementada no pipeline.
- **Synteny/colinearidade**: não avaliada; seria evidência adicional relevante se o projeto avançar para publicação (ex.: GENESPACE/MCScanX).

## 6. Justificativa metodológica

Cada escolha metodológica relevante do pipeline foi checada contra a literatura (busca + verificação da fonte real antes de citar — nenhuma citação abaixo veio de memória sem confirmação). Itens de engenharia pura (estratégia de retry, labels de recurso computacional, `--channel-priority flexible` do conda, ausência de pin de versão nos envs, distância de 20 kb para fusão de loci, score composto `n_evidencias` como soma de 3 flags) não têm — nem precisam de — citação: são decisões pragmáticas de implementação, documentadas como tal por transparência, não por falta de rigor.

### 6.1 Genoma de espécie-proxy (*U. ruziziensis* para *U. brizantha*)

*U. brizantha* é poliploide agâmica sem assembly cromossômico próprio; usar o genoma de *U. ruziziensis* (diploide, parente próximo com assembly disponível) como alvo é consistente com o princípio geral, bem estabelecido em genômica de poliploides, de que genomas de parentes diploides próximos servem de referência para análise gênica/genômica em táxons poliploides sem assembly completo — ex. o uso de genomas de progenitores diploides de trigo hexaploide para reconstruir composição e história do genoma poliploide (Marcussen et al., 2014). A proximidade genômica específica entre *U. ruziziensis* e *U. brizantha*/demais espécies do gênero (incluindo evidência de que cromossomos de *U. decumbens* derivam parcialmente de *U. ruziziensis* ou de *U. brizantha* diploide) está documentada em Ferreira et al. (2021), já citado na proposta do projeto.

**Ressalva honesta**: não encontrei um estudo que faça exatamente o mesmo desenho (genoma de espécie congênere inteiro como substituto de alvo de descoberta gênica em poliploide agâmico sem assembly) — o suporte acima é por analogia a um princípio geral bem documentado (genoma de parente diploide próximo é referência válida na ausência de assembly do próprio poliplóide), não uma réplica direta de metodologia.

### 6.2 Genes de referência (Hd3a, RFT1, ZCN8, FT)

- Kojima, S., Takahashi, Y., Kobayashi, Y., et al. (2002). Hd3a, a Rice Ortholog of the Arabidopsis FT Gene, Promotes Transition to Flowering Downstream of Hd1 under Short-Day Conditions. *Plant and Cell Physiology*, 43, 1096–1105. https://doi.org/10.1093/pcp/pcf156 — caracteriza Hd3a como ortólogo funcional de FT em arroz, base para seu uso como gene-semente.
- Danilevskaya, O. N., Meng, X., Hou, Z., Ananiev, E. V., & Simmons, C. R. (2008). A Genomic and Expression Compendium of the Expanded PEBP Gene Family from Maize. *Plant Physiology*, 146, 250–264. https://doi.org/10.1104/pp.107.109538 — identifica e nomeia a família ZCN em milho (25 genes), incluindo ZCN8.
- Karlgren, A., et al. (2011). Evolution of the PEBP Gene Family in Plants: Functional Diversification in Seed Plant Evolution. *Plant Physiology*, 156(4), 1967–1977. https://doi.org/10.1104/pp.111.176206 — revisão abrangente da família PEBP em plantas (19 cópias em arroz, 24 em milho vs. 6 em Arabidopsis), reforça que Hd3a/RFT1/ZCN8/FT são membros bem caracterizados do clado FT-like nas 3 espécies-modelo usadas.

### 6.3 E-value do TBLASTN (1e-5) e do BLASTN de CDS (1e-3)

Não encontrei um único paper dedicado exclusivamente a validar 1e-5 como corte ótimo para TBLASTN entre táxons divergentes — é prática amplamente adotada (frequentemente citada na faixa 1e-4 a 1e-6 em fluxos de detecção de ortólogos por BLAST recíproco) e consistente com a discussão de escolha de parâmetros de BLAST para detecção de ortólogos em Moreno-Hagelsieb & Latimer (2008, ver §6.5). O corte mais permissivo do BLASTN de CDS (1e-3) reflete a menor sensibilidade conhecida de buscas nucleotídeo×nucleotídeo frente a buscas baseadas em proteína para detectar homologia remota entre táxons divergentes — princípio geral de comparação de sequências (a maior degenerescência do código genético preserva mais sinal em nível de proteína que de nucleotídeo em posições sinônimas), não uma citação específica de threshold. **Marcado honestamente como suporte geral/de boas práticas, não validação direta do valor exato**.

### 6.4 miniprot (alinhamento spliced-aware proteína→genoma)

Li, H. (2023). Protein-to-genome alignment with miniprot. *Bioinformatics*, 39(1), btad014. https://doi.org/10.1093/bioinformatics/btad014 — ferramenta usada exatamente para o propósito do pipeline: mapear proteínas de referência contra um genoma novo respeitando estrutura de éxons/íntrons (splicing) e frameshifts, para anotação de genes codificantes usando genes conhecidos de outras espécies.

### 6.5 Convergência de múltiplas fontes de evidência independentes

Allen, J. E., Pertea, M., & Salzberg, S. L. (2004). Computational Gene Prediction Using Multiple Sources of Evidence. *Genome Research*, 14(1), 142–148. https://doi.org/10.1101/gr.1562804 — demonstra que combinar evidências independentes (preditores *ab initio*, alinhamentos de proteína, EST/cDNA) supera consistentemente qualquer fonte isolada em sensibilidade e especificidade — sustenta diretamente a heurística do pipeline de tratar a convergência TBLASTN+BLASTN_CDS+miniprot como sinal mais forte que qualquer fonte isolada.

### 6.6 HMMER / Pfam (domínio PEBP, PF01161)

- Eddy, S. R. (2011). Accelerated Profile HMM Searches. *PLOS Computational Biology*, 7(10), e1002195. https://doi.org/10.1371/journal.pcbi.1002195 — descreve o HMMER3, ferramenta usada no `HMMSEARCH_PEBP`.
- Mistry, J., Chuguransky, S., Williams, L., Qureshi, M., Salazar, G. A., Sonnhammer, E. L. L., et al. (2021). Pfam: The Protein Families Database in 2021. *Nucleic Acids Research*, 49(D1), D412–D419. https://doi.org/10.1093/nar/gkaa913 — banco de domínios usado como fonte do perfil HMM PF01161.

### 6.7 Reciprocal Best Hit (RBH) como método de ortologia — e suas limitações com parálogos

- Tatusov, R. L., Koonin, E. V., & Lipman, D. J. (1997). A Genomic Perspective on Protein Families. *Science*, 278(5338), 631–637. https://doi.org/10.1126/science.278.5338.631 — trabalho fundacional que popularizou a lógica de melhor-hit-recíproco para inferência de ortologia (base do banco COG).
- Moreno-Hagelsieb, G., & Latimer, K. (2008). Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics*, 24(3), 319–324. https://doi.org/10.1093/bioinformatics/btm585 — discute diretamente a escolha de parâmetros de BLAST (filtragem, algoritmo de alinhamento) que mais afetam a detecção de ortólogos via RBH.
- **Dalquen, D. A., & Dessimoz, C. (2013). Bidirectional Best Hits Miss Many Orthologs in Duplication-Rich Clades such as Plants and Animals. *Genome Biology and Evolution*, 5(10), 1800–1806. https://doi.org/10.1093/gbe/evt132** — mostra que RBH/BBH pode perder até ~60% das relações de ortologia verdadeiras em clados ricos em duplicação, e que quando há duplicação independente (m cópias numa espécie, n em outra), RBH só consegue identificar no máximo min(m,n) pares, deixando os demais parálogos sem confirmação recíproca. **Esta é a citação mais valiosa da varredura**: explica exatamente, com literatura real, por que `locus_015` (parálogo quase idêntico a `locus_009`) só confirma RBH em 1 de 3 espécies (§3.2) — não é falha do pipeline, é uma limitação conhecida e documentada do método RBH clássico diante de duplicação em tandem recente.

### 6.8 Rigor do E-value do RBH (1e-10)

Não encontrei uma citação que valide especificamente 1e-10 como corte ótimo para a etapa de confirmação de ortologia (em oposição aos 1e-5/1e-3 da triagem ampla). A prática de aumentar o rigor estatístico quando o objetivo migra de "triagem ampla" para "confirmação 1:1" é consistente com a discussão de Moreno-Hagelsieb & Latimer (2008, §6.7), mas **marcado honestamente como suporte indireto/parcial**, não validação direta do valor exato 1e-10.

### 6.9 Filogenia: MAFFT, trimAl, ModelFinder, UFBoot2, outgroup MFT

- Katoh, K., & Standley, D. M. (2013). MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. *Molecular Biology and Evolution*, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010
- Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*, 25(15), 1972–1973. https://doi.org/10.1093/bioinformatics/btp348
- Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A., & Jermiin, L. S. (2017). ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods*, 14(6), 587–589. https://doi.org/10.1038/nmeth.4285
- Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q., & Vinh, L. S. (2018). UFBoot2: Improving the Ultrafast Bootstrap Approximation. *Molecular Biology and Evolution*, 35(2), 518–522. https://doi.org/10.1093/molbev/msx281 — importante notar que o suporte UFBoot2 (bootstrap ultrarrápido) se interpreta de forma diferente do bootstrap não-paramétrico clássico (valores UFBoot2 tendem a já ser menos inflacionados que o bootstrap padrão, mas continuam sendo uma aproximação, não o bootstrap exato).
- Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular Biology and Evolution*, 37(5), 1530–1534. https://doi.org/10.1093/molbev/msaa015

**Outgroup MFT** — duas fontes reais e diretamente checadas dão suporte, com uma ressalva importante:
- Karlgren, A., et al. (2011), já citado em §6.2 — mostra que genes MFT-like ocorrem em angiospermas, gimnospermas, licófitas *e* briófitas, enquanto os clados FT-like e TFL1-like estão ausentes em musgos e licófitas — evidência direta de que MFT é a linhagem mais basal da família PEBP.
- de Jesus, D. A., Batista, D. M., Monteiro, E. F., Salzman, S., Carvalho, L. M., Santana, K., & André, T. (2022). Structural changes and adaptive evolutionary constraints in FLOWERING LOCUS T and TERMINAL FLOWER1-like genes of flowering plants. *Frontiers in Genetics*, 13, 954015. https://doi.org/10.3389/fgene.2022.954015 — afirma explicitamente que "genes MFT são ancestrais a FT e TFL1" (citando estudos prévios), sustentando o uso de MFT como outgroup. **Ressalva**: esse mesmo estudo, na prática, usou FT/TFL1 de gimnospermas (não MFT) como outgroup, por falta de sequências MFT disponíveis para as espécies que analisaram — ou seja, apoia o princípio biológico (MFT é basal), mas não é um precedente de bitola idêntica de implementação.
- **Correção**: a versão anterior deste documento (§2.3, primeira edição) atribuía essa lógica de enraizamento a Wickland & Hanzawa (2015). Essa citação é real (Molecular Plant, 8(7), 983–997, DOI 10.1016/j.molp.2015.01.007) e é uma revisão de qualidade sobre a família FT/TFL1, mas seu resumo não faz a afirmação específica sobre MFT ser basal/apropriado como outgroup — mantida como referência geral da biologia de FT/TFL1, mas a citação específica do argumento de enraizamento foi corrigida para Karlgren (2011) e de Jesus et al. (2022) acima.

## Referências

Allen, J. E., Pertea, M., & Salzberg, S. L. (2004). Computational Gene Prediction Using Multiple Sources of Evidence. *Genome Research*, 14(1), 142–148. https://doi.org/10.1101/gr.1562804

Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*, 25(15), 1972–1973. https://doi.org/10.1093/bioinformatics/btp348

Dalquen, D. A., & Dessimoz, C. (2013). Bidirectional Best Hits Miss Many Orthologs in Duplication-Rich Clades such as Plants and Animals. *Genome Biology and Evolution*, 5(10), 1800–1806. https://doi.org/10.1093/gbe/evt132

Danilevskaya, O. N., Meng, X., Hou, Z., Ananiev, E. V., & Simmons, C. R. (2008). A Genomic and Expression Compendium of the Expanded PEBP Gene Family from Maize. *Plant Physiology*, 146, 250–264. https://doi.org/10.1104/pp.107.109538

de Jesus, D. A., Batista, D. M., Monteiro, E. F., Salzman, S., Carvalho, L. M., Santana, K., & André, T. (2022). Structural changes and adaptive evolutionary constraints in FLOWERING LOCUS T and TERMINAL FLOWER1-like genes of flowering plants. *Frontiers in Genetics*, 13, 954015. https://doi.org/10.3389/fgene.2022.954015

Eddy, S. R. (2011). Accelerated Profile HMM Searches. *PLOS Computational Biology*, 7(10), e1002195. https://doi.org/10.1371/journal.pcbi.1002195

Ferreira, R. C. U., et al. (2021). An Overview of the Genetics and Genomics of the Urochloa Species Most Commonly Used in Pastures. *Frontiers in Plant Science*, 12, 770461. https://doi.org/10.3389/fpls.2021.770461

Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q., & Vinh, L. S. (2018). UFBoot2: Improving the Ultrafast Bootstrap Approximation. *Molecular Biology and Evolution*, 35(2), 518–522. https://doi.org/10.1093/molbev/msx281

Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A., & Jermiin, L. S. (2017). ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods*, 14(6), 587–589. https://doi.org/10.1038/nmeth.4285

Karlgren, A., Gyllenstrand, N., Källman, T., Sundström, J. F., Moore, D., Lascoux, M., & Lagercrantz, U. (2011). Evolution of the PEBP Gene Family in Plants: Functional Diversification in Seed Plant Evolution. *Plant Physiology*, 156(4), 1967–1977. https://doi.org/10.1104/pp.111.176206

Katoh, K., & Standley, D. M. (2013). MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. *Molecular Biology and Evolution*, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010

Kojima, S., Takahashi, Y., Kobayashi, Y., et al. (2002). Hd3a, a Rice Ortholog of the Arabidopsis FT Gene, Promotes Transition to Flowering Downstream of Hd1 under Short-Day Conditions. *Plant and Cell Physiology*, 43, 1096–1105. https://doi.org/10.1093/pcp/pcf156

Li, H. (2023). Protein-to-genome alignment with miniprot. *Bioinformatics*, 39(1), btad014. https://doi.org/10.1093/bioinformatics/btad014

Marcussen, T., Sandve, S. R., Heier, L., Spannagl, M., Pfeifer, M., The International Wheat Genome Sequencing Consortium, et al. (2014). Ancient hybridizations among the ancestral genomes of bread wheat. *Science*, 345(6194), 1250092. https://doi.org/10.1126/science.1250092

Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular Biology and Evolution*, 37(5), 1530–1534. https://doi.org/10.1093/molbev/msaa015

Mistry, J., Chuguransky, S., Williams, L., Qureshi, M., Salazar, G. A., Sonnhammer, E. L. L., et al. (2021). Pfam: The Protein Families Database in 2021. *Nucleic Acids Research*, 49(D1), D412–D419. https://doi.org/10.1093/nar/gkaa913

Moreno-Hagelsieb, G., & Latimer, K. (2008). Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics*, 24(3), 319–324. https://doi.org/10.1093/bioinformatics/btm585

Tatusov, R. L., Koonin, E. V., & Lipman, D. J. (1997). A Genomic Perspective on Protein Families. *Science*, 278(5338), 631–637. https://doi.org/10.1126/science.278.5338.631

Wickland, D. P., & Hanzawa, Y. (2015). The FLOWERING LOCUS T/TERMINAL FLOWER 1 Gene Family: Functional Evolution and Molecular Mechanisms. *Molecular Plant*, 8(7), 983–997. https://doi.org/10.1016/j.molp.2015.01.007

## Log de atualizações

| Data | O que mudou |
|---|---|
| 2026-07-10 | Criação do documento. Consolida: descoberta de loci, RBH confirmado (`locus_001` 3/3), expansão com homólogos de *Urochloa*, identificação de TFL1-like/MFT-like reais e configuração do outgroup. Filogenia com outgroup rodando no servidor, resultado pendente. |
| 2026-07-10 | **Filogenia com outgroup MFT concluída.** `locus_001` confirmado dentro do clado FT-like (bootstrap 99), separado de TFL1-like (98) e MFT (raiz) — candidato pronto para reportar. Achado novo: `locus_012` reclassificado de "FT-like fraco" para **MFT-like** (bootstrap 100 com ZCN10 de milho). §3.4, §4 e §5.1 atualizados. |
| 2026-07-10 | **Varredura de literatura adicionada (§6 + Referências).** 20 citações reais verificadas (WebSearch+WebFetch, nenhuma de memória) dando suporte a cada escolha metodológica do pipeline. Destaque: Dalquen & Dessimoz (2013) explica com literatura real por que `locus_015` só confirma RBH em 1/3 espécies (limitação conhecida do RBH clássico com duplicação em tandem, não bug). Corrigida a atribuição da lógica de enraizamento por MFT: de Wickland & Hanzawa (2015, mantida como referência geral) para Karlgren et al. (2011) + de Jesus et al. (2022), que sustentam a afirmação específica. Itens sem citação direta (thresholds exatos de E-value) marcados honestamente como suporte parcial/de boas práticas. |
