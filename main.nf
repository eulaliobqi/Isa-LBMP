#!/usr/bin/env nextflow
// ============================================================================
// Isa-LBMP — Busca de homólogos FT-like (PEBP, PF01161) em Urochloa ruziziensis
// Autor: Eulalio Santos | UFV | Colaboradora: Isa (LBMP)
//
// Plano científico completo:
//   C:\Users\eulal\.claude\plans\objetivo-principal-usar-o-synthetic-peach.md
//
// Passos 1-11 do plano → processos deste workflow:
//   1. BUILD_BLASTDB          — makeblastdb no genoma Urochloa
//   2. TBLASTN                — triagem ampla (proteína×nucleotídeo)
//   2b. BLASTN_CDS            — evidência cruzada (CDS nucleotídeo×nucleotídeo)
//   3. MINIPROT_INDEX/ALIGN   — estrutura gênica spliced-aware + gffread
//   4. CONSOLIDATE_LOCI       — funde TBLASTN + BLASTN_CDS + miniprot em loci candidatos
//   5. HMMSEARCH_PEBP         — confirma domínio PEBP (PF01161)
//   6. DOWNLOAD_REF_PROTEOMES + RECIPROCAL_BEST_HIT — confirma ortologia (RBH)
//   7. PHYLOGENY              — mafft + trimal + iqtree2 (clado FT-like)
//   8. REPORT                 — tabela final + heatmap + diagrama de loci
//
// IMPORTANTE (regras não-negociáveis do setup Eulalio/UFV):
//   - screen -S isa-lbmp ANTES de rodar (SIGTTOU mata processo em background)
//   - mamba, nunca conda
//   - git pull antes de rodar (sincronizar com último fix commitado)
//   - ${System.getenv('HOME')} nos configs, nunca $HOME literal
// ============================================================================

nextflow.enable.dsl = 2

include { BUILD_BLASTDB          } from './modules/local/build_blastdb/main.nf'
include { TBLASTN                } from './modules/local/tblastn/main.nf'
include { BLASTN_CDS             } from './modules/local/blastn_cds/main.nf'
include { MINIPROT_INDEX         } from './modules/local/miniprot_index/main.nf'
include { MINIPROT_ALIGN         } from './modules/local/miniprot_align/main.nf'
include { CONSOLIDATE_LOCI       } from './modules/local/consolidate_loci/main.nf'
include { HMMSEARCH_PEBP         } from './modules/local/hmmsearch_pebp/main.nf'
include { DOWNLOAD_REF_PROTEOMES } from './modules/local/download_ref_proteomes/main.nf'
include { RECIPROCAL_BEST_HIT    } from './modules/local/reciprocal_best_hit/main.nf'
include { PHYLOGENY              } from './modules/local/phylogeny/main.nf'
include { REPORT                 } from './modules/local/report/main.nf'

// ── Validação de parâmetros ──────────────────────────────────────────────────
def validate_params() {
    if (!params.genome) {
        error "ERRO: --genome é obrigatório (FASTA multi-cromossomo do genoma Urochloa)."
    }
    if (!params.query_proteins) {
        error "ERRO: --query_proteins é obrigatório (proteínas FT-like de referência)."
    }
    if (!params.cds_query) {
        error "ERRO: --cds_query é obrigatório (CDS nucleotídeo dos mesmos genes de referência, para BLASTN complementar)."
    }
    if (!params.ref_proteomes || params.ref_proteomes.size() == 0) {
        error "ERRO: params.ref_proteomes vazio — necessário para RECIPROCAL_BEST_HIT."
    }
}

workflow {
    validate_params()

    log.info """
    ═══════════════════════════════════════════════════════════════
     Isa-LBMP — Homólogos FT-like (PEBP) em Urochloa ruziziensis
    ═══════════════════════════════════════════════════════════════
     Genoma alvo        : ${params.genome}
     Proteínas de ref.   : ${params.query_proteins}
     CDS de ref. (BLASTN): ${params.cds_query}
     Referência filogenia: ${params.phylo_reference ?: '(usa query_proteins)'}
     Espécies p/ RBH     : ${params.ref_proteomes.collect { it.name }.join(', ')}
     Saída               : ${params.outdir}
    ═══════════════════════════════════════════════════════════════
    """.stripIndent()

    // Arquivos de entrada reutilizados por múltiplos processos → value channels
    ch_genome = Channel.value(file(params.genome, checkIfExists: true))
    ch_query  = Channel.value(file(params.query_proteins, checkIfExists: true))
    ch_cds    = Channel.value(file(params.cds_query, checkIfExists: true))

    ch_phylo_ref = params.phylo_reference
        ? Channel.value(file(params.phylo_reference, checkIfExists: true))
        : ch_query

    // ── Passos 1-2: triagem BLAST (proteína + CDS complementar) ─────────────
    BUILD_BLASTDB(ch_genome)
    TBLASTN(ch_query, BUILD_BLASTDB.out.db)
    BLASTN_CDS(ch_cds, BUILD_BLASTDB.out.db)

    // ── Passo 3: estrutura gênica spliced-aware (miniprot + gffread) ───────
    MINIPROT_INDEX(ch_genome)
    MINIPROT_ALIGN(MINIPROT_INDEX.out.index, ch_genome, ch_query)

    // ── Passo 4: consolidação de loci candidatos (3 fontes de evidência) ────
    CONSOLIDATE_LOCI(
        TBLASTN.out.tsv,
        BLASTN_CDS.out.tsv,
        MINIPROT_ALIGN.out.gff3,
        MINIPROT_ALIGN.out.proteins
    )

    // ── Passo 5: confirmação do domínio PEBP (PF01161) ──────────────────────
    HMMSEARCH_PEBP(CONSOLIDATE_LOCI.out.candidate_proteins)

    // ── Passo 6: ortologia via Reciprocal Best Hit ───────────────────────────
    ch_ref_species  = Channel.fromList(params.ref_proteomes)
    DOWNLOAD_REF_PROTEOMES(ch_ref_species)

    ch_proteome_faa = DOWNLOAD_REF_PROTEOMES.out.faa.collect()
    ch_proteome_db  = DOWNLOAD_REF_PROTEOMES.out.db.collect()

    RECIPROCAL_BEST_HIT(
        HMMSEARCH_PEBP.out.confirmed_proteins,
        ch_query,
        ch_proteome_faa,
        ch_proteome_db
    )

    // ── Passo 7: filogenia (candidatos confirmados + referência expansível) ─
    PHYLOGENY(
        HMMSEARCH_PEBP.out.confirmed_proteins,
        ch_phylo_ref
    )

    // ── Passo 8: relatório final consolidado ────────────────────────────────
    REPORT(
        CONSOLIDATE_LOCI.out.loci_summary,
        HMMSEARCH_PEBP.out.confirmed_tsv,
        RECIPROCAL_BEST_HIT.out.rbh_summary,
        PHYLOGENY.out.treefile
    )
}

// Forma de atribuição (workflow.onComplete = { ... }), não a forma de statement
// solto (workflow.onComplete { ... }) — a segunda quebra em Nextflow >=24 (strict
// syntax rejeita statement fora de bloco) e, quando movida para dentro do
// workflow{} de entrada, falha em tempo de execução (escopo do onComplete não
// deve ficar aninhado dentro do próprio workflow que ele observa).
workflow.onComplete = {
    log.info """
    ═══════════════════════════════════════════════════════════════
     Isa-LBMP finalizado — status: ${workflow.success ? 'OK' : 'FALHOU'}
     Duração: ${workflow.duration}
     Resultados: ${params.outdir}
    ═══════════════════════════════════════════════════════════════
    """.stripIndent()
}
