process REPORT {
    tag "final_report"
    label 'process_low'
    // matplotlib=3.9.0 não existe de verdade no conda-forge (mais um pin
    // "chutado" errado); pandas=2.2.3/biopython=1.84 já provados em outros
    // módulos (HMMSEARCH_PEBP, CONSOLIDATE_LOCI, RECIPROCAL_BEST_HIT), mantidos.
    conda "conda-forge::pandas=2.2.3 conda-forge::matplotlib conda-forge::biopython=1.84"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path loci_summary, stageAs: 'loci_summary.tsv'
    path pebp_confirmed, stageAs: 'pebp_confirmed.tsv'
    path rbh_summary, stageAs: 'rbh_summary.tsv'
    path treefile, stageAs: 'urochloa_ft_phylo.treefile'

    output:
    path 'final_candidates_table.tsv', emit: table
    path 'figures/bitscore_heatmap.png', emit: heatmap
    path 'figures/loci_diagram.png', emit: diagram

    script:
    """
    mkdir -p figures
    python3 ${projectDir}/bin/build_report.py \\
        --loci-summary loci_summary.tsv \\
        --pebp-confirmed pebp_confirmed.tsv \\
        --rbh-summary rbh_summary.tsv \\
        --treefile urochloa_ft_phylo.treefile \\
        --out-table final_candidates_table.tsv \\
        --out-heatmap figures/bitscore_heatmap.png \\
        --out-diagram figures/loci_diagram.png
    """
}
