process CONSOLIDATE_LOCI {
    tag "consolidate_loci"
    label 'process_low'
    conda "conda-forge::pandas=2.2.3 conda-forge::biopython=1.84"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tblastn_tsv, stageAs: 'tblastn_raw.tsv'
    path blastn_cds_tsv, stageAs: 'blastn_cds_raw.tsv'
    path miniprot_gff3, stageAs: 'miniprot_hits.gff3'
    path miniprot_proteins, stageAs: 'miniprot_predicted_proteins.fasta'

    output:
    path 'loci_summary.tsv', emit: loci_summary
    path 'candidate_proteins.fasta', emit: candidate_proteins

    script:
    """
    python3 ${projectDir}/bin/consolidate_loci.py \\
        --tblastn tblastn_raw.tsv \\
        --blastn-cds blastn_cds_raw.tsv \\
        --miniprot-gff3 miniprot_hits.gff3 \\
        --miniprot-proteins miniprot_predicted_proteins.fasta \\
        --cluster-distance-bp ${params.cluster_distance_bp} \\
        --out-summary loci_summary.tsv \\
        --out-fasta candidate_proteins.fasta
    """
}
