process TBLASTN {
    tag "tblastn_vs_urochloa"
    label 'process_high'
    conda "bioconda::blast=2.16.0"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path query, stageAs: 'query_proteins.fasta'
    path db_files

    output:
    path 'tblastn_raw.tsv', emit: tsv

    script:
    """
    tblastn \\
        -query query_proteins.fasta \\
        -db urochloa_db \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sseq" \\
        -evalue ${params.tblastn_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs ${params.tblastn_max_targets} \\
        -max_hsps ${params.tblastn_max_hsps} \\
        -out tblastn_raw.tsv

    if [ ! -s tblastn_raw.tsv ]; then
        echo "AVISO: tblastn não retornou nenhum hit (e-value ${params.tblastn_evalue})." >&2
    fi
    """
}
