process BLASTN_CDS {
    tag "blastn_cds_vs_urochloa"
    label 'process_high'
    conda "bioconda::blast=2.16.0"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path cds_query, stageAs: 'cds_query.fasta'
    path db_files

    output:
    path 'blastn_cds_raw.tsv', emit: tsv

    script:
    // Complementar ao TBLASTN: BLASTN nucleotídeo x nucleotídeo usando o CDS
    // (sem íntrons) diretamente como query. Sensibilidade menor que TBLASTN entre
    // táxons distantes (arroz/milho/Arabidopsis vs. Urochloa), mas serve como
    // evidência cruzada independente — loci confirmados pelas 3 fontes
    // (TBLASTN + BLASTN_CDS + miniprot) são o sinal mais forte de candidato real.
    """
    blastn \\
        -query cds_query.fasta \\
        -db urochloa_db \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \\
        -evalue ${params.blastn_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs ${params.tblastn_max_targets} \\
        -out blastn_cds_raw.tsv

    if [ ! -s blastn_cds_raw.tsv ]; then
        echo "AVISO: blastn (CDS) não retornou nenhum hit (e-value ${params.blastn_evalue})." >&2
    fi
    """
}
