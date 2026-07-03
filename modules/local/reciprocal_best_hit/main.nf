process RECIPROCAL_BEST_HIT {
    tag "reciprocal_best_hit"
    label 'process_medium'
    conda "bioconda::blast=2.16.0 conda-forge::pandas=2.2.3 conda-forge::biopython=1.84"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path candidates, stageAs: 'candidate_proteins_confirmed.fasta'
    path query, stageAs: 'query_proteins.fasta'
    path proteome_faas
    path proteome_dbs

    output:
    path 'rbh_summary.tsv', emit: rbh_summary

    script:
    def species = params.ref_proteomes.collect { it.name }.join(' ')
    """
    cat candidate_proteins_confirmed.fasta query_proteins.fasta > combined_reference.fasta
    # Sem -parse_seqids: os headers de candidate_proteins_confirmed.fasta (gerados
    # por bin/consolidate_loci.py, formato "locus_NNN|proteina|source=...") passam
    # de 50 caracteres, limite do BLAST para IDs quando -parse_seqids é usado.
    # Não precisa aqui -- combined_reference_db só é alvo de blastp (nunca de
    # blastdbcmd -entry/-entry_batch, que é a única razão de usar -parse_seqids).
    makeblastdb -in combined_reference.fasta -dbtype prot -out combined_reference_db

    for sp in ${species}; do
        # 1) candidatos -> proteoma da espécie (melhor hit direto)
        blastp -query candidate_proteins_confirmed.fasta -db \${sp}_db \\
            -outfmt "6 qseqid sseqid pident length evalue bitscore" \\
            -max_target_seqs 1 -evalue ${params.rbh_evalue} -num_threads ${task.cpus} \\
            -out fwd_\${sp}.tsv

        if [ -s fwd_\${sp}.tsv ]; then
            cut -f2 fwd_\${sp}.tsv | sort -u > fwd_\${sp}_hits.txt
            blastdbcmd -db \${sp}_db -entry_batch fwd_\${sp}_hits.txt -out fwd_\${sp}_hits.fasta

            # 2) melhores hits da espécie -> de volta contra candidatos + queries originais
            blastp -query fwd_\${sp}_hits.fasta -db combined_reference_db \\
                -outfmt "6 qseqid sseqid pident length evalue bitscore" \\
                -max_target_seqs 1 -evalue ${params.rbh_evalue} -num_threads ${task.cpus} \\
                -out rev_\${sp}.tsv
        else
            echo "AVISO: nenhum hit de ${species} encontrado para \${sp} (e-value ${params.rbh_evalue})." >&2
            touch rev_\${sp}.tsv
        fi
    done

    python3 ${projectDir}/bin/reciprocal_best_hit.py \\
        --species ${species} \\
        --out-tsv rbh_summary.tsv
    """
}
