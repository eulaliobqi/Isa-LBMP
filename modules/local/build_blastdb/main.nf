process BUILD_BLASTDB {
    tag "urochloa_db"
    label 'process_medium'
    conda "bioconda::blast=2.16.0"
    publishDir "${params.outdir}/blastdb", mode: 'copy', pattern: 'blastdb_info.txt'

    input:
    path genome, stageAs: 'genome.fasta'

    output:
    path 'urochloa_db.*', emit: db
    path 'blastdb_info.txt', emit: info

    script:
    """
    makeblastdb \\
        -in genome.fasta \\
        -dbtype nucl \\
        -parse_seqids \\
        -title "Urochloa ruziziensis Embrapa_Uruz_1.0 (proxy diploide p/ U. brizantha)" \\
        -out urochloa_db

    blastdbcmd -db urochloa_db -info > blastdb_info.txt
    """
}
