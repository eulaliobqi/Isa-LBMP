process MINIPROT_INDEX {
    tag "miniprot_index"
    label 'process_high'
    conda "bioconda::miniprot=0.13"

    input:
    path genome, stageAs: 'genome.fasta'

    output:
    path 'genome.mpi', emit: index

    script:
    """
    miniprot -t${task.cpus} -d genome.mpi genome.fasta
    """
}
