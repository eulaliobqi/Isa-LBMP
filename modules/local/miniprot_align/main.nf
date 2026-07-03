process MINIPROT_ALIGN {
    tag "miniprot_align"
    label 'process_high'
    conda "bioconda::miniprot=0.13 bioconda::gffread=0.12.7"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path index, stageAs: 'genome.mpi'
    path genome, stageAs: 'genome.fasta'
    path query, stageAs: 'query_proteins.fasta'

    output:
    path 'miniprot_hits.gff3', emit: gff3
    path 'miniprot_predicted_proteins.fasta', emit: proteins

    script:
    // gffread roda no mesmo processo do miniprot para não restagear o genoma
    // (583 MB) em uma segunda task — decisão de design confirmada no plano.
    """
    miniprot -t${task.cpus} --gff genome.mpi query_proteins.fasta > miniprot_hits.gff3

    if [ ! -s miniprot_hits.gff3 ]; then
        echo "AVISO: miniprot não retornou nenhum modelo gênico." >&2
        touch miniprot_predicted_proteins.fasta
    else
        gffread -y miniprot_predicted_proteins.fasta -g genome.fasta miniprot_hits.gff3
    fi
    """
}
