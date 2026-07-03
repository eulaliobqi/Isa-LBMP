process PHYLOGENY {
    tag "phylogeny_ft_pebp"
    label 'process_medium'
    conda "bioconda::mafft=7.525 bioconda::trimal=1.5.0 bioconda::iqtree=2.3.6"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path candidates, stageAs: 'candidate_proteins_confirmed.fasta'
    path phylo_ref, stageAs: 'phylo_reference.fasta'

    output:
    path 'urochloa_ft_phylo.treefile', emit: treefile
    path 'urochloa_ft_phylo.iqtree', emit: report
    path 'trimmed.fasta', emit: trimmed_aln

    script:
    """
    if [ ! -s candidate_proteins_confirmed.fasta ]; then
        echo "ERRO: nenhum candidato confirmado (domínio PEBP) disponível para filogenia." >&2
        exit 1
    fi

    cat candidate_proteins_confirmed.fasta phylo_reference.fasta > combined.fasta
    mafft --auto --thread ${task.cpus} combined.fasta > aligned.fasta
    trimal -in aligned.fasta -out trimmed.fasta -automated1

    iqtree2 \\
        -s trimmed.fasta \\
        -m ${params.iqtree_model} \\
        -bb ${params.iqtree_bootstrap} \\
        -T AUTO \\
        -pre urochloa_ft_phylo \\
        -redo
    """
}
