process PHYLOGENY {
    tag "phylogeny_ft_pebp"
    label 'process_medium'
    // Sem pin de versão: os 3 pacotes juntos com versões exatas fixas geram
    // conflito de dependências transitivas no solver (mamba "is not
    // installable because it conflicts with any installable versions
    // previously reported"); deixar o mamba escolher versões mutuamente
    // compatíveis é mais robusto aqui do que reproducibilidade exata.
    conda "bioconda::mafft bioconda::trimal bioconda::iqtree"
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
