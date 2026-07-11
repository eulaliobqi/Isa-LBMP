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
    // Outgroup MFT-like (basal ao split FT/TFL1 — Karlgren et al. 2011; de Jesus
    // et al. 2022) — tentativa de enraizar a árvore com significado biológico.
    // NOTA: com as 3 sequências MFT-like atuais, o IQ-TREE ainda reporta a árvore
    // como "UNROOTED" (as 3 não formam clado monofilético nesta topologia — ver
    // artigo.md §3.4). A separação FT-like vs. TFL1-like continua robusta mesmo
    // assim; o enraizamento pretendido, especificamente, não foi resolvido.
    // Lista de IDs (primeiro token do header FASTA, separados por vírgula) vem
    // de params.iqtree_outgroup.
    def outgroupFlag = params.iqtree_outgroup ? "-o ${params.iqtree_outgroup}" : ''
    """
    if [ ! -s candidate_proteins_confirmed.fasta ]; then
        echo "ERRO: nenhum candidato confirmado (domínio PEBP) disponível para filogenia." >&2
        exit 1
    fi

    cat candidate_proteins_confirmed.fasta phylo_reference.fasta > combined.fasta
    mafft --auto --thread ${task.cpus} combined.fasta > aligned.fasta
    trimal -in aligned.fasta -out trimmed.fasta -automated1

    # O binário se chama "iqtree2" a partir da v2.0.4, mas builds/versões
    # mais antigas do pacote bioconda "iqtree" só têm o binário "iqtree".
    # Sem pin de versão (necessário pra resolver o conflito de solver),
    # não dá pra saber de antemão qual vai ser instalado.
    IQTREE_BIN=\$(command -v iqtree2 || command -v iqtree)
    if [ -z "\$IQTREE_BIN" ]; then
        echo "ERRO: nem iqtree2 nem iqtree encontrados no PATH." >&2
        exit 1
    fi

    "\$IQTREE_BIN" \\
        -s trimmed.fasta \\
        -m ${params.iqtree_model} \\
        -bb ${params.iqtree_bootstrap} \\
        -T AUTO \\
        -pre urochloa_ft_phylo \\
        ${outgroupFlag} \\
        -redo
    """
}
