process MINIPROT_ALIGN {
    tag "miniprot_align"
    label 'process_high'
    conda "bioconda::miniprot=0.13 conda-forge::biopython=1.84"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path index, stageAs: 'genome.mpi'
    path genome, stageAs: 'genome.fasta'
    path query, stageAs: 'query_proteins.fasta'

    output:
    path 'miniprot_hits.gff3', emit: gff3
    path 'miniprot_predicted_proteins.fasta', emit: proteins

    script:
    // Tradução roda no mesmo processo do miniprot para não restagear o genoma
    // (583 MB) em uma segunda task — decisão de design confirmada no plano.
    //
    // miniprot --gff embute linhas de resumo "##PAF" (extensão própria, fora
    // do padrão GFF3) intercaladas com as features reais -- filtradas antes
    // de gravar o GFF3 final.
    //
    // NÃO usa gffread: miniprot tolera pequenos frameshifts dentro de um
    // éxon individual (comprimento não múltiplo de 3), mas a soma de todos
    // os CDS de um mRNA é múltipla de 3. gffread traduz respeitando a coluna
    // "phase" exon a exon e produz proteína fora de frame (cheia de stop
    // códon/X no meio) com esse padrão -- bin/translate_miniprot_gff.py
    // concatena todos os CDS em ordem genômica e traduz de uma vez, que é o
    // que reproduz corretamente o alinhamento pretendido pelo miniprot
    // (validado localmente com caso sintético antes de subir).
    """
    miniprot -t${task.cpus} --gff genome.mpi query_proteins.fasta \\
        | grep -v '^##PAF' > miniprot_hits.gff3

    if [ ! -s miniprot_hits.gff3 ]; then
        echo "AVISO: miniprot não retornou nenhum modelo gênico." >&2
        touch miniprot_predicted_proteins.fasta
    else
        python3 ${projectDir}/bin/translate_miniprot_gff.py \\
            --gff3 miniprot_hits.gff3 \\
            --genome genome.fasta \\
            --out-fasta miniprot_predicted_proteins.fasta
    fi
    """
}
