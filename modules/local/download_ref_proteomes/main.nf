process DOWNLOAD_REF_PROTEOMES {
    tag "${meta.name}"
    label 'process_low'
    conda "conda-forge::ncbi-datasets-cli=18.31.0 bioconda::blast=2.16.0"
    publishDir "${params.outdir}/reference_proteomes", mode: 'copy', pattern: '*.faa'

    input:
    val meta

    output:
    path "${meta.name}.faa", emit: faa
    path "${meta.name}_db.*", emit: db

    script:
    """
    if [ -f "${params.local_ref_proteomes_dir}/${meta.name}.faa" ]; then
        echo "Usando cache local: ${params.local_ref_proteomes_dir}/${meta.name}.faa" >&2
        cp "${params.local_ref_proteomes_dir}/${meta.name}.faa" ./${meta.name}.faa
    else
        datasets download genome accession ${meta.accession} \\
            --include protein --filename ${meta.name}.zip \\
          || { echo "ERRO: download NCBI falhou para ${meta.accession} (${meta.taxon})," >&2; \\
               echo "provavelmente por bloqueio de rede do servidor." >&2; \\
               echo "Baixe manualmente no Windows: NCBI Datasets -> ${meta.accession} -> Download -> Protein (FASTA)," >&2; \\
               echo "salve como ${meta.name}.faa e transfira via scp para" >&2; \\
               echo "${params.local_ref_proteomes_dir}/${meta.name}.faa antes de re-rodar." >&2; \\
               exit 1; }
        unzip -p ${meta.name}.zip 'ncbi_dataset/data/*/protein.faa' > ${meta.name}.faa
    fi

    makeblastdb -in ${meta.name}.faa -dbtype prot -parse_seqids -out ${meta.name}_db
    """
}
