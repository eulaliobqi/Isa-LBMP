process HMMSEARCH_PEBP {
    tag "hmmsearch_${params.pfam_id}"
    label 'process_medium'
    conda "bioconda::hmmer=3.4"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path candidates, stageAs: 'candidate_proteins.fasta'

    output:
    path 'pebp_domtbl.txt', emit: domtbl
    path 'pebp_confirmed.tsv', emit: confirmed_tsv
    path 'candidate_proteins_confirmed.fasta', emit: confirmed_proteins

    script:
    """
    if [ -f "${params.pfam_a_hmm}" ]; then
        hmmfetch "${params.pfam_a_hmm}" ${params.pfam_id} > ${params.pfam_id}.hmm 2>hmmfetch.log \\
          || { hmmfetch --index "${params.pfam_a_hmm}"; \\
               hmmfetch "${params.pfam_a_hmm}" ${params.pfam_id} > ${params.pfam_id}.hmm; }
    else
        echo "AVISO: Pfam-A.hmm não encontrado em ${params.pfam_a_hmm}." >&2
        echo "Tentando download do HMM de ${params.pfam_id} via InterPro (pode falhar por bloqueio de rede do servidor)..." >&2
        wget -q -O ${params.pfam_id}.hmm.gz \\
            "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${params.pfam_id}/?annotation=hmm" \\
          || { echo "ERRO: download falhou. Baixe ${params.pfam_id}.hmm manualmente no Windows" >&2; \\
               echo "(https://www.ebi.ac.uk/interpro/entry/pfam/${params.pfam_id}/) e transfira via scp," >&2; \\
               echo "ou aponte --pfam_a_hmm para um Pfam-A.hmm já existente no servidor." >&2; \\
               exit 1; }
        gunzip -f ${params.pfam_id}.hmm.gz
    fi

    hmmsearch \\
        -E ${params.hmmer_evalue} \\
        --domtblout pebp_domtbl.txt \\
        --cpu ${task.cpus} \\
        --noali \\
        ${params.pfam_id}.hmm \\
        candidate_proteins.fasta > hmmsearch.log

    python3 ${projectDir}/bin/parse_hmmer_domtbl.py \\
        --domtbl pebp_domtbl.txt \\
        --candidates candidate_proteins.fasta \\
        --out-tsv pebp_confirmed.tsv \\
        --out-fasta candidate_proteins_confirmed.fasta
    """
}
