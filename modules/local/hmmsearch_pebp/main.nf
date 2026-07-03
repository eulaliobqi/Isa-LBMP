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
    // Extração via awk (leitura direta do texto plano do Pfam-A.hmm, sem usar
    // hmmfetch/.ssi) -- evita depender/mutar o índice do arquivo compartilhado
    // com o Kerson-paper (hmmfetch --index falhava com "SSI index already
    // exists" em runs repetidos, e não temos como saber se a .ssi existente
    // foi construída com o mesmo formato de accession). Aceita ACC com ou
    // sem sufixo de versão (ex. PF01161 ou PF01161.24).
    """
    if [ -f "${params.pfam_a_hmm}" ]; then
        awk -v acc="${params.pfam_id}" '
            /^HMMER/ { buf = \$0 "\\n"; keep = 0; next }
            /^ACC/   { if (\$2 == acc || \$2 ~ ("^" acc "\\\\.")) keep = 1 }
            { buf = buf \$0 "\\n" }
            /^\\/\\/\$/ { if (keep) printf "%s", buf; buf = ""; keep = 0 }
        ' "${params.pfam_a_hmm}" > ${params.pfam_id}.hmm

        if [ ! -s ${params.pfam_id}.hmm ]; then
            echo "ERRO: ${params.pfam_id} não encontrado em ${params.pfam_a_hmm} (nenhum bloco ACC correspondente)." >&2
            exit 1
        fi
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
