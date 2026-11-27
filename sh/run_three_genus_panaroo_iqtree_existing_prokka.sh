#!/usr/bin/env bash

set -euo pipefail

THREADS=${THREADS:-8}

PROJECT_DIR=$(pwd)

PROKKA_DIR="${PROJECT_DIR}/01_prokka"
RESULT_DIR="${PROJECT_DIR}/three_genus_results"

mkdir -p "${RESULT_DIR}"

echo "=========================================="
echo "Three genus Panaroo + IQ-TREE pipeline"
echo "Using existing Prokka results"
echo "Threads: ${THREADS}"
echo "=========================================="


for x in panaroo mafft iqtree2; do
    conda run -n genome_compare bash -c "command -v ${x}" >/dev/null 2>&1 || {
        echo "Missing ${x} in genome_compare environment"
        exit 1
    }
done


prepare_gffs(){

    GROUP=$1

    mkdir -p "${RESULT_DIR}/${GROUP}/gffs"

    echo "Preparing ${GROUP} GFF files"

    case ${GROUP} in

    Genus1)
        PATTERN="1CP|A1_final"
        ;;

    Genus2)
        PATTERN="3CP|A3_final"
        ;;

    Genus3)
        PATTERN="^G1_final|^GCP|^GLR134387"
        ;;

    esac


    find "${PROKKA_DIR}" -mindepth 2 -maxdepth 2 \
        -name "*.gff" | \
        while read gff
        do
            name=$(basename "${gff}" .gff)

            if echo "${name}" | grep -Eq "${PATTERN}"
            then
                cp "${gff}" "${RESULT_DIR}/${GROUP}/gffs/"
            fi
        done


    echo "${GROUP} GFF number:"
    ls "${RESULT_DIR}/${GROUP}/gffs/"*.gff | wc -l
}


run_panaroo(){

    GROUP=$1

    WORK="${RESULT_DIR}/${GROUP}"

    mkdir -p "${WORK}/panaroo"


    echo "=========================================="
    echo "Running Panaroo ${GROUP}"
    echo "=========================================="


    conda run -n genome_compare panaroo \
        -i ${WORK}/gffs/*.gff \
        -o ${WORK}/panaroo \
        --clean-mode moderate \
        --remove-invalid-genes \
        -a core \
        --aligner mafft \
        --codons \
        --core_threshold 0.95 \
        -t ${THREADS}


}


run_iqtree(){

    GROUP=$1

    WORK="${RESULT_DIR}/${GROUP}"

    ALIGN=$(find "${WORK}/panaroo" \
        -name "*core_gene_alignment*.aln" | head -1)


    if [ ! -s "${ALIGN}" ]
    then
        echo "No core alignment found for ${GROUP}"
        exit 1
    fi


    mkdir -p "${WORK}/phylogeny"

    cp "${ALIGN}" \
       "${WORK}/phylogeny/core_alignment.aln"


    echo "Running IQ-TREE ${GROUP}"


    conda run -n genome_compare iqtree2 \
        -s "${WORK}/phylogeny/core_alignment.aln" \
        -m MFP \
        -B 1000 \
        -alrt 1000 \
        -T ${THREADS} \
        --prefix "${WORK}/phylogeny/${GROUP}"

}


for GROUP in Genus1 Genus2 Genus3
do

    prepare_gffs ${GROUP}

    run_panaroo ${GROUP}

    run_iqtree ${GROUP}

done


echo "=========================================="
echo "Pipeline finished successfully"
echo "=========================================="
