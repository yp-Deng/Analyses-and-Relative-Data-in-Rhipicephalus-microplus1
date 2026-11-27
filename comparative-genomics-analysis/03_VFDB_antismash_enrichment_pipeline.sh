#!/usr/bin/env bash

#############################################################
# VFDB + antiSMASH comparative genomics pipeline
#
# Input:
# three_genus_results/
#
# Output:
# virulence_BGC_results/
#
# Requirements:
# diamond
# VFDB.dmnd
# antismash
#
#############################################################


set -euo pipefail


THREADS=8


BASE=$(pwd)


INPUT_DIR="${BASE}/three_genus_results"

OUT="${BASE}/virulence_BGC_results"


mkdir -p ${OUT}

mkdir -p ${OUT}/VFDB
mkdir -p ${OUT}/antiSMASH
mkdir -p ${OUT}/VF_summary
mkdir -p ${OUT}/BGC_summary



echo "========================================"
echo "VFDB + antiSMASH pipeline started"
echo "Threads: ${THREADS}"
echo "========================================"



#############################################################
# Step 1
# VFDB diamond annotation
#############################################################


echo ""
echo "===== Step 1 VFDB annotation ====="
echo ""


for G in Genus1 Genus2 Genus3

do


echo "Running VFDB for ${G}"


diamond blastp \
-q ${INPUT_DIR}/${G}/panaroo/pan_genome_reference.fa \
-d VFDB \
-o ${OUT}/VFDB/${G}_VFDB.tsv \
-e 1e-5 \
-k 1 \
--threads ${THREADS} \
--outfmt 6


done



#############################################################
# Step 2
# Count VF genes
#############################################################


echo ""
echo "===== Step 2 VF gene counting ====="
echo ""



for G in Genus1 Genus2 Genus3

do


awk '
{
print $1
}
' ${OUT}/VFDB/${G}_VFDB.tsv | sort -u \
> ${OUT}/VF_summary/${G}_VF_genes.txt



echo "${G} VF gene number:"
cat ${OUT}/VF_summary/${G}_VF_genes.txt | wc -l


done



#############################################################
# Step 3
# Combine VF gene matrix
#############################################################


echo ""
echo "===== Step 3 VF matrix ====="
echo ""



python3 <<EOF


import os


groups=["Genus1","Genus2","Genus3"]


all_genes=set()


data={}



for g in groups:

    file="virulence_BGC_results/VF_summary/"+g+"_VF_genes.txt"

    genes=set(open(file).read().splitlines())

    data[g]=genes

    all_genes.update(genes)



with open(
"virulence_BGC_results/VF_summary/VF_gene_matrix.csv",
"w"
) as out:


    out.write(
    "VF_gene,"+
    ",".join(groups)+
    "\\n"
    )


    for gene in sorted(all_genes):

        row=[gene]

        for g in groups:

            if gene in data[g]:

                row.append("1")

            else:

                row.append("0")


        out.write(",".join(row)+"\n")

EOF



#############################################################
# Step 4
# antiSMASH
#############################################################


echo ""
echo "===== Step 4 antiSMASH ====="
echo ""



for G in Genus1 Genus2 Genus3

do


echo "Running antiSMASH ${G}"


mkdir -p ${OUT}/antiSMASH/${G}



antismash \
${INPUT_DIR}/${G}/panaroo/pan_genome_reference.fa \
--output-dir ${OUT}/antiSMASH/${G} \
--cpus ${THREADS}



done




#############################################################
# Step 5
# BGC summary
#############################################################


echo ""
echo "===== Step 5 BGC summary ====="
echo ""



echo "Genus,BGC_number" \
> ${OUT}/BGC_summary/BGC_counts.csv



for G in Genus1 Genus2 Genus3

do


COUNT=$(find ${OUT}/antiSMASH/${G} \
-name "*.gbk" | wc -l)


echo "${G},${COUNT}" \
>> ${OUT}/BGC_summary/BGC_counts.csv


done




echo ""
echo "========================================"
echo "Pipeline finished"
echo ""
echo "Results:"
echo "${OUT}"
echo ""
echo "========================================"
