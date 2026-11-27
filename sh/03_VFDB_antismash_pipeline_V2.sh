#!/usr/bin/env bash

set -euo pipefail

THREADS=${THREADS:-8}
BASE=$(pwd)

PAN_DIR="${BASE}/three_genus_results"
GENOME_DIR="${BASE}/genomes"
OUT="${BASE}/virulence_BGC_results_v2"

mkdir -p ${OUT}/VFDB ${OUT}/VF_summary ${OUT}/antiSMASH ${OUT}/BGC_summary

echo "VFDB + antiSMASH pipeline V2 started"

if [ ! -f VFDB.dmnd ]; then
 echo "ERROR: VFDB.dmnd not found"
 exit 1
fi

command -v diamond >/dev/null || { echo "diamond missing"; exit 1; }
command -v antismash >/dev/null || { echo "antismash missing"; exit 1; }

for G in Genus1 Genus2 Genus3
do
 diamond blastp  -q ${PAN_DIR}/${G}/panaroo/pan_genome_reference.fa  -d VFDB  -e 1e-5  -k 1  --threads ${THREADS}  --outfmt 6  -o ${OUT}/VFDB/${G}_VFDB.tsv

 awk '{print $1}' ${OUT}/VFDB/${G}_VFDB.tsv | sort -u > ${OUT}/VF_summary/${G}_VF_genes.txt
done

python3 <<'PY'
import os
groups=["Genus1","Genus2","Genus3"]
data={}
allgenes=set()
for g in groups:
    f=f"virulence_BGC_results_v2/VF_summary/{g}_VF_genes.txt"
    genes=set(open(f).read().splitlines())
    data[g]=genes
    allgenes.update(genes)

with open("virulence_BGC_results_v2/VF_summary/VF_gene_matrix.csv","w") as out:
    out.write("VF_gene,"+",".join(groups)+"\n")
    for gene in sorted(allgenes):
        out.write(gene+","+",".join(["1" if gene in data[g] else "0" for g in groups])+"\n")
PY

for genome in ${GENOME_DIR}/*.fasta
do
 name=$(basename ${genome} .fasta)
 antismash ${genome} --output-dir ${OUT}/antiSMASH/${name} --cpus ${THREADS}
done

echo "Genome,BGC_number" > ${OUT}/BGC_summary/BGC_counts.csv

for d in ${OUT}/antiSMASH/*
do
 name=$(basename ${d})
 count=$(find ${d} -name "*.gbk" | wc -l)
 echo "${name},${count}" >> ${OUT}/BGC_summary/BGC_counts.csv
done

echo "Pipeline finished successfully"
