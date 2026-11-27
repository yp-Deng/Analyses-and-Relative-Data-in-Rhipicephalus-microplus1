#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-8}
ANTISMASH_MODE=${ANTISMASH_MODE:-none}

PROJECT_DIR=$(pwd)

GENOME_DIR="${PROJECT_DIR}/genomes"
PROKKA_DIR="${PROJECT_DIR}/01_prokka"
PANAROO_DIR="${PROJECT_DIR}/02_panaroo"
PHYLO_DIR="${PROJECT_DIR}/03_phylogeny"
ANTISMASH_DIR="${PROJECT_DIR}/04_antismash"
SUMMARY_DIR="${PROJECT_DIR}/05_summary"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$PHYLO_DIR" "$ANTISMASH_DIR" "$SUMMARY_DIR" "$LOG_DIR"

echo "Comparative genomics pipeline v3"

for x in prokka panaroo mafft; do
  conda run -n genome_compare bash -lc "command -v $x" >/dev/null 2>&1 || {
    echo "[ERROR] missing $x"; exit 1;
  }
done

if conda run -n genome_compare bash -lc "command -v iqtree2" >/dev/null 2>&1; then
 IQTREE=iqtree2
else
 IQTREE=iqtree
fi

echo "STEP 1 Prokka check"

mapfile -t FASTAS < <(find "$GENOME_DIR" -maxdepth 1 -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) | sort)

for fasta in "${FASTAS[@]}"; do
 sample=$(basename "$fasta")
 sample=${sample%.*}

 if [[ -s "$PROKKA_DIR/$sample/$sample.gff" &&
       -s "$PROKKA_DIR/$sample/$sample.gbk" &&
       -s "$PROKKA_DIR/$sample/$sample.faa" ]]; then
   echo "[SKIP] Prokka $sample"
 else
   echo "[ERROR] Prokka incomplete: $sample"
   exit 1
 fi
done

echo "STEP 2 Panaroo check"

PANAROO_DNA="$PANAROO_DIR/combined_DNA_CDS.fasta"
PANAROO_GPA="$PANAROO_DIR/gene_presence_absence.csv"

if [[ ! -s "$PANAROO_DNA" || ! -s "$PANAROO_GPA" ]]; then
 echo "[ERROR] Panaroo result incomplete"
 exit 1
fi

echo "[SKIP] Panaroo complete"

ALIGN="$PHYLO_DIR/core_alignment_for_IQTREE.aln"

echo "STEP 3 MAFFT"

if [[ -s "$ALIGN" ]]; then
 echo "[SKIP] MAFFT alignment"
else
 conda run -n genome_compare mafft \
 --auto \
 --thread "$THREADS" \
 "$PANAROO_DNA" \
 > "$ALIGN" 2> "$LOG_DIR/mafft.log"
 echo "[OK] MAFFT finished"
fi

echo "STEP 4 IQ-TREE"

TREE_PREFIX="$PHYLO_DIR/core_genome"

if [[ -s "${TREE_PREFIX}.treefile" ]]; then
 echo "[SKIP] IQ-TREE"
else
 conda run -n genome_compare "$IQTREE" \
 -s "$ALIGN" \
 -m MFP \
 -B 1000 \
 -alrt 1000 \
 -T "$THREADS" \
 --prefix "$TREE_PREFIX" \
 > "$LOG_DIR/iqtree.log" 2>&1
 echo "[OK] IQ-TREE finished"
fi

echo "STEP 5 antiSMASH"

if [[ "$ANTISMASH_MODE" == "none" ]]; then
 echo "[SKIP] antiSMASH"
else
 for gbk in "$PROKKA_DIR"/*/*.gbk; do
   sample=$(basename "$gbk" .gbk)
   out="$ANTISMASH_DIR/$sample"

   if [[ -s "$out/index.html" ]]; then
     echo "[SKIP] antiSMASH $sample"
   else
     rm -rf "$out"
     conda run -n antismash antismash \
       "$gbk" \
       --output-dir "$out" \
       --cpus "$THREADS" \
       > "$LOG_DIR/antismash_${sample}.log" 2>&1
   fi
 done
fi

cp "$PANAROO_GPA" "$SUMMARY_DIR/gene_presence_absence.csv"
cp "$ALIGN" "$SUMMARY_DIR/core_alignment_for_IQTREE.aln"
cp "${TREE_PREFIX}.treefile" "$SUMMARY_DIR/core_genome.treefile"

echo "Completed $(date)" > "$SUMMARY_DIR/PIPELINE_COMPLETE.txt"

echo "Pipeline completed successfully"
