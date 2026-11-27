#!/usr/bin/env bash
set -euo pipefail

############################################################
# Prokka -> Panaroo/MAFFT -> IQ-TREE -> antiSMASH
#
# INPUT:
#   genomes/*.fna   or genomes/*.fa or genomes/*.fasta
#
# PROJECT:
#   comparative_genomics/
#     genomes/
#     01_prokka/
#     02_panaroo/
#     03_phylogeny/
#     04_antismash/
#     05_summary/
#
# Usage:
#   cd ~/comparative_genomics
#   bash run_comparative_genomics.sh
############################################################

THREADS=${THREADS:-8}
CORE_THRESHOLD=${CORE_THRESHOLD:-0.99}
CLEAN_MODE=${CLEAN_MODE:-strict}

PROJECT_DIR="$(pwd)"
GENOME_DIR="${PROJECT_DIR}/genomes"
PROKKA_DIR="${PROJECT_DIR}/01_prokka"
PANAROO_DIR="${PROJECT_DIR}/02_panaroo"
PHYLO_DIR="${PROJECT_DIR}/03_phylogeny"
ANTISMASH_DIR="${PROJECT_DIR}/04_antismash"
SUMMARY_DIR="${PROJECT_DIR}/05_summary"

mkdir -p \
  "$GENOME_DIR" \
  "$PROKKA_DIR" \
  "$PANAROO_DIR" \
  "$PHYLO_DIR" \
  "$ANTISMASH_DIR" \
  "$SUMMARY_DIR"

# ----------------------------------------------------------
# Helper: find IQ-TREE executable inside conda env
# ----------------------------------------------------------
find_iqtree() {
  if conda run -n genome_compare bash -lc 'command -v iqtree3' >/dev/null 2>&1; then
    echo "iqtree3"
  elif conda run -n genome_compare bash -lc 'command -v iqtree2' >/dev/null 2>&1; then
    echo "iqtree2"
  elif conda run -n genome_compare bash -lc 'command -v iqtree' >/dev/null 2>&1; then
    echo "iqtree"
  else
    echo ""
  fi
}

IQTREE_BIN="$(find_iqtree)"

if [ -z "$IQTREE_BIN" ]; then
  echo "[ERROR] IQ-TREE executable not found in genome_compare environment."
  exit 1
fi

# ----------------------------------------------------------
# 0. Verify inputs
# ----------------------------------------------------------
mapfile -t FASTAS < <(
  find "$GENOME_DIR" -maxdepth 1 -type f \
    \( -iname "*.fna" -o -iname "*.fa" -o -iname "*.fasta" \) \
    | sort
)

if [ "${#FASTAS[@]}" -lt 3 ]; then
  echo "[ERROR] Need at least 3 genome FASTA files in: $GENOME_DIR"
  exit 1
fi

echo "[INFO] Found ${#FASTAS[@]} genome FASTA files."

# ----------------------------------------------------------
# 1. Prokka annotation
# ----------------------------------------------------------
echo
echo "============================================================"
echo "STEP 1. Prokka annotation"
echo "============================================================"

for FASTA in "${FASTAS[@]}"; do
  BASE="$(basename "$FASTA")"
  SAMPLE="${BASE%.*}"
  SAFE_SAMPLE="$(echo "$SAMPLE" | sed 's/[^A-Za-z0-9_]/_/g')"

  OUT="${PROKKA_DIR}/${SAFE_SAMPLE}"

  if [ -f "${OUT}/${SAFE_SAMPLE}.gff" ]; then
    echo "[SKIP] Existing Prokka result: $SAFE_SAMPLE"
    continue
  fi

  # locus tags must be simple and unique
  LOCUS="$(echo "$SAFE_SAMPLE" | tr -cd '[:alnum:]' | cut -c1-8 | tr '[:lower:]' '[:upper:]')"
  if [ -z "$LOCUS" ]; then
    LOCUS="GENOME"
  fi

  echo "[RUN] Prokka: $SAFE_SAMPLE"

  conda run -n genome_compare prokka \
    --outdir "$OUT" \
    --prefix "$SAFE_SAMPLE" \
    --locustag "$LOCUS" \
    --kingdom Bacteria \
    --cpus "$THREADS" \
    "$FASTA"
done

# ----------------------------------------------------------
# 2. Collect GFFs
# ----------------------------------------------------------
mapfile -t GFFS < <(
  find "$PROKKA_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.gff" | sort
)

if [ "${#GFFS[@]}" -ne "${#FASTAS[@]}" ]; then
  echo "[ERROR] Number of GFF files does not equal number of input genomes."
  echo "FASTA: ${#FASTAS[@]}"
  echo "GFF:   ${#GFFS[@]}"
  exit 1
fi

printf "%s\n" "${GFFS[@]}" > "${SUMMARY_DIR}/prokka_gff_files.txt"

# ----------------------------------------------------------
# 3. Panaroo pangenome + MAFFT core alignment
# ----------------------------------------------------------
echo
echo "============================================================"
echo "STEP 2. Panaroo pangenome and core alignment"
echo "============================================================"

if [ ! -f "${PANAROO_DIR}/gene_presence_absence.csv" ]; then

  conda run -n genome_compare panaroo \
    -i "${GFFS[@]}" \
    -o "$PANAROO_DIR" \
    --clean-mode "$CLEAN_MODE" \
    --remove-invalid-genes \
    -a core \
    --aligner mafft \
    --codons \
    --core_threshold "$CORE_THRESHOLD" \
    -t "$THREADS"

else
  echo "[SKIP] Existing Panaroo result detected."
fi

# ----------------------------------------------------------
# 4. Select core alignment
# ----------------------------------------------------------
FILTERED_ALIGN="${PANAROO_DIR}/core_gene_alignment_filtered.aln"
RAW_ALIGN="${PANAROO_DIR}/core_gene_alignment.aln"

if [ -s "$FILTERED_ALIGN" ]; then
  CORE_ALIGN="$FILTERED_ALIGN"
  echo "[INFO] Using filtered Panaroo core alignment:"
  echo "       $CORE_ALIGN"
elif [ -s "$RAW_ALIGN" ]; then
  CORE_ALIGN="$RAW_ALIGN"
  echo "[WARN] Filtered alignment absent; using raw core alignment:"
  echo "       $CORE_ALIGN"
else
  echo "[ERROR] No Panaroo core alignment found."
  exit 1
fi

cp "$CORE_ALIGN" "${PHYLO_DIR}/core_alignment_for_IQTREE.aln"

# ----------------------------------------------------------
# 5. Optional explicit MAFFT example
# ----------------------------------------------------------
# Panaroo already performed the core-gene alignments above.
#
# If you later have an unaligned homologous FASTA:
#
# conda run -n genome_compare mafft \
#   --auto \
#   --thread "$THREADS" \
#   homologs.fasta \
#   > homologs.aligned.fasta
#
# DO NOT re-align the already concatenated Panaroo core alignment.

# ----------------------------------------------------------
# 6. IQ-TREE phylogenomics
# ----------------------------------------------------------
echo
echo "============================================================"
echo "STEP 3. IQ-TREE phylogenomics"
echo "============================================================"

TREE_PREFIX="${PHYLO_DIR}/core_genome"

if [ ! -f "${TREE_PREFIX}.treefile" ]; then

  conda run -n genome_compare "$IQTREE_BIN" \
    -s "${PHYLO_DIR}/core_alignment_for_IQTREE.aln" \
    -m MFP \
    -B 1000 \
    -alrt 1000 \
    -T AUTO \
    --prefix "$TREE_PREFIX"

else
  echo "[SKIP] Existing IQ-TREE result detected."
fi

# ----------------------------------------------------------
# 7. antiSMASH for every genome
# ----------------------------------------------------------
echo
echo "============================================================"
echo "STEP 4. antiSMASH"
echo "============================================================"

mapfile -t GBKS < <(
  find "$PROKKA_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.gbk" | sort
)

for GBK in "${GBKS[@]}"; do
  BASE="$(basename "$GBK")"
  SAMPLE="${BASE%.gbk}"
  OUT="${ANTISMASH_DIR}/${SAMPLE}"

  if find "$OUT" -type f -name "*.region*.gbk" 2>/dev/null | grep -q .; then
    echo "[SKIP] Existing antiSMASH result: $SAMPLE"
    continue
  fi

  mkdir -p "$OUT"

  echo "[RUN] antiSMASH: $SAMPLE"

  conda run -n antismash antismash \
    "$GBK" \
    --output-dir "$OUT" \
    --output-basename "$SAMPLE" \
    --cpus "$THREADS" \
    --knownclusterblast
done

# ----------------------------------------------------------
# 8. Software versions
# ----------------------------------------------------------
{
  echo "DATE: $(date)"
  echo
  echo "Prokka:"
  conda run -n genome_compare prokka --version || true
  echo
  echo "Panaroo:"
  conda run -n genome_compare panaroo --version || true
  echo
  echo "MAFFT:"
  conda run -n genome_compare mafft --version || true
  echo
  echo "IQ-TREE:"
  conda run -n genome_compare "$IQTREE_BIN" --version || true
  echo
  echo "antiSMASH:"
  conda run -n antismash antismash --version || true
} > "${SUMMARY_DIR}/software_versions.txt" 2>&1

# ----------------------------------------------------------
# 9. Copy key files for Windows/R
# ----------------------------------------------------------
cp "${PANAROO_DIR}/gene_presence_absence.csv" \
   "${SUMMARY_DIR}/gene_presence_absence.csv"

if [ -f "${PANAROO_DIR}/summary_statistics.txt" ]; then
  cp "${PANAROO_DIR}/summary_statistics.txt" \
     "${SUMMARY_DIR}/panaroo_summary_statistics.txt"
fi

cp "${PHYLO_DIR}/core_alignment_for_IQTREE.aln" \
   "${SUMMARY_DIR}/core_alignment_for_IQTREE.aln"

cp "${TREE_PREFIX}.treefile" \
   "${SUMMARY_DIR}/core_genome.treefile"

# Copy antiSMASH region table list
find "$ANTISMASH_DIR" -type f -name "*.region*.gbk" | sort \
  > "${SUMMARY_DIR}/antismash_region_files.txt"

echo
echo "============================================================"
echo "PIPELINE FINISHED"
echo "============================================================"
echo "Key results:"
echo "  Panaroo:   $PANAROO_DIR"
echo "  Tree:      ${TREE_PREFIX}.treefile"
echo "  antiSMASH: $ANTISMASH_DIR"
echo "  Summary:   $SUMMARY_DIR"
echo
echo "For Windows, copy the entire project folder or at least:"
echo "  02_panaroo/"
echo "  03_phylogeny/"
echo "  04_antismash/"
echo "  05_summary/"
echo "============================================================"
