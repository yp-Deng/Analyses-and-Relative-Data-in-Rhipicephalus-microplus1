#!/usr/bin/env bash
set -euo pipefail

echo "============================================================"
echo "Comparative genomics software installer for Ubuntu / WSL2"
echo "============================================================"

# 1. Basic Ubuntu utilities
sudo apt update
sudo apt install -y \
  wget curl git unzip bzip2 gzip ca-certificates \
  build-essential

# 2. Install Miniforge if conda is absent
if ! command -v conda >/dev/null 2>&1; then
  echo "[INFO] Conda not found. Installing Miniforge..."
  cd "$HOME"
  wget -O Miniforge3.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

  bash Miniforge3.sh -b -p "$HOME/miniforge3"

  source "$HOME/miniforge3/etc/profile.d/conda.sh"
  "$HOME/miniforge3/bin/conda" init bash

  echo "[INFO] Miniforge installed in $HOME/miniforge3"
else
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi

# Make sure conda works in this shell
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

# 3. Recommended Bioconda channel configuration
conda config --add channels bioconda || true
conda config --add channels conda-forge || true
conda config --set channel_priority strict

# 4. Install phylogenomics/pangenome tools in one environment
echo "[INFO] Creating genome_compare environment..."
conda create -y -n genome_compare \
  prokka \
  panaroo \
  mafft \
  iqtree

# 5. Install antiSMASH separately to reduce dependency conflicts
echo "[INFO] Creating antismash environment..."
conda create -y -n antismash antismash

# 6. Verify software
echo
echo "============================================================"
echo "Verification"
echo "============================================================"

echo "[genome_compare]"
conda run -n genome_compare prokka --version || true
conda run -n genome_compare panaroo --version || true
conda run -n genome_compare mafft --version || true

if conda run -n genome_compare iqtree3 --version >/dev/null 2>&1; then
  conda run -n genome_compare iqtree3 --version
elif conda run -n genome_compare iqtree2 --version >/dev/null 2>&1; then
  conda run -n genome_compare iqtree2 --version
else
  conda run -n genome_compare iqtree --version || true
fi

echo
echo "[antismash]"
conda run -n antismash antismash --version || true

echo
echo "============================================================"
echo "Installation finished."
echo "Restart Ubuntu or run:"
echo "  source ~/.bashrc"
echo
echo "Then activate:"
echo "  conda activate genome_compare"
echo "============================================================"
