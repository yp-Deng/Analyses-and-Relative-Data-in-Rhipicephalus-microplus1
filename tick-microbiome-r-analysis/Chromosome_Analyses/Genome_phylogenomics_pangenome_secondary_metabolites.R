############################################################
# Three-genome comparative genomics master workflow
#
# Analyses:
#  1. Genome annotation harmonization (Prokka; optional)
#  2. Pangenome analysis with Panaroo
#  3. Core-genome phylogenomics
#       Panaroo orthology
#       MAFFT core-gene alignment
#       IQ-TREE maximum-likelihood tree
#  4. Pangenome statistics and plots in R
#  5. Secondary metabolite/BGC analysis with antiSMASH
#  6. antiSMASH BGC comparison and plots in R
#
# INPUT:
#   genomes/
#      Genome1.fna
#      Genome2.fna
#      Genome3.fna
#
# For a publication-quality phylogeny, ADD close relatives/type
# strains and at least one appropriate outgroup to genomes/.
#
# EXTERNAL SOFTWARE REQUIRED IN PATH:
#   prokka
#   panaroo
#   mafft
#   iqtree3 or iqtree2 or iqtree
#   antismash
#
# Suggested conda setup (run in terminal, NOT inside R):
#
# mamba create -n comparative_genomics \
#   -c conda-forge -c bioconda \
#   prokka panaroo mafft iqtree antismash
#
# conda activate comparative_genomics
#
# NOTE:
# Prokka is used here because Panaroo explicitly supports Prokka-style
# GFF3 input. If you already have mutually consistent Prokka GFF/GBK
# annotations, set run_annotation <- FALSE.
############################################################

rm(list = ls())
set.seed(20260813)

############################
# 0. USER SETTINGS
############################

genome_dir <- "genomes"
annotation_dir <- "01_annotation"
panaroo_dir <- "02_panaroo"
phylo_dir <- "03_phylogenomics"
antismash_dir <- "04_antismash"
summary_dir <- "05_comparative_summary"

dirs <- c(
  annotation_dir,
  panaroo_dir,
  phylo_dir,
  antismash_dir,
  summary_dir
)

invisible(
  lapply(
    dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# Set TRUE if starting from chromosome FASTA files.
run_annotation <- TRUE

# Set TRUE to overwrite previous Prokka results.
force_annotation <- FALSE

# Panaroo cleaning.
panaroo_clean_mode <- "strict"

# For only 3 genomes, "core" = present in all genomes is recommended.
# Panaroo's threshold is a fraction; 0.99 means effectively all 3.
core_threshold <- 0.99

# Number of CPUs.
threads <- max(
  1,
  parallel::detectCores() - 1
)

# IQ-TREE support.
ufboot_replicates <- 1000
alrt_replicates <- 1000

# antiSMASH options.
# Default antiSMASH detection plus KnownClusterBlast for comparison
# with characterized BGCs.
antismash_knownclusterblast <- TRUE

############################
# 1. R PACKAGES
############################

pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "scales",
  "patchwork",
  "ape",
  "ggtree"
)

# ggtree is Bioconductor; install separately if needed.
cran_pkgs <- setdiff(
  pkgs,
  "ggtree"
)

need_cran <- cran_pkgs[
  !vapply(
    cran_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(need_cran) > 0) {
  install.packages(need_cran)
}

if (!requireNamespace("ggtree", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    "ggtree",
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(scales)
  library(patchwork)
  library(ape)
  library(ggtree)
})


############################
# 2. FIGURE STYLE
############################

nature_palette <- c(
  "#E64B35",
  "#4DBBD5",
  "#00A087",
  "#3C5488",
  "#F39B7F",
  "#8491B4",
  "#91D1C2",
  "#7E6148"
)

theme_nature <- function(base_size = 11) {

  theme_classic(base_size = base_size) +

    theme(
      axis.text = element_text(
        colour = "black"
      ),

      axis.title = element_text(
        colour = "black"
      ),

      plot.title = element_text(
        face = "bold",
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.45
      )
    )
}


save_fig <- function(
    p,
    filename,
    width = 7,
    height = 5
) {

  ggsave(
    file.path(
      summary_dir,
      paste0(filename, ".pdf")
    ),
    p,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    file.path(
      summary_dir,
      paste0(filename, ".png")
    ),
    p,
    width = width,
    height = height,
    dpi = 600
  )
}


############################################################
# PART I. SOFTWARE AND INPUT CHECK
############################################################

############################
# 3. CHECK EXTERNAL TOOLS
############################

find_executable <- function(candidates) {

  hits <- Sys.which(
    candidates
  )

  hits <- hits[
    hits != ""
  ]

  if (length(hits) == 0) {
    return(NA_character_)
  }

  unname(
    hits[1]
  )
}


prokka_bin <- find_executable(
  c("prokka")
)

panaroo_bin <- find_executable(
  c("panaroo")
)

mafft_bin <- find_executable(
  c("mafft")
)

iqtree_bin <- find_executable(
  c(
    "iqtree3",
    "iqtree2",
    "iqtree"
  )
)

antismash_bin <- find_executable(
  c("antismash")
)


software_status <- data.frame(
  Software = c(
    "Prokka",
    "Panaroo",
    "MAFFT",
    "IQ-TREE",
    "antiSMASH"
  ),

  Executable = c(
    prokka_bin,
    panaroo_bin,
    mafft_bin,
    iqtree_bin,
    antismash_bin
  ),

  Available = !is.na(
    c(
      prokka_bin,
      panaroo_bin,
      mafft_bin,
      iqtree_bin,
      antismash_bin
    )
  )
)

write.csv(
  software_status,
  file.path(
    summary_dir,
    "00_software_status.csv"
  ),
  row.names = FALSE
)

print(
  software_status
)

required_now <- c(
  panaroo_bin,
  mafft_bin,
  iqtree_bin,
  antismash_bin
)

if (run_annotation) {
  required_now <- c(
    required_now,
    prokka_bin
  )
}

if (any(is.na(required_now))) {

  stop(
    "\nOne or more required external programs are not in PATH.\n",
    "See 00_software_status.csv and install/activate the conda environment.\n"
  )
}


############################
# 4. FIND CHROMOSOME FASTA FILES
############################

if (!dir.exists(genome_dir)) {
  dir.create(
    genome_dir,
    recursive = TRUE
  )

  stop(
    "Created folder 'genomes'. Put chromosome FASTA files there and rerun."
  )
}


genome_files <- list.files(
  genome_dir,
  pattern = "\\.(fna|fa|fasta)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(genome_files) < 3) {

  stop(
    "At least 3 chromosome FASTA files are required. Found: ",
    length(genome_files)
  )
}


genome_names <- tools::file_path_sans_ext(
  basename(
    genome_files
  )
)

# Sanitize names for downstream command-line software.
safe_genome_names <- make.names(
  genome_names,
  unique = TRUE
)

safe_genome_names <- gsub(
  "\\.",
  "_",
  safe_genome_names
)

genome_map <- data.frame(
  Original = genome_names,
  Genome = safe_genome_names,
  FASTA = normalizePath(
    genome_files,
    winslash = "/"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  genome_map,
  file.path(
    summary_dir,
    "00_genome_name_map.csv"
  ),
  row.names = FALSE
)


############################################################
# PART II. CONSISTENT ANNOTATION
############################################################

############################
# 5. RUN PROKKA
############################

if (run_annotation) {

  for (i in seq_len(nrow(genome_map))) {

    genome <- genome_map$Genome[i]
    fasta <- genome_map$FASTA[i]

    out <- file.path(
      annotation_dir,
      genome
    )

    gff_expected <- file.path(
      out,
      paste0(
        genome,
        ".gff"
      )
    )

    if (
      file.exists(
        gff_expected
      ) &&
      !force_annotation
    ) {

      message(
        "Skipping annotation already present: ",
        genome
      )

      next
    }

    args <- c(
      "--outdir",
      out,

      "--prefix",
      genome,

      "--locustag",
      toupper(
        substr(
          gsub(
            "[^A-Za-z0-9]",
            "",
            genome
          ),
          1,
          8
        )
      ),

      "--kingdom",
      "Bacteria",

      "--cpus",
      threads,

      "--force",

      fasta
    )

    message(
      "\nRunning Prokka for ",
      genome
    )

    status <- system2(
      prokka_bin,
      args = args
    )

    if (status != 0) {
      stop(
        "Prokka failed for ",
        genome
      )
    }
  }
}


############################
# 6. COLLECT GFF / GBK / FAA
############################

gff_files <- file.path(
  annotation_dir,
  genome_map$Genome,
  paste0(
    genome_map$Genome,
    ".gff"
  )
)

gbk_files <- file.path(
  annotation_dir,
  genome_map$Genome,
  paste0(
    genome_map$Genome,
    ".gbk"
  )
)

faa_files <- file.path(
  annotation_dir,
  genome_map$Genome,
  paste0(
    genome_map$Genome,
    ".faa"
  )
)

if (!all(file.exists(gff_files))) {

  stop(
    "Not all GFF files were found. ",
    "If run_annotation=FALSE, place Prokka outputs in 01_annotation/<Genome>/."
  )
}

if (!all(file.exists(gbk_files))) {
  stop(
    "Not all GBK files were found."
  )
}

if (!all(file.exists(faa_files))) {
  stop(
    "Not all FAA files were found."
  )
}


############################################################
# PART III. PANGENOME WITH PANAROO
############################################################

############################
# 7. RUN PANAROO
############################

panaroo_core_alignment <- file.path(
  panaroo_dir,
  "core_gene_alignment.aln"
)

if (!file.exists(panaroo_core_alignment)) {

  args <- c(
    "-i",
    gff_files,

    "-o",
    panaroo_dir,

    "--clean-mode",
    panaroo_clean_mode,

    "-a",
    "core",

    "--aligner",
    "mafft",

    "--core_threshold",
    core_threshold,

    "--codons",

    "-t",
    threads
  )

  message(
    "\nRunning Panaroo pangenome + core-gene alignment..."
  )

  status <- system2(
    panaroo_bin,
    args = args
  )

  if (status != 0) {
    stop(
      "Panaroo failed."
    )
  }
}


############################################################
# PART IV. EXPLICIT SEQUENCE ALIGNMENT CODE
############################################################

############################
# 8. ALIGN AN INDIVIDUAL HOMOLOGOUS GENE SET WITH MAFFT
############################
# Panaroo already uses MAFFT above for core genes.
# This section provides a reusable explicit alignment function for
# ANY FASTA containing homologous nucleotide or protein sequences.
#
# Example:
# align_with_mafft(
#   input_fasta = "my_gene_sequences.fasta",
#   output_fasta = "my_gene_sequences.aligned.fasta"
# )

align_with_mafft <- function(
    input_fasta,
    output_fasta,
    strategy = "auto",
    threads = threads
) {

  if (!file.exists(input_fasta)) {
    stop(
      "Input FASTA not found: ",
      input_fasta
    )
  }

  args <- c(
    paste0(
      "--",
      strategy
    ),

    "--thread",
    threads,

    input_fasta
  )

  stderr_file <- paste0(
    output_fasta,
    ".mafft.log"
  )

  status <- system2(
    mafft_bin,
    args = args,
    stdout = output_fasta,
    stderr = stderr_file
  )

  if (status != 0) {

    stop(
      "MAFFT failed for ",
      input_fasta,
      ". See ",
      stderr_file
    )
  }

  invisible(
    output_fasta
  )
}


############################
# 9. OPTIONAL: REALIGN CONCATENATED CORE ALIGNMENT
############################
# Usually NOT necessary because Panaroo already aligned individual
# core genes before concatenation.
#
# Keep the Panaroo-generated core_gene_alignment.aln for tree inference.

file.copy(
  panaroo_core_alignment,
  file.path(
    phylo_dir,
    "core_gene_alignment.aln"
  ),
  overwrite = TRUE
)

core_alignment <- file.path(
  phylo_dir,
  "core_gene_alignment.aln"
)


############################################################
# PART V. CORE-GENOME PHYLOGENOMICS
############################################################

############################
# 10. IQ-TREE ML PHYLOGENY
############################
#
# -m MFP = ModelFinder Plus
# -B 1000 = ultrafast bootstrap
# -alrt 1000 = SH-like aLRT
# -T AUTO = automatic CPU allocation
#
# Because Panaroo --codons produced nucleotide codon alignments,
# IQ-TREE will infer the best nucleotide model automatically.

tree_prefix <- file.path(
  phylo_dir,
  "core_genome_IQTREE"
)

tree_file <- paste0(
  tree_prefix,
  ".treefile"
)

if (!file.exists(tree_file)) {

  args <- c(
    "-s",
    core_alignment,

    "-m",
    "MFP",

    "-B",
    ufboot_replicates,

    "-alrt",
    alrt_replicates,

    "-T",
    "AUTO",

    "--prefix",
    tree_prefix
  )

  message(
    "\nRunning IQ-TREE core-genome phylogeny..."
  )

  status <- system2(
    iqtree_bin,
    args = args
  )

  if (status != 0) {
    stop(
      "IQ-TREE failed."
    )
  }
}


############################
# 11. TREE PLOT IN R
############################

if (file.exists(tree_file)) {

  tr <- ape::read.tree(
    tree_file
  )

  p_tree <- ggtree::ggtree(
    tr
  ) +

    ggtree::geom_tiplab(
      size = 3.5,
      fontface = "italic"
    ) +

    ggtree::geom_text2(
      aes(
        subset = !isTip,
        label = label
      ),
      hjust = -0.15,
      size = 3
    ) +

    ggplot2::labs(
      title =
        "Core-genome maximum-likelihood phylogeny"
    ) +

    theme_tree2()

  ggsave(
    file.path(
      summary_dir,
      "01_core_genome_phylogeny.pdf"
    ),
    p_tree,
    width = 7,
    height = 5,
    device = cairo_pdf
  )

  ggsave(
    file.path(
      summary_dir,
      "01_core_genome_phylogeny.png"
    ),
    p_tree,
    width = 7,
    height = 5,
    dpi = 600
  )
}


############################################################
# PART VI. PANGENOME ANALYSIS IN R
############################################################

############################
# 12. READ PANAROO GENE PRESENCE/ABSENCE
############################

pa_file <- file.path(
  panaroo_dir,
  "gene_presence_absence.csv"
)

if (!file.exists(pa_file)) {
  stop(
    "Panaroo gene_presence_absence.csv not found."
  )
}

pa <- readr::read_csv(
  pa_file,
  show_col_types = FALSE
)

write.csv(
  pa,
  file.path(
    summary_dir,
    "02_pangenome_gene_presence_absence_copy.csv"
  ),
  row.names = FALSE
)


############################
# 13. IDENTIFY GENOME COLUMNS
############################

# Panaroo metadata columns usually end before genome-specific columns.
# Use known genome names to find the actual sample columns.

genome_cols <- intersect(
  genome_map$Genome,
  names(
    pa
  )
)

if (length(genome_cols) != nrow(genome_map)) {

  # Try partial matching if Panaroo sanitized names.
  candidate_cols <- names(
    pa
  )[
    names(pa) %in%
      basename(
        tools::file_path_sans_ext(
          gff_files
        )
      )
  ]

  genome_cols <- candidate_cols
}

if (length(genome_cols) < 3) {

  stop(
    "Could not identify genome columns in Panaroo gene_presence_absence.csv.\n",
    "Detected: ",
    paste(
      genome_cols,
      collapse = ", "
    )
  )
}


############################
# 14. BINARY PRESENCE MATRIX
############################

presence <- pa %>%
  select(
    all_of(
      genome_cols
    )
  ) %>%
  mutate(
    across(
      everything(),
      ~ as.integer(
        !is.na(.x) &
          trimws(
            as.character(.x)
          ) != ""
      )
    )
  ) %>%
  as.data.frame()

rownames(
  presence
) <- if (
  "Gene" %in%
    names(pa)
) {
  make.unique(
    as.character(
      pa$Gene
    )
  )
} else {
  paste0(
    "Cluster_",
    seq_len(
      nrow(pa)
    )
  )
}


############################
# 15. PANGENOME CATEGORIES
############################

n_genomes <- ncol(
  presence
)

n_present <- rowSums(
  presence
)

pan_category <- case_when(

  n_present ==
    n_genomes ~
    "Core",

  n_present ==
    n_genomes - 1 ~
    "Shell",

  TRUE ~
    "Cloud"
)

# For exactly 3 genomes:
# Core  = 3/3
# Shell = 2/3
# Cloud = 1/3

pan_summary <- data.frame(
  Cluster = rownames(
    presence
  ),

  N_genomes =
    n_present,

  Category =
    pan_category
)

if ("Gene" %in% names(pa)) {
  pan_summary$Gene <-
    pa$Gene
}

if ("Annotation" %in% names(pa)) {
  pan_summary$Annotation <-
    pa$Annotation
}

write.csv(
  pan_summary,
  file.path(
    summary_dir,
    "02_pangenome_cluster_categories.csv"
  ),
  row.names = FALSE
)


############################
# 16. PANGENOME CATEGORY COUNTS
############################

pan_counts <- pan_summary %>%
  count(
    Category,
    name = "Gene_clusters"
  )

write.csv(
  pan_counts,
  file.path(
    summary_dir,
    "02_pangenome_category_counts.csv"
  ),
  row.names = FALSE
)


p_pan_counts <- ggplot(
  pan_counts,
  aes(
    x = Category,
    y = Gene_clusters,
    fill = Category
  )
) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label = Gene_clusters
    ),
    vjust = -0.35,
    size = 4
  ) +

  scale_fill_manual(
    values = c(
      "Core" = "#E64B35",
      "Shell" = "#4DBBD5",
      "Cloud" = "#00A087"
    )
  ) +

  labs(
    x = NULL,
    y = "Gene clusters",
    title = "Pangenome composition"
  ) +

  theme_nature() +

  theme(
    legend.position = "none"
  )

save_fig(
  p_pan_counts,
  "02_pangenome_composition",
  width = 5.5,
  height = 4.6
)


############################
# 17. GENOME-SPECIFIC ACCESSORY GENES
############################

specific_list <- lapply(
  genome_cols,
  function(g) {

    others <- setdiff(
      genome_cols,
      g
    )

    idx <-
      presence[
        ,
        g
      ] == 1 &
      rowSums(
        presence[
          ,
          others,
          drop = FALSE
        ]
      ) == 0

    data.frame(
      Genome = g,
      Cluster = rownames(
        presence
      )[
        idx
      ],
      stringsAsFactors = FALSE
    )
  }
)

genome_specific <- bind_rows(
  specific_list
)

# Join annotations if possible.
annotation_lookup <- pan_summary %>%
  select(
    any_of(
      c(
        "Cluster",
        "Gene",
        "Annotation"
      )
    )
  )

genome_specific <- genome_specific %>%
  left_join(
    annotation_lookup,
    by = "Cluster"
  )

write.csv(
  genome_specific,
  file.path(
    summary_dir,
    "02_genome_specific_accessory_clusters.csv"
  ),
  row.names = FALSE
)


specific_counts <- genome_specific %>%
  count(
    Genome,
    name = "Specific_clusters"
  )

p_specific <- ggplot(
  specific_counts,
  aes(
    x = Genome,
    y = Specific_clusters,
    fill = Genome
  )
) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label = Specific_clusters
    ),
    vjust = -0.35
  ) +

  scale_fill_manual(
    values =
      rep(
        nature_palette,
        length.out =
          length(
            unique(
              specific_counts$Genome
            )
          )
      )
  ) +

  labs(
    x = NULL,
    y = "Genome-specific gene clusters",
    title = "Genome-specific accessory repertoire"
  ) +

  theme_nature() +

  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

save_fig(
  p_specific,
  "02_genome_specific_accessory_genes",
  width = max(
    5.5,
    nrow(
      specific_counts
    ) * 0.75
  ),
  height = 4.8
)


############################
# 18. PANGENOME ACCUMULATION
############################
# With only 3 genomes this is DESCRIPTIVE only.
# It should not be used to claim an open/closed pangenome.
#
# Compute all permutations when n <= 8.

pan_accumulate <- function(
    presence_matrix
) {

  genomes <- colnames(
    presence_matrix
  )

  n <- length(
    genomes
  )

  if (n <= 8) {

    # All permutations.
    perms <- gtools::permutations(
      n = n,
      r = n,
      v = genomes
    )

  } else {

    # Random permutations for larger data sets.
    perms <- replicate(
      1000,
      sample(
        genomes
      ),
      simplify = FALSE
    )

    perms <- do.call(
      rbind,
      perms
    )
  }

  out <- vector(
    "list",
    nrow(
      perms
    )
  )

  for (i in seq_len(nrow(perms))) {

    ord <- perms[
      i,
    ]

    rows <- vector(
      "list",
      n
    )

    for (k in seq_len(n)) {

      sub <- presence_matrix[
        ,
        ord[
          seq_len(k)
        ],
        drop = FALSE
      ]

      rows[[k]] <- data.frame(
        Permutation = i,
        N_genomes = k,
        Pangenome =
          sum(
            rowSums(
              sub
            ) > 0
          ),
        Core =
          sum(
            rowSums(
              sub
            ) == k
          )
      )
    }

    out[[i]] <- bind_rows(
      rows
    )
  }

  bind_rows(
    out
  )
}


# gtools is used only here.
if (!requireNamespace("gtools", quietly = TRUE)) {
  install.packages(
    "gtools"
  )
}

accum <- pan_accumulate(
  as.matrix(
    presence
  )
)

write.csv(
  accum,
  file.path(
    summary_dir,
    "02_pangenome_accumulation_permutations.csv"
  ),
  row.names = FALSE
)


accum_summary <- accum %>%
  group_by(
    N_genomes
  ) %>%
  summarise(
    Pan_mean =
      mean(
        Pangenome
      ),

    Pan_min =
      min(
        Pangenome
      ),

    Pan_max =
      max(
        Pangenome
      ),

    Core_mean =
      mean(
        Core
      ),

    Core_min =
      min(
        Core
      ),

    Core_max =
      max(
        Core
      ),

    .groups = "drop"
  )


p_accum <- ggplot(
  accum_summary,
  aes(
    x = N_genomes
  )
) +

  geom_ribbon(
    aes(
      ymin = Pan_min,
      ymax = Pan_max
    ),
    alpha = 0.15
  ) +

  geom_line(
    aes(
      y = Pan_mean,
      linetype = "Pangenome"
    ),
    linewidth = 1
  ) +

  geom_point(
    aes(
      y = Pan_mean
    ),
    size = 2.5
  ) +

  geom_ribbon(
    aes(
      ymin = Core_min,
      ymax = Core_max
    ),
    alpha = 0.15
  ) +

  geom_line(
    aes(
      y = Core_mean,
      linetype = "Core genome"
    ),
    linewidth = 1
  ) +

  geom_point(
    aes(
      y = Core_mean
    ),
    size = 2.5
  ) +

  scale_linetype_manual(
    values = c(
      "Pangenome" = "solid",
      "Core genome" = "dashed"
    )
  ) +

  scale_x_continuous(
    breaks = seq_len(
      n_genomes
    )
  ) +

  labs(
    x = "Number of genomes",
    y = "Gene clusters",
    linetype = NULL,
    title = "Pangenome accumulation",
    subtitle =
      "Descriptive only when based on three genomes"
  ) +

  theme_nature()

save_fig(
  p_accum,
  "02_pangenome_accumulation",
  width = 6.5,
  height = 4.8
)


############################################################
# PART VII. SECONDARY METABOLITE ANALYSIS WITH ANTISMASH
############################################################

############################
# 19. RUN ANTISMASH
############################
# antiSMASH ideally uses annotated GenBank input.
#
# Each genome gets its own output directory.

for (i in seq_along(gbk_files)) {

  genome <- genome_map$Genome[i]
  gbk <- normalizePath(
    gbk_files[i],
    winslash = "/"
  )

  out <- file.path(
    antismash_dir,
    genome
  )

  json_expected <- file.path(
    out,
    paste0(
      genome,
      ".json"
    )
  )

  # antiSMASH output basename behavior can differ,
  # so also detect any JSON.
  existing_json <- list.files(
    out,
    pattern = "\\.json$",
    full.names = TRUE
  )

  if (length(existing_json) > 0) {

    message(
      "Skipping existing antiSMASH result: ",
      genome
    )

    next
  }

  args <- c(
    gbk,

    "--output-dir",
    out,

    "--output-basename",
    genome,

    "--cpus",
    threads
  )

  if (antismash_knownclusterblast) {

    args <- c(
      args,
      "--knownclusterblast"
    )
  }

  message(
    "\nRunning antiSMASH for ",
    genome
  )

  status <- system2(
    antismash_bin,
    args = args
  )

  if (status != 0) {
    stop(
      "antiSMASH failed for ",
      genome
    )
  }
}


############################################################
# PART VIII. PARSE ANTISMASH REGION GENBANK FILES
############################################################

############################
# 20. REGION GBK PARSER
############################
# antiSMASH produces per-region GenBank files such as:
#   Genome.region001.gbk
#
# This parser reads the top-level "region" feature and extracts:
#   location
#   product(s)
#
# It does NOT infer the final chemical compound.
# antiSMASH BGC boundaries are predictions.

extract_quoted_qualifier <- function(
    block,
    qualifier
) {

  pat <- paste0(
    "/",
    qualifier,
    "=\""
  )

  idx <- grep(
    pat,
    block,
    fixed = TRUE
  )

  if (length(idx) == 0) {
    return(
      NA_character_
    )
  }

  line <- block[
    idx[1]
  ]

  value <- sub(
    paste0(
      ".*",
      pat
    ),
    "",
    line,
    fixed = TRUE
  )

  # Handle multi-line quoted qualifiers.
  j <- idx[1]

  while (
    !grepl(
      "\"$",
      value
    ) &&
      j < length(block)
  ) {

    j <- j + 1

    value <- paste0(
      value,
      " ",
      trimws(
        block[j]
      )
    )
  }

  value <- sub(
    "\"$",
    "",
    value
  )

  value
}


parse_antismash_region <- function(
    file,
    genome
) {

  x <- readLines(
    file,
    warn = FALSE
  )

  # Find exact region feature line.
  region_idx <- grep(
    "^     region\\s+",
    x
  )

  if (length(region_idx) == 0) {

    return(
      data.frame(
        Genome = genome,
        Region_file = basename(file),
        Region_location = NA_character_,
        BGC_type = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  i <- region_idx[1]

  location <- trimws(
    sub(
      "^     region\\s+",
      "",
      x[i]
    )
  )

  # Find next feature beginning with five spaces + non-space,
  # excluding continuation/qualifier lines.
  next_feature <- which(
    seq_along(x) > i &
      grepl(
        "^     [A-Za-z_][A-Za-z0-9_]*\\s+",
        x
      )
  )

  if (length(next_feature) > 0) {
    end <- min(
      next_feature
    ) - 1
  } else {
    end <- length(
      x
    )
  }

  block <- x[
    i:end
  ]

  product <- extract_quoted_qualifier(
    block,
    "product"
  )

  # Some antiSMASH versions may store products in /rules.
  rules <- extract_quoted_qualifier(
    block,
    "rules"
  )

  if (
    is.na(product) ||
    product == ""
  ) {
    product <- rules
  }

  data.frame(
    Genome = genome,
    Region_file = basename(file),
    Region_location = location,
    BGC_type = product,
    stringsAsFactors = FALSE
  )
}


############################
# 21. PARSE ALL REGIONS
############################

bgc_list <- list()

counter <- 1

for (genome in genome_map$Genome) {

  out <- file.path(
    antismash_dir,
    genome
  )

  region_files <- list.files(
    out,
    pattern = "region[0-9]+\\.gbk$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(region_files) == 0) {

    bgc_list[[counter]] <- data.frame(
      Genome = genome,
      Region_file = NA_character_,
      Region_location = NA_character_,
      BGC_type = NA_character_,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1

    next
  }

  for (rf in region_files) {

    bgc_list[[counter]] <-
      parse_antismash_region(
        rf,
        genome
      )

    counter <- counter + 1
  }
}


bgc_regions <- bind_rows(
  bgc_list
)

write.csv(
  bgc_regions,
  file.path(
    summary_dir,
    "03_antismash_BGC_regions.csv"
  ),
  row.names = FALSE
)


############################
# 22. BGC COUNTS PER GENOME
############################

bgc_counts <- bgc_regions %>%
  filter(
    !is.na(
      Region_file
    )
  ) %>%
  count(
    Genome,
    name = "BGC_regions"
  ) %>%
  tidyr::complete(
    Genome =
      genome_map$Genome,
    fill = list(
      BGC_regions = 0
    )
  )


write.csv(
  bgc_counts,
  file.path(
    summary_dir,
    "03_antismash_BGC_counts.csv"
  ),
  row.names = FALSE
)


p_bgc_total <- ggplot(
  bgc_counts,
  aes(
    x = Genome,
    y = BGC_regions,
    fill = Genome
  )
) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label = BGC_regions
    ),
    vjust = -0.35
  ) +

  scale_fill_manual(
    values = rep(
      nature_palette,
      length.out =
        nrow(
          bgc_counts
        )
    )
  ) +

  labs(
    x = NULL,
    y = "Predicted BGC regions",
    title = "Secondary-metabolite biosynthetic potential"
  ) +

  theme_nature() +

  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

save_fig(
  p_bgc_total,
  "03_antismash_total_BGCs",
  width = max(
    5.5,
    nrow(
      bgc_counts
    ) * 0.75
  ),
  height = 4.8
)


############################
# 23. SPLIT HYBRID/MULTIPLE BGC TYPE LABELS
############################
# antiSMASH region product strings can contain multiple product classes.
# We retain both:
#   - complete region-level product string
#   - individual components split on commas/semicolons

bgc_types <- bgc_regions %>%
  filter(
    !is.na(
      BGC_type
    ),
    BGC_type != ""
  ) %>%
  separate_rows(
    BGC_type,
    sep = "[,;]"
  ) %>%
  mutate(
    BGC_type =
      trimws(
        BGC_type
      )
  ) %>%
  filter(
    BGC_type != ""
  )


bgc_type_counts <- bgc_types %>%
  count(
    Genome,
    BGC_type,
    name = "N_regions"
  )

write.csv(
  bgc_type_counts,
  file.path(
    summary_dir,
    "03_antismash_BGC_type_counts.csv"
  ),
  row.names = FALSE
)


############################
# 24. BGC TYPE HEATMAP
############################

if (nrow(bgc_type_counts) > 0) {

  p_bgc_heat <- ggplot(
    bgc_type_counts,
    aes(
      x = Genome,
      y = BGC_type,
      fill = N_regions
    )
  ) +

    geom_tile(
      colour = "white",
      linewidth = 0.4
    ) +

    geom_text(
      aes(
        label = N_regions
      ),
      size = 3.5
    ) +

    scale_fill_gradient(
      low = "white",
      high = "#E64B35"
    ) +

    labs(
      x = NULL,
      y = NULL,
      fill = "Regions",
      title = "Predicted biosynthetic gene-cluster types"
    ) +

    theme_nature() +

    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      )
    )

  save_fig(
    p_bgc_heat,
    "03_antismash_BGC_type_heatmap",
    width = max(
      6,
      length(
        unique(
          bgc_type_counts$Genome
        )
      ) * 0.8
    ),
    height = max(
      4.8,
      length(
        unique(
          bgc_type_counts$BGC_type
        )
      ) * 0.38
    )
  )
}


############################################################
# PART IX. INTEGRATED GENOME COMPARISON SUMMARY
############################################################

############################
# 25. PROKKA CDS COUNTS
############################

prokka_stats <- purrr::map_dfr(
  seq_along(
    genome_map$Genome
  ),
  function(i) {

    genome <- genome_map$Genome[i]

    tsv <- file.path(
      annotation_dir,
      genome,
      paste0(
        genome,
        ".tsv"
      )
    )

    if (!file.exists(tsv)) {

      return(
        data.frame(
          Genome = genome,
          Total_CDS = NA_integer_
        )
      )
    }

    x <- readr::read_tsv(
      tsv,
      show_col_types = FALSE
    )

    if (
      "ftype" %in%
        names(
          x
        )
    ) {

      n_cds <- sum(
        x$ftype == "CDS",
        na.rm = TRUE
      )

    } else {

      n_cds <- NA_integer_
    }

    data.frame(
      Genome = genome,
      Total_CDS = n_cds
    )
  }
)


############################
# 26. CORE / ACCESSORY COUNTS BY GENOME
############################

core_idx <- n_present ==
  n_genomes

shell_idx <- n_present ==
  n_genomes - 1

cloud_idx <- n_present <
  n_genomes - 1


per_genome_pan <- purrr::map_dfr(
  genome_cols,
  function(g) {

    data.frame(
      Genome = g,

      Core_clusters =
        sum(
          presence[
            core_idx,
            g
          ] == 1
        ),

      Shell_clusters =
        sum(
          presence[
            shell_idx,
            g
          ] == 1
        ),

      Cloud_clusters =
        sum(
          presence[
            cloud_idx,
            g
          ] == 1
        )
    )
  }
)


############################
# 27. FINAL SUMMARY TABLE
############################

final_summary <- genome_map %>%
  select(
    Genome
  ) %>%
  left_join(
    prokka_stats,
    by = "Genome"
  ) %>%
  left_join(
    per_genome_pan,
    by = "Genome"
  ) %>%
  left_join(
    specific_counts,
    by = "Genome"
  ) %>%
  left_join(
    bgc_counts,
    by = "Genome"
  )

write.csv(
  final_summary,
  file.path(
    summary_dir,
    "04_integrated_genome_comparison_summary.csv"
  ),
  row.names = FALSE
)


############################
# 28. SESSION / COMMAND DOCUMENTATION
############################

capture.output(
  sessionInfo(),
  file = file.path(
    summary_dir,
    "99_R_sessionInfo.txt"
  )
)


version_commands <- list(
  prokka = c(
    prokka_bin,
    "--version"
  ),

  panaroo = c(
    panaroo_bin,
    "--version"
  ),

  mafft = c(
    mafft_bin,
    "--version"
  ),

  iqtree = c(
    iqtree_bin,
    "--version"
  ),

  antismash = c(
    antismash_bin,
    "--version"
  )
)


version_output <- purrr::imap_chr(
  version_commands,
  function(cmd, nm) {

    res <- tryCatch(
      system2(
        cmd[1],
        args = cmd[-1],
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) {
        paste(
          "ERROR:",
          e$message
        )
      }
    )

    paste0(
      nm,
      ": ",
      paste(
        res,
        collapse = " "
      )
    )
  }
)

writeLines(
  version_output,
  file.path(
    summary_dir,
    "99_external_software_versions.txt"
  )
)


############################
# 29. ANALYSIS NOTES
############################

notes <- c(

  "Genome phylogenomics:",
  paste0(
    "- Panaroo core threshold = ",
    core_threshold,
    "."
  ),
  "- Panaroo core genes were aligned with MAFFT using codon-aware output.",
  paste0(
    "- IQ-TREE used ModelFinder Plus, ",
    ufboot_replicates,
    " ultrafast bootstrap replicates and ",
    alrt_replicates,
    " SH-aLRT replicates."
  ),

  "",

  "Pangenome:",
  "- With exactly three genomes, Core=3/3, Shell=2/3 and Cloud=1/3.",
  "- Genome-specific clusters are clusters present in one genome and absent from all others.",
  "- Pangenome accumulation with only three genomes is descriptive; do not infer open/closed pangenome status.",

  "",

  "Secondary metabolism:",
  "- antiSMASH detects biosynthetic gene-cluster regions, not experimentally confirmed metabolites.",
  "- Region boundaries are predictions.",
  "- Direct comparisons require the same antiSMASH version and settings for all genomes.",
  "- KnownClusterBlast similarities, when enabled, are supporting homology evidence rather than proof of metabolite production.",

  "",

  "Phylogenomic design:",
  "- A tree containing only three focal genomes has limited phylogenetic information.",
  "- For publication, add close relatives/type strains and an appropriate outgroup to genomes/ and rerun the identical workflow."
)

writeLines(
  notes,
  file.path(
    summary_dir,
    "00_analysis_notes.txt"
  )
)


cat(
  "\n====================================================\n"
)

cat(
  "Comparative genomics workflow completed.\n"
)

cat(
  "Genomes analysed:",
  nrow(
    genome_map
  ),
  "\n"
)

cat(
  "Core alignment:",
  core_alignment,
  "\n"
)

cat(
  "ML tree:",
  tree_file,
  "\n"
)

cat(
  "Panaroo results:",
  panaroo_dir,
  "\n"
)

cat(
  "antiSMASH results:",
  antismash_dir,
  "\n"
)

cat(
  "R summaries:",
  summary_dir,
  "\n"
)

cat(
  "====================================================\n"
)
