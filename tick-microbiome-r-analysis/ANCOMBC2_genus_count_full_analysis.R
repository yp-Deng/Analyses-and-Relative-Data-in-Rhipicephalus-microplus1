############################################################
# ANCOM-BC2 differential abundance analysis
# Rhipicephalus microplus genus count data
#
# Input:
#   genus raw counts: e397a227-ad13-4a30-811f-54da5a2cf7cf.xlsx
#   metadata:         da4d90c8-27d4-4d3b-b78a-40c8d776b5d5.xlsx
#
# Dataset checked before script generation:
#   69 samples
#   1902 genus features
#   count data are non-negative integers
#   SampleID matches metadata
#
# Main biological question:
#   Which genera are differentially abundant among Clades A, B and C
#   after adjustment for developmental stage and measured environment?
#
# Outputs:
#   QC tables + plots
#   ANCOM-BC2 primary coefficient results
#   Global Clade test
#   Pairwise Clade comparisons
#   Structural-zero results
#   Significant-genus effect plots
#   Volcano plots
#   Heatmap of robust clade-associated genera
#   Summary tables for manuscript use
#
# Palette: NPG/Nature-style
############################################################

rm(list = ls())
set.seed(20260812)

############################
# 0. Input / output
############################

count_file <- "Genus_Count.xlsx"
metadata_file <- "da4d90c8-27d4-4d3b-b78a-40c8d776b5d5.xlsx"

outdir <- "ANCOMBC2_genus_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################
# 1. Packages
############################

cran_pkgs <- c(
  "readxl", "dplyr", "tidyr", "tibble", "ggplot2",
  "ggrepel", "scales", "pheatmap", "patchwork", "car"
)

need_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(need_cran) > 0) {
  install.packages(need_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "ANCOMBC", "phyloseq", "BiocParallel"
)

need_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(need_bioc) > 0) {
  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(pheatmap)
  library(patchwork)
  library(car)
  library(ANCOMBC)
  library(phyloseq)
  library(BiocParallel)
})

############################
# 2. Nature-style plotting
############################

nature_clade <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)

nature_contrast <- c(
  "B vs A" = "#4DBBD5",
  "C vs A" = "#00A087",
  "C vs B" = "#E64B35"
)

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.45
      )
    )
}

############################
# 3. Read data
############################

counts0 <- readxl::read_excel(count_file)
meta0   <- readxl::read_excel(metadata_file)

stopifnot("SampleID" %in% names(counts0))
stopifnot("SampleID" %in% names(meta0))

# Rename metadata columns
meta <- meta0 %>%
  rename(
    Region      = `Region `,
    Humidity    = `Humidity(%)`,
    Temperature = `Temperature(℃)`,
    Longitude   = `Longtitude (E)`,
    Latitude    = `Latitude (N)`,
    Altitude    = `Altitude (m)`
  ) %>%
  mutate(
    SampleID = as.character(SampleID),
    Clade = factor(Clade, levels = c("A", "B", "C")),
    Stage = factor(Stage, levels = c("Adult", "Larva", "Nymph")),
    Region = factor(Region),
    Host = factor(Host)
  )

counts0$SampleID <- as.character(counts0$SampleID)

############################
# 4. Sample matching
############################

if (!setequal(counts0$SampleID, meta$SampleID)) {

  only_count <- setdiff(counts0$SampleID, meta$SampleID)
  only_meta  <- setdiff(meta$SampleID, counts0$SampleID)

  writeLines(
    c(
      paste("Only in count table:", paste(only_count, collapse = ", ")),
      paste("Only in metadata:", paste(only_meta, collapse = ", "))
    ),
    file.path(outdir, "ERROR_sample_mismatch.txt")
  )

  stop("SampleID mismatch. See ERROR_sample_mismatch.txt")
}

# Reorder count table to metadata
counts0 <- counts0[
  match(meta$SampleID, counts0$SampleID),
]

stopifnot(all(counts0$SampleID == meta$SampleID))

############################
# 5. Convert counts and QC
############################

count_df <- counts0 %>%
  column_to_rownames("SampleID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

if (any(is.na(count_df))) {
  stop("NA detected in genus count matrix.")
}

if (any(count_df < 0)) {
  stop("Negative values detected in genus count matrix.")
}

if (any(abs(as.matrix(count_df) -
            round(as.matrix(count_df))) > 1e-8)) {
  stop(
    "Non-integer values detected. ANCOM-BC2 requires raw counts."
  )
}

# Remove taxa with zero total counts
count_df <- count_df[
  ,
  colSums(count_df) > 0,
  drop = FALSE
]

library_size <- rowSums(count_df)

qc <- data.frame(
  SampleID = rownames(count_df),
  LibrarySize = library_size,
  Clade = meta$Clade,
  Stage = meta$Stage
)

write.csv(
  qc,
  file.path(outdir, "00_library_size_QC.csv"),
  row.names = FALSE
)

cat(
  "\nSamples:", nrow(count_df),
  "\nGenera before prevalence filtering:", ncol(count_df),
  "\nLibrary size range:",
  min(library_size), "-", max(library_size),
  "\n"
)

############################
# 6. Prevalence QC
############################

prevalence <- colMeans(count_df > 0)

prev_df <- data.frame(
  Genus = names(prevalence),
  Prevalence = as.numeric(prevalence),
  TotalCount = colSums(count_df)
) %>%
  arrange(desc(Prevalence), desc(TotalCount))

write.csv(
  prev_df,
  file.path(outdir, "00_genus_prevalence.csv"),
  row.names = FALSE
)

# ANCOM-BC2 will also perform prv_cut = 0.10 internally.
cat(
  "Genera with prevalence >=10%:",
  sum(prevalence >= 0.10),
  "\n"
)

############################
# 7. QC figures
############################

p_lib <- ggplot(
  qc,
  aes(x = Clade, y = LibrarySize, fill = Clade)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(colour = Clade),
    width = 0.10,
    size = 2.0,
    alpha = 0.80
  ) +
  scale_fill_manual(values = nature_clade) +
  scale_colour_manual(values = nature_clade) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Tick clade",
    y = "Total genus counts per sample",
    title = "Sequencing library size"
  ) +
  theme_nature() +
  theme(legend.position = "none")

ggsave(
  file.path(outdir, "00A_library_size_by_clade.pdf"),
  p_lib,
  width = 4.8,
  height = 4.5,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "00A_library_size_by_clade.png"),
  p_lib,
  width = 4.8,
  height = 4.5,
  dpi = 600
)

p_prev <- ggplot(
  prev_df,
  aes(x = Prevalence)
) +
  geom_histogram(
    bins = 30,
    fill = "#3C5488",
    colour = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0.10,
    linetype = 2,
    colour = "#E64B35",
    linewidth = 0.7
  ) +
  labs(
    x = "Genus prevalence",
    y = "Number of genera",
    title = "Genus prevalence distribution"
  ) +
  theme_nature()

ggsave(
  file.path(outdir, "00B_genus_prevalence.pdf"),
  p_prev,
  width = 5.0,
  height = 4.5,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "00B_genus_prevalence.png"),
  p_prev,
  width = 5.0,
  height = 4.5,
  dpi = 600
)

############################
# 8. Design diagnostics
############################

# IMPORTANT:
# Clade A contains adults only, whereas B/C include multiple stages.
# Therefore do NOT fit Clade*Stage as the main model.

clade_stage <- table(meta$Clade, meta$Stage)

write.csv(
  as.data.frame.matrix(clade_stage),
  file.path(outdir, "00_Clade_by_Stage.csv")
)

# Environmental correlations
env_vars <- c(
  "Humidity", "Temperature",
  "Longitude", "Latitude", "Altitude"
)

env_cor <- cor(
  meta[, env_vars],
  use = "pairwise.complete.obs"
)

write.csv(
  env_cor,
  file.path(outdir, "00_environment_correlation.csv")
)

# Scale continuous environmental covariates.
# This changes units, not biological meaning.
meta <- meta %>%
  mutate(
    Humidity_z    = as.numeric(scale(Humidity)),
    Temperature_z = as.numeric(scale(Temperature)),
    Longitude_z   = as.numeric(scale(Longitude)),
    Latitude_z    = as.numeric(scale(Latitude)),
    Altitude_z    = as.numeric(scale(Altitude))
  )

# VIF-like diagnostic for environment.
# This is a design diagnostic, not an ANCOM-BC2 result.
env_vif_fit <- lm(
  Temperature_z ~ Humidity_z +
    Longitude_z +
    Latitude_z +
    Altitude_z,
  data = meta
)

env_vif <- car::vif(env_vif_fit)

write.csv(
  data.frame(
    Variable = names(env_vif),
    VIF = as.numeric(env_vif)
  ),
  file.path(outdir, "00_environment_VIF.csv"),
  row.names = FALSE
)

############################
# 9. Construct phyloseq
############################

# phyloseq OTU table = taxa x samples
otu_mat <- t(as.matrix(count_df))

OTU <- phyloseq::otu_table(
  otu_mat,
  taxa_are_rows = TRUE
)

META <- phyloseq::sample_data(
  data.frame(
    meta,
    row.names = meta$SampleID
  )
)

ps <- phyloseq::phyloseq(
  OTU,
  META
)

############################
# 10. ANCOM-BC2 model
############################

# Main model:
# Clade = variable of interest.
# Stage + measured environmental variables = adjustment covariates.
#
# group = "Clade" enables:
#   structural-zero detection
#   global test
#   pairwise multi-group comparisons
#
# A is the reference level because:
# factor(meta$Clade, levels = c("A","B","C"))
#
# NOTE:
# pseudo_sens = TRUE is intentionally retained because current
# ANCOM-BC2 guidance recommends sensitivity analysis to reduce
# false positive discoveries.

set.seed(20260812)

ancom_out <- ANCOMBC::ancombc2(
  data = ps,
  tax_level = NULL,
  fix_formula =
    "Clade + Stage + Temperature_z + Humidity_z + Altitude_z + Latitude_z + Longitude_z",
  rand_formula = NULL,
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 0,
  s0_perc = 0.05,
  group = "Clade",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE
)

saveRDS(
  ancom_out,
  file.path(outdir, "01_ANCOMBC2_full_object.rds")
)

############################
# 11. Export primary results
############################

res_primary <- as.data.frame(ancom_out$res)

# Make sure taxon/genus is explicit
if (!"taxon" %in% names(res_primary)) {
  res_primary$taxon <- rownames(res_primary)
}

write.csv(
  res_primary,
  file.path(outdir, "01_ANCOMBC2_primary_coefficients.csv"),
  row.names = FALSE
)

############################
# 12. Global Clade test
############################

if (!is.null(ancom_out$res_global)) {

  res_global <- as.data.frame(ancom_out$res_global)

  if (!"taxon" %in% names(res_global)) {
    res_global$taxon <- rownames(res_global)
  }

  write.csv(
    res_global,
    file.path(outdir, "02_ANCOMBC2_global_Clade_test.csv"),
    row.names = FALSE
  )

} else {

  res_global <- NULL
  warning("No res_global object was returned.")
}

############################
# 13. Pairwise Clade results
############################

if (!is.null(ancom_out$res_pair)) {

  res_pair <- as.data.frame(ancom_out$res_pair)

  if (!"taxon" %in% names(res_pair)) {
    res_pair$taxon <- rownames(res_pair)
  }

  write.csv(
    res_pair,
    file.path(outdir, "03_ANCOMBC2_pairwise_Clade.csv"),
    row.names = FALSE
  )

} else {

  res_pair <- NULL
  warning("No res_pair object was returned.")
}

############################
# 14. Structural zeros
############################

if (!is.null(ancom_out$zero_ind)) {

  zero_ind <- as.data.frame(ancom_out$zero_ind)

  if (!"taxon" %in% names(zero_ind)) {
    zero_ind$taxon <- rownames(zero_ind)
  }

  write.csv(
    zero_ind,
    file.path(outdir, "04_ANCOMBC2_structural_zeros.csv"),
    row.names = FALSE
  )
}

############################
# 15. Robust result extractor
############################
# ANCOMBC output column names can vary slightly among package
# versions. These helper functions detect columns automatically.

find_first_col <- function(
    df,
    exact = character(),
    patterns = character()) {

  for (x in exact) {
    if (x %in% names(df)) return(x)
  }

  for (p in patterns) {
    hit <- grep(
      p,
      names(df),
      value = TRUE,
      ignore.case = TRUE
    )
    if (length(hit) > 0) return(hit[1])
  }

  return(NA_character_)
}

############################
# 16. Plot primary B vs A / C vs A
############################

primary_contrasts <- c(
  "CladeB" = "B vs A",
  "CladeC" = "C vs A"
)

primary_list <- list()

for (coef_name in names(primary_contrasts)) {

  label <- primary_contrasts[[coef_name]]

  lfc_col <- find_first_col(
    res_primary,
    exact = paste0("lfc_", coef_name)
  )

  se_col <- find_first_col(
    res_primary,
    exact = paste0("se_", coef_name)
  )

  q_col <- find_first_col(
    res_primary,
    exact = paste0("q_", coef_name)
  )

  p_col <- find_first_col(
    res_primary,
    exact = paste0("p_", coef_name)
  )

  diff_col <- find_first_col(
    res_primary,
    exact = paste0("diff_", coef_name)
  )

  if (is.na(lfc_col)) next

  tmp <- data.frame(
    Genus = res_primary$taxon,
    Contrast = label,
    LFC = res_primary[[lfc_col]],
    SE = if (!is.na(se_col)) res_primary[[se_col]] else NA,
    P = if (!is.na(p_col)) res_primary[[p_col]] else NA,
    Q = if (!is.na(q_col)) res_primary[[q_col]] else NA,
    stringsAsFactors = FALSE
  )

  if (!is.na(diff_col)) {
    tmp$Significant <- as.logical(
      res_primary[[diff_col]]
    )
  } else {
    tmp$Significant <- !is.na(tmp$Q) & tmp$Q < 0.05
  }

  primary_list[[label]] <- tmp
}

primary_long <- bind_rows(primary_list)

write.csv(
  primary_long,
  file.path(outdir, "05_primary_Clade_effects_tidy.csv"),
  row.names = FALSE
)

############################
# 17. Forest plot
############################

if (nrow(primary_long) > 0) {

  forest_dat <- primary_long %>%
    filter(
      Significant %in% TRUE,
      is.finite(LFC)
    ) %>%
    group_by(Contrast) %>%
    slice_max(
      order_by = abs(LFC),
      n = 15,
      with_ties = FALSE
    ) %>%
    ungroup()

  if (nrow(forest_dat) > 0) {

    forest_dat <- forest_dat %>%
      mutate(
        GenusPlot = reorder(
          paste(Contrast, Genus, sep = "__"),
          LFC
        )
      )

    p_forest <- ggplot(
      forest_dat,
      aes(
        x = LFC,
        y = GenusPlot,
        colour = Contrast
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = 2,
        colour = "grey45",
        linewidth = 0.45
      ) +
      geom_errorbarh(
        aes(
          xmin = LFC - 1.96 * SE,
          xmax = LFC + 1.96 * SE
        ),
        height = 0.15,
        linewidth = 0.50,
        na.rm = TRUE
      ) +
      geom_point(size = 2.8) +
      facet_wrap(
        ~Contrast,
        scales = "free_y"
      ) +
      scale_colour_manual(
        values = nature_contrast
      ) +
      scale_y_discrete(
        labels = function(x) sub("^.*__", "", x)
      ) +
      labs(
        x = "ANCOM-BC2 log fold change",
        y = NULL,
        colour = NULL,
        title =
          "Clade-associated genera after multivariable adjustment"
      ) +
      theme_nature() +
      theme(
        legend.position = "none"
      )

    ggsave(
      file.path(
        outdir,
        "05A_ANCOMBC2_effect_size_forest.pdf"
      ),
      p_forest,
      width = 9,
      height = 7,
      device = cairo_pdf
    )

    ggsave(
      file.path(
        outdir,
        "05A_ANCOMBC2_effect_size_forest.png"
      ),
      p_forest,
      width = 9,
      height = 7,
      dpi = 600
    )
  }
}

############################
# 18. Volcano plots
############################

if (nrow(primary_long) > 0) {

  volcano_dat <- primary_long %>%
    mutate(
      Qplot = ifelse(
        is.na(Q),
        1,
        pmax(Q, 1e-300)
      ),
      Significance = ifelse(
        Significant %in% TRUE,
        "FDR < 0.05",
        "Not significant"
      )
    )

  p_volcano <- ggplot(
    volcano_dat,
    aes(
      x = LFC,
      y = -log10(Qplot)
    )
  ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      colour = "grey45"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      colour = "grey45"
    ) +
    geom_point(
      aes(
        colour = Significance
      ),
      size = 2,
      alpha = 0.75
    ) +
    facet_wrap(
      ~Contrast,
      scales = "free"
    ) +
    scale_colour_manual(
      values = c(
        "FDR < 0.05" = "#E64B35",
        "Not significant" = "grey75"
      )
    ) +
    labs(
      x = "ANCOM-BC2 log fold change",
      y = expression(-log[10](q)),
      colour = NULL,
      title = "ANCOM-BC2 differential abundance"
    ) +
    theme_nature()

  # Label strongest significant genera
  label_dat <- volcano_dat %>%
    filter(Significant %in% TRUE) %>%
    group_by(Contrast) %>%
    slice_max(
      order_by = -log10(Qplot),
      n = 10,
      with_ties = FALSE
    ) %>%
    ungroup()

  if (nrow(label_dat) > 0) {
    p_volcano <- p_volcano +
      ggrepel::geom_text_repel(
        data = label_dat,
        aes(label = Genus),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.35,
        min.segment.length = 0
      )
  }

  ggsave(
    file.path(
      outdir,
      "05B_ANCOMBC2_volcano.pdf"
    ),
    p_volcano,
    width = 9,
    height = 5,
    device = cairo_pdf
  )

  ggsave(
    file.path(
      outdir,
      "05B_ANCOMBC2_volcano.png"
    ),
    p_volcano,
    width = 9,
    height = 5,
    dpi = 600
  )
}

############################
# 19. Global-test figure
############################

if (!is.null(res_global)) {

  # Detect relevant columns
  qg <- find_first_col(
    res_global,
    exact = c("q", "q_val", "q_global"),
    patterns = c("^q_", "q")
  )

  pg <- find_first_col(
    res_global,
    exact = c("p", "p_val", "p_global"),
    patterns = c("^p_", "p")
  )

  dg <- find_first_col(
    res_global,
    exact = c("diff", "diff_global"),
    patterns = c("^diff")
  )

  wg <- find_first_col(
    res_global,
    exact = c("W", "W_global"),
    patterns = c("^W")
  )

  global_plot <- data.frame(
    Genus = res_global$taxon,
    Q = if (!is.na(qg)) res_global[[qg]] else NA,
    P = if (!is.na(pg)) res_global[[pg]] else NA,
    W = if (!is.na(wg)) res_global[[wg]] else NA
  )

  if (!is.na(dg)) {
    global_plot$Significant <-
      as.logical(res_global[[dg]])
  } else {
    global_plot$Significant <-
      !is.na(global_plot$Q) &
      global_plot$Q < 0.05
  }

  write.csv(
    global_plot,
    file.path(
      outdir,
      "06_global_Clade_test_tidy.csv"
    ),
    row.names = FALSE
  )

  global_sig <- global_plot %>%
    filter(Significant %in% TRUE) %>%
    arrange(Q)

  if (nrow(global_sig) > 0) {

    nshow <- min(30, nrow(global_sig))

    global_top <- global_sig %>%
      slice_head(n = nshow) %>%
      mutate(
        Genus = factor(
          Genus,
          levels = rev(Genus)
        ),
        neglogQ = -log10(
          pmax(Q, 1e-300)
        )
      )

    p_global <- ggplot(
      global_top,
      aes(
        x = neglogQ,
        y = Genus
      )
    ) +
      geom_col(
        fill = "#3C5488",
        width = 0.70
      ) +
      geom_vline(
        xintercept = -log10(0.05),
        linetype = 2,
        colour = "#E64B35"
      ) +
      labs(
        x = expression(-log[10](q)),
        y = NULL,
        title =
          "Global ANCOM-BC2 test for Clade"
      ) +
      theme_nature()

    ggsave(
      file.path(
        outdir,
        "06_ANCOMBC2_global_Clade.pdf"
      ),
      p_global,
      width = 6.5,
      height = 7,
      device = cairo_pdf
    )

    ggsave(
      file.path(
        outdir,
        "06_ANCOMBC2_global_Clade.png"
      ),
      p_global,
      width = 6.5,
      height = 7,
      dpi = 600
    )
  }
}

############################
# 20. Pairwise-result parsing
############################

# This section exports pairwise output as returned by the installed
# ANCOMBC version. It also attempts to reshape standard ANCOM-BC2
# pairwise columns for plotting.

if (!is.null(res_pair)) {

  pair_col_names <- names(res_pair)

  writeLines(
    pair_col_names,
    file.path(
      outdir,
      "07_pairwise_column_names.txt"
    )
  )

  # Save all significant indicators, if present
  diff_cols <- grep(
    "^diff_",
    pair_col_names,
    value = TRUE
  )

  if (length(diff_cols) > 0) {

    pair_sig_summary <- data.frame(
      ContrastColumn = diff_cols,
      N_significant = sapply(
        diff_cols,
        function(x) {
          sum(
            res_pair[[x]] %in% TRUE,
            na.rm = TRUE
          )
        }
      )
    )

    write.csv(
      pair_sig_summary,
      file.path(
        outdir,
        "07_pairwise_significant_counts.csv"
      ),
      row.names = FALSE
    )
  }
}

############################
# 21. Heatmap of robust genera
############################

# Define robust candidate genera as:
# significant in at least one primary Clade contrast.
sig_genera <- primary_long %>%
  filter(Significant %in% TRUE) %>%
  pull(Genus) %>%
  unique()

if (length(sig_genera) > 0) {

  # Limit to most informative 40 genera for readability.
  ranking <- primary_long %>%
    filter(
      Significant %in% TRUE,
      Genus %in% sig_genera
    ) %>%
    group_by(Genus) %>%
    summarise(
      max_abs_LFC = max(
        abs(LFC),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(max_abs_LFC))

  heat_genera <- head(
    ranking$Genus,
    40
  )

  heat_counts <- count_df[
    ,
    intersect(
      heat_genera,
      colnames(count_df)
    ),
    drop = FALSE
  ]

  # Relative abundance only for visualization.
  rel <- heat_counts /
    rowSums(count_df)

  # log10 transformation + row z-score
  heat_mat <- t(
    log10(rel + 1e-8)
  )

  heat_z <- t(
    scale(
      t(heat_mat)
    )
  )

  annotation_col <- data.frame(
    Clade = meta$Clade,
    Stage = meta$Stage
  )

  rownames(annotation_col) <-
    meta$SampleID

  # Order samples by Clade then Stage
  ord <- order(
    meta$Clade,
    meta$Stage
  )

  heat_z <- heat_z[
    ,
    ord,
    drop = FALSE
  ]

  annotation_col <- annotation_col[
    ord,
    ,
    drop = FALSE
  ]

  ann_colors <- list(
    Clade = nature_clade,
    Stage = c(
      "Adult" = "#3C5488",
      "Larva" = "#F39B7F",
      "Nymph" = "#8491B4"
    )
  )

  pdf(
    file.path(
      outdir,
      "08_ANCOMBC2_significant_genera_heatmap.pdf"
    ),
    width = 11,
    height = 9
  )

  pheatmap::pheatmap(
    heat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    fontsize_row = 8,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    border_color = NA,
    color = colorRampPalette(
      c(
        "#3C5488",
        "white",
        "#E64B35"
      )
    )(100),
    main =
      "ANCOM-BC2 significant genera"
  )

  dev.off()
}

############################
# 22. Manuscript summary
############################

summary_tab <- primary_long %>%
  group_by(Contrast) %>%
  summarise(
    Tested = sum(!is.na(LFC)),
    Significant_FDR_0.05 =
      sum(Significant %in% TRUE),
    Positive_LFC =
      sum(
        Significant %in% TRUE &
        LFC > 0,
        na.rm = TRUE
      ),
    Negative_LFC =
      sum(
        Significant %in% TRUE &
        LFC < 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write.csv(
  summary_tab,
  file.path(
    outdir,
    "09_ANCOMBC2_manuscript_summary.csv"
  ),
  row.names = FALSE
)

############################
# 23. Top significant taxa tables
############################

if (nrow(primary_long) > 0) {

  top_sig <- primary_long %>%
    filter(Significant %in% TRUE) %>%
    arrange(
      Contrast,
      Q,
      desc(abs(LFC))
    )

  write.csv(
    top_sig,
    file.path(
      outdir,
      "09_ANCOMBC2_all_significant_Clade_genera.csv"
    ),
    row.names = FALSE
  )

  top20_each <- top_sig %>%
    group_by(Contrast) %>%
    slice_head(n = 20) %>%
    ungroup()

  write.csv(
    top20_each,
    file.path(
      outdir,
      "09_ANCOMBC2_top20_each_contrast.csv"
    ),
    row.names = FALSE
  )
}

############################
# 24. Save package/version information
############################

capture.output(
  packageVersion("ANCOMBC"),
  file = file.path(
    outdir,
    "99_ANCOMBC_version.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)

############################
# 25. Final messages
############################

cat("\n============================================\n")
cat("ANCOM-BC2 analysis completed.\n")
cat("Output directory:", outdir, "\n")
cat("Reference Clade: A\n")
cat("Prevalence filter: >=10%\n")
cat("FDR method: Benjamini-Hochberg\n")
cat("Sensitivity analysis: ON\n")
cat("Global Clade test: ON\n")
cat("Pairwise Clade tests: ON\n")
cat("Structural-zero detection: ON\n")
cat("============================================\n")
