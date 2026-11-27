############################################################
# R. microplus genus-level microbiome analyses
# Data:
#   1) genus abundance: 69 samples x 1903 genera, relative abundance (rows sum ~ 1)
#   2) metadata: 69 samples
#
# Analyses:
#   A. PERMANOVA: sequential + marginal + sensitivity (Clade last)
#   B. PERMDISP
#   C. ANCOM-BC2 (publication-grade requires RAW COUNTS)
#   D. MaAsLin2
#   E. LEfSe
#
# Figures use a Nature/NPG-like palette:
#   A = #E64B35, B = #4DBBD5, C = #00A087
#
# IMPORTANT:
# - Current uploaded genus table is RELATIVE ABUNDANCE, not raw counts.
# - PERMANOVA/PERMDISP/MaAsLin2/LEfSe can use the current table.
# - ANCOM-BC2 should be run on raw genus counts for the final manuscript.
############################################################

rm(list = ls())
set.seed(20260812)

############################
# 0. File paths
############################
# Change these two paths to your local filenames.
genus_file <- "Genus.xlsx"
genus_count_file <- "Genus_Count.xlsx"
metadata_file <- "metadata.xlsx"

# OPTIONAL: raw genus count table for publication-grade ANCOM-BC2.
# Format must be identical to genus_file: SampleID + genus columns.



outdir <- "Genus_analysis_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################
# 1. Install/load packages
############################
cran_pkgs <- c(
  "readxl", "dplyr", "tidyr", "tibble", "ggplot2", "vegan",
  "ggrepel", "patchwork", "scales", "car"
)

to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "ANCOMBC", "phyloseq", "Maaslin2", "lefser",
  "SummarizedExperiment", "S4Vectors"
)
to_install_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(ggrepel)
  library(patchwork)
})

############################
# 2. Read and clean data
############################
genus_raw <- readxl::read_excel(genus_file)
meta_raw  <- readxl::read_excel(metadata_file)

stopifnot("SampleID" %in% names(genus_raw))
stopifnot("SampleID" %in% names(meta_raw))

# Rename metadata columns to R-safe, concise names.
meta <- meta_raw %>%
  rename(
    Region      = `Region`,
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
    Region = factor(Region)
  )

genus_raw$SampleID <- as.character(genus_raw$SampleID)

# Check sample matching.
if (!setequal(genus_raw$SampleID, meta$SampleID)) {
  stop("SampleID mismatch between genus abundance and metadata.")
}

# Reorder genus table exactly to metadata.
genus_raw <- genus_raw[match(meta$SampleID, genus_raw$SampleID), ]
stopifnot(all(genus_raw$SampleID == meta$SampleID))

# Convert abundance columns to numeric matrix.
genus <- genus_raw %>%
  column_to_rownames("SampleID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# Remove all-zero taxa if any.
genus <- genus[, colSums(genus, na.rm = TRUE) > 0, drop = FALSE]

# Basic checks.
row_sums <- rowSums(genus)
cat("\nRow-sum summary of genus table:\n")
print(summary(row_sums))

is_relative <- all(abs(row_sums - 1) < 0.01)
cat("\nDetected relative-abundance table:", is_relative, "\n")
cat("Samples:", nrow(genus), "\n")
cat("Genera retained:", ncol(genus), "\n\n")

write.csv(
  data.frame(SampleID = rownames(genus), row_sum = row_sums),
  file.path(outdir, "00_genus_row_sum_check.csv"),
  row.names = FALSE
)

############################
# 3. Shared style and filters
############################

# Nature/NPG-like palette.
nature_clade <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)

nature_stage <- c(
  "Adult" = "#3C5488",
  "Larva" = "#F39B7F",
  "Nymph" = "#8491B4"
)

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 1, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )
}

# For differential abundance analyses:
# keep genera present (>0) in >=10% of samples.
prevalence <- colMeans(genus > 0)
genus_da <- genus[, prevalence >= 0.10, drop = FALSE]

# Re-normalize after prevalence filtering for methods that use relative abundance.
genus_da_rel <- genus_da / rowSums(genus_da)

cat("Genera after >=10% prevalence filtering:", ncol(genus_da_rel), "\n")

write.csv(
  data.frame(Genus = names(prevalence), Prevalence = prevalence),
  file.path(outdir, "00_genus_prevalence.csv"),
  row.names = FALSE
)

############################
# 4. Metadata diagnostics
############################

# Clade x Stage table: important because Clade A contains adults only.
tab_clade_stage <- table(meta$Clade, meta$Stage)
print(tab_clade_stage)
write.csv(as.data.frame.matrix(tab_clade_stage),
          file.path(outdir, "00_Clade_by_Stage.csv"))

# Correlation among environmental variables.
env_vars <- c("Humidity", "Temperature", "Longitude", "Latitude", "Altitude")
env_cor <- cor(meta[, env_vars], use = "pairwise.complete.obs")
write.csv(env_cor, file.path(outdir, "00_environment_correlation.csv"))

# Standardize continuous covariates.
meta <- meta %>%
  mutate(
    Humidity_z    = as.numeric(scale(Humidity)),
    Temperature_z = as.numeric(scale(Temperature)),
    Longitude_z   = as.numeric(scale(Longitude)),
    Latitude_z    = as.numeric(scale(Latitude)),
    Altitude_z    = as.numeric(scale(Altitude))
  )

# VIF diagnostic using an ordinary model only as a design-matrix diagnostic.
vif_fit <- lm(
  Temperature_z ~ Humidity_z + Longitude_z + Latitude_z + Altitude_z,
  data = meta
)
vif_values <- car::vif(vif_fit)
write.csv(
  data.frame(variable = names(vif_values), VIF = as.numeric(vif_values)),
  file.path(outdir, "00_environment_VIF_diagnostic.csv"),
  row.names = FALSE
)

############################
# A. PERMANOVA + PERMDISP
############################

# Use the full genus relative-abundance table for community-level Bray-Curtis.
# Relative abundances are appropriate here; do NOT CLR-transform before Bray-Curtis.
bray <- vegan::vegdist(genus, method = "bray")

# A1. Sequential PERMANOVA (order-dependent).
# Biological-first order: Clade -> Stage -> environment.
perm_seq <- vegan::adonis2(
  bray ~ Clade + Stage + Temperature_z + Humidity_z +
    Altitude_z + Latitude_z + Longitude_z,
  data = meta,
  permutations = 9999,
  by = "terms"
)

# A2. Marginal PERMANOVA:
# Each term is tested after accounting for all other terms in the same model.
perm_margin <- vegan::adonis2(
  bray ~ Clade + Stage + Temperature_z + Humidity_z +
    Altitude_z + Latitude_z + Longitude_z,
  data = meta,
  permutations = 9999,
  by = "margin"
)

# A3. Sensitivity sequential test with Clade entered LAST.
# This directly asks how much variation Clade adds after Stage + measured environment.
perm_seq_clade_last <- vegan::adonis2(
  bray ~ Stage + Temperature_z + Humidity_z +
    Altitude_z + Latitude_z + Longitude_z + Clade,
  data = meta,
  permutations = 9999,
  by = "terms"
)

capture.output(perm_seq,
               file = file.path(outdir, "01_PERMANOVA_sequential.txt"))
capture.output(perm_margin,
               file = file.path(outdir, "01_PERMANOVA_marginal.txt"))
capture.output(perm_seq_clade_last,
               file = file.path(outdir, "01_PERMANOVA_sequential_Clade_last.txt"))

# Helper to tidy adonis2.
tidy_adonis <- function(x, test_name) {
  z <- as.data.frame(x)
  z$Term <- rownames(z)
  rownames(z) <- NULL

  pcol <- grep("^Pr", names(z), value = TRUE)[1]
  r2col <- grep("^R2$", names(z), value = TRUE)

  out <- z %>%
    filter(!Term %in% c("Residual", "Total")) %>%
    transmute(
      Test = test_name,
      Term = Term,
      Df = Df,
      R2 = if (length(r2col) == 1) .data[[r2col]] else NA_real_,
      F = F,
      P = .data[[pcol]]
    )
  out
}

perm_df <- bind_rows(
  tidy_adonis(perm_seq, "Sequential"),
  tidy_adonis(perm_margin, "Marginal")
)

write.csv(perm_df,
          file.path(outdir, "01_PERMANOVA_sequential_and_marginal.csv"),
          row.names = FALSE)

# A4. PCoA for visualization.
pcoa <- cmdscale(bray, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  SampleID = rownames(genus),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2]
) %>%
  left_join(meta, by = "SampleID")

var_exp <- 100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0])

p_pcoa <- ggplot(pcoa_df, aes(PCoA1, PCoA2, colour = Clade)) +
  geom_point(size = 3, alpha = 0.88) +
  stat_ellipse(aes(fill = Clade),
               geom = "polygon", alpha = 0.10, colour = NA,
               type = "t", level = 0.95) +
  scale_colour_manual(values = nature_clade) +
  scale_fill_manual(values = nature_clade) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", var_exp[1]),
    y = sprintf("PCoA2 (%.1f%%)", var_exp[2]),
    title = "Genus-level Bray-Curtis PCoA"
  ) +
  theme_nature()

ggsave(file.path(outdir, "01A_PCoA_Clade.pdf"), p_pcoa,
       width = 6.2, height = 5.0, device = cairo_pdf)
ggsave(file.path(outdir, "01A_PCoA_Clade.png"), p_pcoa,
       width = 6.2, height = 5.0, dpi = 600)

# A5. R2 plot: sequential vs marginal.
perm_plot_df <- perm_df %>%
  mutate(
    Term = recode(
      Term,
      "Temperature_z" = "Temperature",
      "Humidity_z" = "Humidity",
      "Altitude_z" = "Altitude",
      "Latitude_z" = "Latitude",
      "Longitude_z" = "Longitude"
    ),
    sig = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ ""
    )
  )

p_perm <- ggplot(perm_plot_df,
                 aes(x = reorder(Term, R2), y = R2, fill = Test)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sig),
    position = position_dodge(width = 0.75),
    hjust = -0.15, size = 4
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Sequential" = "#3C5488",
    "Marginal" = "#E64B35"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = expression(R^2),
    fill = NULL,
    title = "PERMANOVA effect sizes"
  ) +
  theme_nature()

ggsave(file.path(outdir, "01B_PERMANOVA_R2.pdf"), p_perm,
       width = 6.5, height = 4.8, device = cairo_pdf)
ggsave(file.path(outdir, "01B_PERMANOVA_R2.png"), p_perm,
       width = 6.5, height = 4.8, dpi = 600)

# A6. PERMDISP for Clade.
bd_clade <- vegan::betadisper(bray, group = meta$Clade, type = "median")
bd_perm <- vegan::permutest(bd_clade, permutations = 9999)
bd_tukey <- TukeyHSD(bd_clade)

capture.output(bd_perm,
               file = file.path(outdir, "01_PERMDISP_Clade_permutation.txt"))
capture.output(bd_tukey,
               file = file.path(outdir, "01_PERMDISP_Clade_TukeyHSD.txt"))

disp_df <- data.frame(
  SampleID = meta$SampleID,
  Clade = meta$Clade,
  DistanceToMedian = bd_clade$distances
)

# Overall PERMDISP P-value
bd_p <- bd_perm$tab[1, "Pr(>F)"]

p_disp <- ggplot(disp_df, aes(Clade, DistanceToMedian, fill = Clade)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.55, alpha = 0.75) +
  geom_jitter(aes(colour = Clade), width = 0.10,
              size = 2.1, alpha = 0.80) +
  scale_fill_manual(values = nature_clade) +
  scale_colour_manual(values = nature_clade) +
  labs(
    x = "Tick clade",
    y = "Distance to group median",
    title = sprintf("PERMDISP among clades (P = %.3g)", bd_p)
  ) +
  theme_nature() +
  theme(legend.position = "none")

ggsave(file.path(outdir, "01C_PERMDISP_Clade.pdf"), p_disp,
       width = 4.8, height = 4.6, device = cairo_pdf)
ggsave(file.path(outdir, "01C_PERMDISP_Clade.png"), p_disp,
       width = 4.8, height = 4.6, dpi = 600)

# Combined figure suitable for supplementary/main figure assembly.
p_permanova_combined <- p_pcoa + p_perm + p_disp +
  plot_annotation(tag_levels = "A")
ggsave(file.path(outdir, "01D_PERMANOVA_PERMDISP_combined.pdf"),
       p_permanova_combined, width = 16, height = 5.0, device = cairo_pdf)
ggsave(file.path(outdir, "01D_PERMANOVA_PERMDISP_combined.png"),
       p_permanova_combined, width = 16, height = 5.0, dpi = 600)

############################
# B. ANCOM-BC2
############################
# IMPORTANT:
# ANCOM-BC2 is a count-based bias-correction method.
# The uploaded genus table is relative abundance (rows sum ~1).
# For the final manuscript, supply RAW genus counts in genus_count_file.
#
# This block will SKIP if genus_count_file is NULL.
# Do not multiply relative abundances by an arbitrary constant and treat them
# as true counts for the final publication.
outdir <- "ANCOMBC2_genus_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# B1. Nature-style plotting
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
# B2. Read data
############################

counts0 <- readxl::read_excel(count_file)
meta0   <- readxl::read_excel(metadata_file)

stopifnot("SampleID" %in% names(counts0))
stopifnot("SampleID" %in% names(meta0))

# Rename metadata columns
meta <- meta0 %>%
  rename(
    Region      = `Region`,
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

#B3. Sample matching
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
# B4. Convert counts and QC
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
# B5. Prevalence QC
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
# B6. QC figures
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
  scale_y_continuous(labels = ) +
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
# B7. Design diagnostics
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
# B8. Construct phyloseq
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
# B9. ANCOM-BC2 model
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
# B10. Export primary results
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
# B11. Global Clade test
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
# B12. Pairwise Clade results
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
# B13. Structural zeros
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
# B14. Robust result extractor
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
# B15. Plot primary B vs A / C vs A
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
# B16. Forest plot
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
# B17. Volcano plots
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
# B18. Global-test figure
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
# B19. Pairwise-result parsing
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
# B20. Heatmap of robust genera
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
# B21. ANCOMBC Summary
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
# B22. Top significant taxa tables
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
# B23. Save package/version information
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
# B24. Messages Summary
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

############################
# C. MaAsLin2
############################
suppressPackageStartupMessages(library(Maaslin2))

# Current table is already relative abundance/TSS-normalized.
# Therefore use normalization = "NONE".
# LOG transformation is appropriate for positive skew; MaAsLin2 handles zeros internally.
maaslin_meta <- meta %>%
  select(
    SampleID, Clade, Stage,
    Temperature_z, Humidity_z, Altitude_z, Latitude_z, Longitude_z
  ) %>%
  column_to_rownames("SampleID")

# Input samples x features.
maaslin_data <- genus_da_rel
maaslin_data <- maaslin_data[rownames(maaslin_meta), , drop = FALSE]

maaslin_dir <- file.path(outdir, "03_MaAsLin2")
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)

fit_maaslin <- Maaslin2::Maaslin2(
  input_data = maaslin_data,
  input_metadata = maaslin_meta,
  output = maaslin_dir,
  fixed_effects = c(
    "Clade", "Stage",
    "Temperature_z", "Humidity_z", "Altitude_z", "Latitude_z", "Longitude_z"
  ),
  random_effects = NULL,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  standardize = FALSE,
  min_abundance = 0,
  min_prevalence = 0,
  max_significance = 0.05,
  reference = c("Clade,A", "Stage,Adult"),
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

maaslin_res <- fit_maaslin$results
write.csv(maaslin_res,
          file.path(outdir, "03_MaAsLin2_all_results.csv"),
          row.names = FALSE)

# Focus figure on adjusted Clade associations.
# In MaAsLin2 results:
# metadata = variable, value = factor level, coef = effect vs reference.
maaslin_clade <- maaslin_res %>%
  filter(metadata == "Clade", value %in% c("B", "C")) %>%
  mutate(
    qval_plot = pmax(qval, 1e-300),
    Significant = qval < 0.05,
    Contrast = paste0("Clade ", value, " vs A")
  )

write.csv(maaslin_clade,
          file.path(outdir, "03_MaAsLin2_Clade_results.csv"),
          row.names = FALSE)

maaslin_top <- maaslin_clade %>%
  filter(qval < 0.05) %>%
  group_by(Contrast) %>%
  slice_max(order_by = abs(coef), n = 15, with_ties = FALSE) %>%
  ungroup()

if (nrow(maaslin_top) > 0) {

  maaslin_top <- maaslin_top %>%
    mutate(feature = reorder(feature, coef))

  p_maaslin <- ggplot(
    maaslin_top,
    aes(x = coef, y = feature, colour = value)
  ) +
    geom_vline(xintercept = 0, linetype = 2,
               linewidth = 0.45, colour = "grey45") +
    geom_errorbar(
      aes(xmin = coef - 1.96 * stderr,
          xmax = coef + 1.96 * stderr),
      height = 0.18, linewidth = 0.5
    ) +
    geom_point(size = 2.8) +
    facet_wrap(~Contrast, scales = "free_y") +
    scale_colour_manual(values = c("B" = "#4DBBD5", "C" = "#00A087")) +
    labs(
      x = "MaAsLin2 coefficient (adjusted association)",
      y = NULL,
      colour = "Clade",
      title = "Clade-associated genera after multivariable adjustment"
    ) +
    theme_nature()

  ggsave(file.path(outdir, "03_MaAsLin2_Clade_top_genera.pdf"),
         p_maaslin, width = 9, height = 6.5, device = cairo_pdf)
  ggsave(file.path(outdir, "03_MaAsLin2_Clade_top_genera.png"),
         p_maaslin, width = 9, height = 6.5, dpi = 600)

} else {

  # Generate a diagnostic volcano-style plot even if no q<0.05 hits.
  p_maaslin <- ggplot(
    maaslin_clade,
    aes(x = coef, y = -log10(qval_plot), colour = value)
  ) +
    geom_hline(yintercept = -log10(0.05),
               linetype = 2, linewidth = 0.45) +
    geom_vline(xintercept = 0,
               linetype = 2, linewidth = 0.45) +
    geom_point(alpha = 0.75, size = 2) +
    scale_colour_manual(values = c("B" = "#4DBBD5", "C" = "#00A087")) +
    labs(
      x = "MaAsLin2 coefficient",
      y = expression(-log[10](q)),
      colour = "Clade",
      title = "MaAsLin2 clade associations"
    ) +
    theme_nature()

  ggsave(file.path(outdir, "03_MaAsLin2_Clade_volcano.pdf"),
         p_maaslin, width = 6.5, height = 5.0, device = cairo_pdf)
  ggsave(file.path(outdir, "03_MaAsLin2_Clade_volcano.png"),
         p_maaslin, width = 6.5, height = 5.0, dpi = 600)
}

############################
# D. LEfSe
############################
############################################################
# LEfSe pairwise analysis for Clade A/B/C
# Current lefser requires a dichotomous two-level class
############################################################

suppressPackageStartupMessages({
  library(lefser)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

# ----------------------------------------------------------
# 1. Nature/NPG-like palette
# ----------------------------------------------------------

nature_clade <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text =
        element_text(
          colour = "black"
        ),
      axis.title =
        element_text(
          colour = "black"
        ),
      plot.title =
        element_text(
          face = "bold"
        ),
      legend.title =
        element_text(
          face = "bold"
        ),
      panel.border =
        element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.45
        )
      
    )
}


# ----------------------------------------------------------
# 2. Output directory
# ----------------------------------------------------------
lefse_dir <- file.path(
  getwd(),
  "Genus_analysis_results",
  "04_LEfSe_pairwise"
)

dir.create(
  lefse_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ----------------------------------------------------------
# 3. Prepare genus abundance
# ----------------------------------------------------------
# genus_da_rel should already be:
# samples × genera
# and each sample should sum approximately to 1.
summary(
  rowSums(
    genus_da_rel
  )
)

# Make sure rows correspond exactly to metadata
genus_da_rel <-
  genus_da_rel[
    meta$SampleID,
    ,
    drop = FALSE
  ]

# 4. Define reusable LEfSe function
run_pairwise_lefse <- function(
    group1,
    group2,
    abundance_data,
    metadata,
    lda_cutoff = 2,
    output_dir = lefse_dir
) {
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Running LEfSe:",
    group1,
    "vs",
    group2,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  ##########################################################
  # 1. Select samples
  ##########################################################
  
  keep <-
    metadata$Clade %in%
    c(
      group1,
      group2
    )
  
  meta_sub <-
    metadata[
      keep,
      ,
      drop = FALSE
    ]
  
  abundance_sub <-
    abundance_data[
      meta_sub$SampleID,
      ,
      drop = FALSE
    ]
  
  meta_sub$Clade <-
    droplevels(
      factor(
        meta_sub$Clade,
        levels = c(
          group1,
          group2
        )
      )
    )
  
  cat(
    "\nSample number:\n"
  )
  
  print(
    table(
      meta_sub$Clade
    )
  )
  
  
  ##########################################################
  # 2. Feature filtering
  ##########################################################
  
  prevalence <-
    colMeans(
      abundance_sub > 0
    )
  
  mean_abundance <-
    colMeans(
      abundance_sub
    )
  
  keep_taxa <-
    prevalence >= 0.20 &
    mean_abundance >= 0.0001
  
  abundance_sub <-
    abundance_sub[
      ,
      keep_taxa,
      drop = FALSE
    ]
  
  cat(
    "\nAfter prevalence/abundance filtering:",
    ncol(abundance_sub),
    "genera\n"
  )
  
  
  ##########################################################
  # 3. Remove zero-variance genera
  ##########################################################
  
  taxa_variance <-
    apply(
      abundance_sub,
      2,
      var
    )
  
  abundance_sub <-
    abundance_sub[
      ,
      taxa_variance > 0,
      drop = FALSE
    ]
  
  
  ##########################################################
  # 4. Remove completely duplicated genus profiles
  ##########################################################
  
  genus_profiles <-
    as.data.frame(
      t(
        as.matrix(
          abundance_sub
        )
      )
    )
  
  duplicated_profiles <-
    duplicated(
      genus_profiles
    )
  
  cat(
    "Duplicated abundance profiles removed:",
    sum(duplicated_profiles),
    "\n"
  )
  
  genus_profiles <-
    genus_profiles[
      !duplicated_profiles,
      ,
      drop = FALSE
    ]
  
  abundance_sub <-
    as.data.frame(
      t(
        genus_profiles
      )
    )
  
  
  ##########################################################
  # 5. Build SummarizedExperiment
  ##########################################################
  
  lefse_mat <-
    t(
      as.matrix(
        abundance_sub
      )
    )
  
  lefse_se <-
    SummarizedExperiment::SummarizedExperiment(
      
      assays =
        list(
          relative_abundance =
            lefse_mat
        ),
      
      colData =
        S4Vectors::DataFrame(
          
          Clade =
            meta_sub$Clade,
          
          row.names =
            meta_sub$SampleID
        )
    )
  
  
  ##########################################################
  # 6. Confirm relative abundance
  ##########################################################
  
  lefse_se <-
    lefser::relativeAb(
      lefse_se
    )
  
  
  ##########################################################
  # 7. LEfSe
  ##########################################################
  
  set.seed(
    20260812
  )
  
  lefse_res <-
    lefser::lefser(
      
      lefse_se,
      
      classCol =
        "Clade",
      
      lda.threshold =
        lda_cutoff,
      
      method =
        "BH"
    )
  
  
  ##########################################################
  # 8. Convert result
  ##########################################################
  
  lefse_df <-
    as.data.frame(
      lefse_res
    )
  
  if (
    "features" %in%
    names(lefse_df)
  ) {
    
    lefse_df$Genus <-
      lefse_df$features
    
  } else {
    
    lefse_df$Genus <-
      rownames(
        lefse_df
      )
  }
  
  lefse_df$Comparison <-
    paste0(
      group1,
      " vs ",
      group2
    )
  
  
  ##########################################################
  # 9. Determine enriched group
  ##########################################################
  
  means <-
    abundance_sub %>%
    
    rownames_to_column(
      "SampleID"
    ) %>%
    
    left_join(
      
      meta_sub %>%
        select(
          SampleID,
          Clade
        ),
      
      by =
        "SampleID"
    ) %>%
    
    pivot_longer(
      
      -c(
        SampleID,
        Clade
      ),
      
      names_to =
        "Genus",
      
      values_to =
        "Abundance"
    ) %>%
    
    group_by(
      Genus,
      Clade
    ) %>%
    
    summarise(
      
      MeanAbundance =
        mean(
          Abundance
        ),
      
      .groups =
        "drop"
    ) %>%
    
    group_by(
      Genus
    ) %>%
    
    slice_max(
      
      order_by =
        MeanAbundance,
      
      n =
        1,
      
      with_ties =
        FALSE
    ) %>%
    
    ungroup() %>%
    
    rename(
      
      EnrichedClade =
        Clade
    )
  
  
  lefse_df <-
    lefse_df %>%
    
    left_join(
      
      means %>%
        select(
          Genus,
          EnrichedClade
        ),
      
      by =
        "Genus"
    )
  
  
  ##########################################################
  # 10. Export
  ##########################################################
  
  comparison_name <-
    paste0(
      group1,
      "_vs_",
      group2
    )
  
  write.csv(
    
    lefse_df,
    
    file.path(
      
      output_dir,
      
      paste0(
        "LEfSe_",
        comparison_name,
        "_results.csv"
      )
    ),
    
    row.names =
      FALSE
  )
  
  
  ##########################################################
  # 11. LDA plot
  ##########################################################
  
  if (
    nrow(lefse_df) > 0
  ) {
    
    lefse_plot <-
      lefse_df %>%
      
      mutate(
        
        abs_score =
          abs(
            scores
          )
        
      ) %>%
      
      slice_max(
        
        order_by =
          abs_score,
        
        n =
          30,
        
        with_ties =
          FALSE
        
      ) %>%
      
      arrange(
        abs_score
      ) %>%
      
      mutate(
        
        Genus =
          factor(
            Genus,
            levels =
              unique(
                Genus
              )
          )
      )
    
    
    p_lda <-
      ggplot(
        
        lefse_plot,
        
        aes(
          
          x =
            abs_score,
          
          y =
            Genus,
          
          fill =
            EnrichedClade
        )
      ) +
      
      geom_col(
        
        width =
          0.70,
        
        colour =
          "black",
        
        linewidth =
          0.25
      ) +
      
      scale_fill_manual(
        
        values =
          nature_clade
      ) +
      
      labs(
        
        x =
          "LEfSe LDA score",
        
        y =
          NULL,
        
        fill =
          "Enriched in",
        
        title =
          paste0(
            "LEfSe: ",
            group1,
            " vs ",
            group2
          )
      ) +
      
      theme_nature()
    
    
    ggsave(
      
      file.path(
        
        output_dir,
        
        paste0(
          "LEfSe_",
          comparison_name,
          "_LDA.pdf"
        )
      ),
      
      p_lda,
      
      width =
        7,
      
      height =
        7,
      
      device =
        cairo_pdf
    )
    
    
    ggsave(
      
      file.path(
        
        output_dir,
        
        paste0(
          "LEfSe_",
          comparison_name,
          "_LDA.png"
        )
      ),
      
      p_lda,
      
      width =
        7,
      
      height =
        7,
      
      dpi =
        600
    )
  }
  
  
  return(
    lefse_df
  )
}

############################################################
# 5. Run all three pairwise comparisons
############################################################

lefse_AB <- run_pairwise_lefse(
  group1 = "A",
  group2 = "B",
  abundance_data = genus_da_rel,
  metadata = meta
)

lefse_AC <- run_pairwise_lefse(
  group1 = "A",
  group2 = "C",
  abundance_data = genus_da_rel,
  metadata = meta
)

lefse_BC <- run_pairwise_lefse(
  group1 = "B",
  group2 = "C",
  abundance_data = genus_da_rel,
  metadata = meta
)

# 6. Combine results
lefse_all <-
  bind_rows(
    lefse_AB,
    lefse_AC,
    lefse_BC
  )

write.csv(
  lefse_all,
  file.path(
    lefse_dir,
    "LEfSe_all_pairwise_results.csv"
  ),
  row.names =
    FALSE
)


# 7. Summary of recurring biomarkers
lefse_frequency <-
  lefse_all %>%
  group_by(
    Genus
  ) %>%
  summarise(
    Comparisons =
      n_distinct(
        Comparison
      ),
    Max_LDA =
      max(
        abs(
          scores
        ),
        na.rm =
          TRUE
      ),
    .groups =
      "drop"
  ) %>%
  arrange(
    desc(
      Comparisons
    ),
    desc(
      Max_LDA
    )
    
  )

write.csv(
  lefse_frequency,
  file.path(
    lefse_dir,
    "LEfSe_biomarker_stability.csv"
  ),
  row.names =
    FALSE
)

############################
# 5. Cross-method consensus table
############################
# Useful for selecting robust clade-associated genera:
# ANCOM-BC2 (once raw counts are supplied) + MaAsLin2 + LEfSe.

maaslin_sig_names <- maaslin_clade %>%
  filter(qval < 0.05) %>%
  pull(feature) %>%
  unique()

lefse_sig_names <- unique(lefse_frequency$Genus)

consensus <- data.frame(
  Genus = union(maaslin_sig_names, lefse_sig_names)
) %>%
  mutate(
    MaAsLin2_q_lt_0.05 = Genus %in% maaslin_sig_names,
    LEfSe_LDA_ge_2 = Genus %in% lefse_sig_names,
    Methods_supported = MaAsLin2_q_lt_0.05 + LEfSe_LDA_ge_2
  ) %>%
  arrange(desc(Methods_supported), Genus)

write.csv(consensus,
          file.path(outdir, "05_consensus_MaAsLin2_LEfSe.csv"),
          row.names = FALSE)

############################
# 6. Session information
############################
capture.output(sessionInfo(),
               file = file.path(outdir, "99_sessionInfo.txt"))

message("\nAll analyses finished.")
message("Results folder: ", normalizePath(outdir))
