############################################################
# Family-level microbiome spatial/environmental analysis
# Improved publication-ready workflow
#
# INPUT
#   Family.xlsx   : SampleID + family relative abundance columns
#   metadata.xlsx : SampleID + Clade + Stage + Temperature +
#                   Humidity + Altitude + Longitude + Latitude
#
# ANALYSES
#   1. Data QC and alignment
#   2. Bray-Curtis dissimilarity (relative abundance)
#   3. Sequential + marginal PERMANOVA
#   4. PERMDISP
#   5. dbRDA + term-wise permutation tests + VIF diagnostics
#   6. Variation partitioning:
#        Host = Clade + Stage
#        Environment = Temperature + Humidity + Altitude
#        Space = Longitude + Latitude
#      Response = Hellinger-transformed family abundance
#   7. Overall and clade-specific Mantel distance-decay
#   8. Partial Mantel controlling for clade membership
#   9. Publication-style figures and machine-readable outputs
#
# IMPORTANT STATISTICAL NOTES
# - Bray-Curtis is computed from relative-abundance data directly.
# - Hellinger transformation is used for variation partitioning/RDA.
# - Pairwise distance-decay regression lines are visualization only;
#   inferential significance is taken from Mantel permutation tests.
# - Partial Mantel uses a binary same/different-clade distance matrix;
#   A/B/C are NOT treated as ordinal numeric values.
############################################################

rm(list = ls())
set.seed(20260813)

############################
# 0. File paths
############################

family_file <- "Family.xlsx"
metadata_file <- "metadata.xlsx"

outdir <- "Family_level_ecology_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################
# 1. Packages
############################

cran_pkgs <- c(
  "readxl", "dplyr", "tidyr", "tibble", "ggplot2",
  "ggrepel", "vegan", "geosphere", "patchwork", "scales"
)

need <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(need) > 0) {
  install.packages(need)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(geosphere)
  library(patchwork)
  library(scales)
})

############################
# 2. Nature/NPG-like style
############################

clade_colors <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)

group_colors <- c(
  "Overall" = "black",
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)

vpa_colors <- c(
  "Host (pure)" = "#3C5488",
  "Environment (pure)" = "#E64B35",
  "Spatial (pure)" = "#00A087",
  "Shared" = "#8491B4",
  "Residual" = "grey80"
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

      plot.subtitle = element_text(
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.45
      ),

      strip.background = element_blank(),

      strip.text = element_text(
        face = "bold",
        colour = "black"
      )
    )
}


save_fig <- function(
    plot,
    filename,
    width = 7,
    height = 5
) {

  ggsave(
    file.path(
      outdir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    file.path(
      outdir,
      paste0(filename, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = 600
  )
}


############################
# 3. Read data
############################

family_raw <- readxl::read_excel(
  family_file
)

meta_raw <- readxl::read_excel(
  metadata_file
)

if (!"SampleID" %in% names(family_raw)) {
  stop("Family.xlsx must contain a SampleID column.")
}

if (!"SampleID" %in% names(meta_raw)) {
  stop("metadata.xlsx must contain a SampleID column.")
}


############################
# 4. Standardize metadata column names
############################
# Compatible with the metadata structure used in the previous analyses.
# If columns are already simplified, they are retained.

meta <- meta_raw

rename_if_present <- function(
    df,
    old,
    new
) {

  if (
    old %in% names(df) &&
    !new %in% names(df)
  ) {
    names(df)[names(df) == old] <- new
  }

  df
}

meta <- rename_if_present(
  meta,
  "Region",
  "Region"
)

meta <- rename_if_present(
  meta,
  "Humidity(%)",
  "Humidity"
)

meta <- rename_if_present(
  meta,
  "Temperature(℃)",
  "Temperature"
)

meta <- rename_if_present(
  meta,
  "Longtitude (E)",
  "Longitude"
)

meta <- rename_if_present(
  meta,
  "Longitude (E)",
  "Longitude"
)

meta <- rename_if_present(
  meta,
  "Latitude (N)",
  "Latitude"
)

meta <- rename_if_present(
  meta,
  "Altitude (m)",
  "Altitude"
)


required_meta <- c(
  "SampleID",
  "Clade",
  "Stage",
  "Temperature",
  "Humidity",
  "Altitude",
  "Longitude",
  "Latitude"
)

missing_meta <- setdiff(
  required_meta,
  names(meta)
)

if (length(missing_meta) > 0) {

  stop(
    "Missing required metadata columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}


############################
# 5. Convert metadata types
############################

meta <- meta %>%
  mutate(
    SampleID = as.character(SampleID),

    Clade = factor(
      Clade,
      levels = c(
        "A",
        "B",
        "C"
      )
    ),

    Stage = factor(Stage),

    across(
      c(
        Temperature,
        Humidity,
        Altitude,
        Longitude,
        Latitude
      ),
      as.numeric
    )
  )


############################
# 6. Prepare family abundance matrix
############################

family_raw$SampleID <- as.character(
  family_raw$SampleID
)

# Check duplicates
if (anyDuplicated(family_raw$SampleID) > 0) {
  stop("Duplicated SampleID found in Family.xlsx.")
}

if (anyDuplicated(meta$SampleID) > 0) {
  stop("Duplicated SampleID found in metadata.xlsx.")
}

# Use common samples only, but report any mismatch.
only_family <- setdiff(
  family_raw$SampleID,
  meta$SampleID
)

only_meta <- setdiff(
  meta$SampleID,
  family_raw$SampleID
)

if (
  length(only_family) > 0 ||
  length(only_meta) > 0
) {

  warning(
    "Sample mismatch detected. ",
    "Analysis will use common SampleIDs only."
  )

  writeLines(
    c(
      paste0(
        "Only in Family.xlsx: ",
        paste(
          only_family,
          collapse = ", "
        )
      ),

      paste0(
        "Only in metadata.xlsx: ",
        paste(
          only_meta,
          collapse = ", "
        )
      )
    ),
    file.path(
      outdir,
      "00_sample_mismatch.txt"
    )
  )
}

common_samples <- intersect(
  meta$SampleID,
  family_raw$SampleID
)

meta <- meta %>%
  filter(
    SampleID %in% common_samples
  )

family_raw <- family_raw %>%
  filter(
    SampleID %in% common_samples
  )

# Reorder family table to metadata.
family_raw <- family_raw[
  match(
    meta$SampleID,
    family_raw$SampleID
  ),
]

stopifnot(
  all(
    meta$SampleID ==
      family_raw$SampleID
  )
)

# REMOVE SampleID before numerical transformation.
family_mat <- family_raw %>%
  select(
    -SampleID
  ) %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  ) %>%
  as.data.frame()

rownames(family_mat) <- meta$SampleID

if (anyNA(family_mat)) {
  stop(
    "NA values found in family abundance matrix after numeric conversion."
  )
}

if (any(family_mat < 0)) {
  stop(
    "Negative abundance values found in Family.xlsx."
  )
}

# Remove all-zero taxa.
family_mat <- family_mat[
  ,
  colSums(
    family_mat
  ) > 0,
  drop = FALSE
]

# Check sample totals.
family_sum <- rowSums(
  family_mat
)

write.csv(
  data.frame(
    SampleID = meta$SampleID,
    FamilySum = family_sum,
    Clade = meta$Clade,
    Stage = meta$Stage
  ),
  file.path(
    outdir,
    "00_family_abundance_QC.csv"
  ),
  row.names = FALSE
)

cat(
  "\nFamily abundance row-sum summary:\n"
)

print(
  summary(
    family_sum
  )
)


############################
# 7. Convert to relative abundance if needed
############################
# If rows already sum ~1 or ~100, normalize anyway to a consistent 0-1 scale.
# This has no effect on relative composition within each sample.

if (any(family_sum <= 0)) {
  stop(
    "At least one sample has zero total family abundance."
  )
}

family_rel <- family_mat /
  rowSums(
    family_mat
  )

# Hellinger transformation used for RDA/VPA response.
family_hel <- vegan::decostand(
  family_rel,
  method = "hellinger"
)


############################
# 8. Scale continuous metadata for multivariable models
############################

meta <- meta %>%
  mutate(
    Temperature_z =
      as.numeric(
        scale(
          Temperature
        )
      ),

    Humidity_z =
      as.numeric(
        scale(
          Humidity
        )
      ),

    Altitude_z =
      as.numeric(
        scale(
          Altitude
        )
      ),

    Longitude_z =
      as.numeric(
        scale(
          Longitude
        )
      ),

    Latitude_z =
      as.numeric(
        scale(
          Latitude
        )
      )
  )


############################
# 9. Environmental correlation matrix
############################

env_cor <- cor(
  meta[
    ,
    c(
      "Temperature",
      "Humidity",
      "Altitude",
      "Longitude",
      "Latitude"
    )
  ],
  use = "pairwise.complete.obs",
  method = "spearman"
)

write.csv(
  env_cor,
  file.path(
    outdir,
    "00_environment_Spearman_correlations.csv"
  )
)


############################################################
# PART I. BRAY-CURTIS
############################################################

# Bray-Curtis directly on relative abundance.
bc_dist <- vegan::vegdist(
  family_rel,
  method = "bray"
)


############################################################
# PART II. PERMANOVA
############################################################

# Sequential test: order-dependent.
set.seed(20260813)

adonis_seq <- vegan::adonis2(
  bc_dist ~
    Clade +
    Stage +
    Temperature_z +
    Humidity_z +
    Altitude_z +
    Latitude_z +
    Longitude_z,

  data = meta,

  permutations = 9999,

  by = "terms"
)

# Marginal test: each term tested after accounting for all others.
set.seed(20260813)

adonis_margin <- vegan::adonis2(
  bc_dist ~
    Clade +
    Stage +
    Temperature_z +
    Humidity_z +
    Altitude_z +
    Latitude_z +
    Longitude_z,

  data = meta,

  permutations = 9999,

  by = "margin"
)

write.csv(
  as.data.frame(
    adonis_seq
  ),
  file.path(
    outdir,
    "01_PERMANOVA_sequential.csv"
  )
)

write.csv(
  as.data.frame(
    adonis_margin
  ),
  file.path(
    outdir,
    "01_PERMANOVA_marginal.csv"
  )
)

capture.output(
  adonis_seq,
  file = file.path(
    outdir,
    "01_PERMANOVA_sequential.txt"
  )
)

capture.output(
  adonis_margin,
  file = file.path(
    outdir,
    "01_PERMANOVA_marginal.txt"
  )
)


############################################################
# PART III. PERMDISP
############################################################

bd_clade <- vegan::betadisper(
  bc_dist,
  group = meta$Clade,
  type = "median"
)

set.seed(20260813)

bd_clade_test <- vegan::permutest(
  bd_clade,
  permutations = 9999
)

bd_clade_tukey <- TukeyHSD(
  bd_clade
)

capture.output(
  bd_clade_test,
  file = file.path(
    outdir,
    "02_PERMDISP_Clade.txt"
  )
)

capture.output(
  bd_clade_tukey,
  file = file.path(
    outdir,
    "02_PERMDISP_Clade_Tukey.csv.txt"
  )
)

disp_df <- data.frame(
  SampleID = meta$SampleID,
  Clade = meta$Clade,
  DistanceToMedian =
    bd_clade$distances
)

p_disp <- ggplot(
  disp_df,
  aes(
    x = Clade,
    y = DistanceToMedian,
    fill = Clade
  )
) +

  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.5
  ) +

  geom_jitter(
    aes(
      colour = Clade
    ),
    width = 0.08,
    size = 2,
    alpha = 0.8
  ) +

  scale_fill_manual(
    values = clade_colors
  ) +

  scale_colour_manual(
    values = clade_colors
  ) +

  labs(
    x = "Tick clade",
    y = "Distance to group median",
    title = "Family-level multivariate dispersion"
  ) +

  theme_nature() +

  theme(
    legend.position = "none"
  )

save_fig(
  p_disp,
  "02_PERMDISP_Clade",
  width = 5,
  height = 4.5
)


############################################################
# PART IV. dbRDA
############################################################

# Bray-Curtis dbRDA.
# Use scaled continuous predictors.
db_model <- vegan::capscale(
  bc_dist ~
    Clade +
    Stage +
    Temperature_z +
    Humidity_z +
    Altitude_z +
    Latitude_z +
    Longitude_z,

  data = meta,

  add = "lingoes"
)

# Overall model test
set.seed(20260813)

db_overall <- anova(
  db_model,
  permutations = 9999
)

# Sequential term tests
set.seed(20260813)

db_terms <- anova(
  db_model,
  by = "terms",
  permutations = 9999
)

# Marginal term tests
set.seed(20260813)

db_margin <- anova(
  db_model,
  by = "margin",
  permutations = 9999
)

write.csv(
  as.data.frame(
    db_overall
  ),
  file.path(
    outdir,
    "03_dbRDA_overall.csv"
  )
)

write.csv(
  as.data.frame(
    db_terms
  ),
  file.path(
    outdir,
    "03_dbRDA_terms_sequential.csv"
  )
)

write.csv(
  as.data.frame(
    db_margin
  ),
  file.path(
    outdir,
    "03_dbRDA_terms_marginal.csv"
  )
)


############################
# 10. dbRDA VIF
############################
# Large VIF values indicate unstable interpretation of individual predictors.
# Report rather than automatically deleting variables.

db_vif <- vegan::vif.cca(
  db_model
)

write.csv(
  data.frame(
    Variable = names(
      db_vif
    ),
    VIF = as.numeric(
      db_vif
    )
  ),
  file.path(
    outdir,
    "03_dbRDA_VIF.csv"
  ),
  row.names = FALSE
)


############################
# 11. dbRDA plotting scores
############################

site_scores <- as.data.frame(
  vegan::scores(
    db_model,
    display = "sites",
    choices = 1:2,
    scaling = 2
  )
)

if (ncol(site_scores) < 2) {
  stop(
    "dbRDA returned fewer than two constrained axes."
  )
}

colnames(
  site_scores
)[1:2] <- c(
  "Axis1",
  "Axis2"
)

site_scores$SampleID <- rownames(
  site_scores
)

plot_data <- site_scores %>%
  left_join(
    meta,
    by = "SampleID"
  )

bp_scores <- as.data.frame(
  vegan::scores(
    db_model,
    display = "bp",
    choices = 1:2,
    scaling = 2
  )
)

if (nrow(bp_scores) > 0) {

  colnames(
    bp_scores
  )[1:2] <- c(
    "Axis1",
    "Axis2"
  )

  bp_scores$Variable <- rownames(
    bp_scores
  )

  # Keep only continuous environmental/spatial vectors.
  allowed_bp <- c(
    "Temperature_z",
    "Humidity_z",
    "Altitude_z",
    "Latitude_z",
    "Longitude_z"
  )

  bp_scores <- bp_scores %>%
    filter(
      Variable %in% allowed_bp
    ) %>%
    mutate(
      Variable = recode(
        Variable,
        Temperature_z = "Temperature",
        Humidity_z = "Humidity",
        Altitude_z = "Altitude",
        Latitude_z = "Latitude",
        Longitude_z = "Longitude"
      )
    )
}


############################
# 12. Axis labels
############################

eig <- eigenvals(
  db_model,
  model = "constrained"
)

if (length(eig) >= 2) {

  constrained_total <- sum(
    eig[
      eig > 0
    ]
  )

  axis1_pct <-
    100 * eig[1] /
    constrained_total

  axis2_pct <-
    100 * eig[2] /
    constrained_total

} else {

  axis1_pct <- NA_real_
  axis2_pct <- NA_real_
}

xlab <- ifelse(
  is.finite(
    axis1_pct
  ),
  paste0(
    "dbRDA1 (",
    round(
      axis1_pct,
      1
    ),
    "% of constrained variation)"
  ),
  "dbRDA1"
)

ylab <- ifelse(
  is.finite(
    axis2_pct
  ),
  paste0(
    "dbRDA2 (",
    round(
      axis2_pct,
      1
    ),
    "% of constrained variation)"
  ),
  "dbRDA2"
)


############################
# 13. dbRDA plot
############################

p_db <- ggplot(
  plot_data,
  aes(
    x = Axis1,
    y = Axis2,
    colour = Clade
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey75"
  ) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey75"
  ) +

  stat_ellipse(
    aes(
      fill = Clade
    ),
    type = "t",
    geom = "polygon",
    alpha = 0.10,
    colour = NA,
    level = 0.95
  ) +

  geom_point(
    size = 3,
    alpha = 0.9
  )

if (nrow(bp_scores) > 0) {

  # Data-driven arrow scaling.
  site_range <- max(
    abs(
      c(
        plot_data$Axis1,
        plot_data$Axis2
      )
    ),
    na.rm = TRUE
  )

  bp_range <- max(
    abs(
      c(
        bp_scores$Axis1,
        bp_scores$Axis2
      )
    ),
    na.rm = TRUE
  )

  arrow_mult <- ifelse(
    bp_range > 0,
    0.65 * site_range / bp_range,
    1
  )

  p_db <- p_db +

    geom_segment(
      data = bp_scores,
      aes(
        x = 0,
        y = 0,
        xend = Axis1 *
          arrow_mult,
        yend = Axis2 *
          arrow_mult
      ),
      inherit.aes = FALSE,
      arrow = arrow(
        length = unit(
          0.18,
          "cm"
        )
      ),
      linewidth = 0.65,
      colour = "black"
    ) +

    ggrepel::geom_text_repel(
      data = bp_scores,
      aes(
        x = Axis1 *
          arrow_mult *
          1.12,
        y = Axis2 *
          arrow_mult *
          1.12,
        label = Variable
      ),
      inherit.aes = FALSE,
      size = 3.6,
      colour = "black",
      min.segment.length = 0
    )
}

p_db <- p_db +

  scale_colour_manual(
    values = clade_colors
  ) +

  scale_fill_manual(
    values = clade_colors
  ) +

  labs(
    x = xlab,
    y = ylab,
    colour = "Clade",
    fill = "Clade",
    title = "Family-level dbRDA"
  ) +

  theme_nature()

save_fig(
  p_db,
  "03_dbRDA",
  width = 7.2,
  height = 5.8
)


############################################################
# PART V. VARIATION PARTITIONING
############################################################

# Response:
# Hellinger-transformed family relative abundances.
#
# X1 Host:
#   Clade + Stage
#
# X2 Environment:
#   Temperature + Humidity + Altitude
#
# X3 Spatial:
#   Longitude + Latitude

# Use model matrices to handle factors explicitly.
host_mm <- model.matrix(
  ~ Clade + Stage,
  data = meta
)[
  ,
  -1,
  drop = FALSE
]

env_mm <- scale(
  as.matrix(
    meta[
      ,
      c(
        "Temperature",
        "Humidity",
        "Altitude"
      )
    ]
  )
)

space_mm <- scale(
  as.matrix(
    meta[
      ,
      c(
        "Longitude",
        "Latitude"
      )
    ]
  )
)

colnames(
  env_mm
) <- c(
  "Temperature",
  "Humidity",
  "Altitude"
)

colnames(
  space_mm
) <- c(
  "Longitude",
  "Latitude"
)


############################
# 14. varpart
############################

vpa <- vegan::varpart(
  family_hel,
  host_mm,
  env_mm,
  space_mm
)

capture.output(
  vpa,
  file = file.path(
    outdir,
    "04_VPA_full_output.txt"
  )
)


############################
# 15. Extract VPA fractions automatically
############################

if (
  is.null(
    vpa$part$indfract
  )
) {

  stop(
    "Could not extract individual fractions from varpart object."
  )
}

indfract <- as.data.frame(
  vpa$part$indfract
)

write.csv(
  cbind(
    Fraction = rownames(
      indfract
    ),
    indfract
  ),
  file.path(
    outdir,
    "04_VPA_individual_fractions.csv"
  ),
  row.names = FALSE
)


# For three explanatory sets, vegan labels individual fractions:
# [a] = pure X1 (Host)
# [b] = pure X2 (Environment)
# [c] = pure X3 (Space)
# [d]-[g] = shared fractions
# [h] = residual
#
# Use adjusted R-squared column.

adj_col <- grep(
  "Adj",
  names(
    indfract
  ),
  value = TRUE
)[1]

if (is.na(adj_col)) {

  stop(
    "Adjusted R-squared column not found in varpart individual fractions."
  )
}

get_fraction <- function(
    pattern
) {

  idx <- grep(
    pattern,
    rownames(
      indfract
    ),
    fixed = TRUE
  )

  if (length(idx) == 0) {
    return(
      NA_real_
    )
  }

  as.numeric(
    indfract[
      idx[1],
      adj_col
    ]
  )
}

pure_host <- get_fraction(
  "[a]"
)

pure_env <- get_fraction(
  "[b]"
)

pure_space <- get_fraction(
  "[c]"
)

shared_parts <- c(
  get_fraction("[d]"),
  get_fraction("[e]"),
  get_fraction("[f]"),
  get_fraction("[g]")
)

shared_total <- sum(
  shared_parts,
  na.rm = TRUE
)

# Sum all explained adjusted fractions.
explained_total <- sum(
  c(
    pure_host,
    pure_env,
    pure_space,
    shared_parts
  ),
  na.rm = TRUE
)

residual <- 1 -
  explained_total

vpa_df <- data.frame(
  Component = c(
    "Host (pure)",
    "Environment (pure)",
    "Spatial (pure)",
    "Shared",
    "Residual"
  ),

  Adjusted_R2 = c(
    pure_host,
    pure_env,
    pure_space,
    shared_total,
    residual
  )
) %>%
  mutate(
    Percent =
      100 *
      Adjusted_R2
  )

write.csv(
  vpa_df,
  file.path(
    outdir,
    "04_VPA_summary.csv"
  ),
  row.names = FALSE
)


############################
# 16. Test pure VPA fractions
############################
# Because varpart itself partitions variation but does not provide
# permutation P values for every pure fraction.

# Pure host | environment + space
rda_host_pure <- vegan::rda(
  family_hel ~
    Clade +
    Stage +
    Condition(
      Temperature_z,
      Humidity_z,
      Altitude_z,
      Longitude_z,
      Latitude_z
    ),
  data = meta
)

set.seed(20260813)

test_host_pure <- anova(
  rda_host_pure,
  permutations = 9999
)


# Pure environment | host + space
rda_env_pure <- vegan::rda(
  family_hel ~
    Temperature_z +
    Humidity_z +
    Altitude_z +
    Condition(
      Clade,
      Stage,
      Longitude_z,
      Latitude_z
    ),
  data = meta
)

set.seed(20260813)

test_env_pure <- anova(
  rda_env_pure,
  permutations = 9999
)


# Pure spatial | host + environment
rda_space_pure <- vegan::rda(
  family_hel ~
    Longitude_z +
    Latitude_z +
    Condition(
      Clade,
      Stage,
      Temperature_z,
      Humidity_z,
      Altitude_z
    ),
  data = meta
)

set.seed(20260813)

test_space_pure <- anova(
  rda_space_pure,
  permutations = 9999
)


capture.output(
  test_host_pure,
  file = file.path(
    outdir,
    "04_VPA_test_pure_host.txt"
  )
)

capture.output(
  test_env_pure,
  file = file.path(
    outdir,
    "04_VPA_test_pure_environment.txt"
  )
)

capture.output(
  test_space_pure,
  file = file.path(
    outdir,
    "04_VPA_test_pure_spatial.txt"
  )
)


############################
# 17. VPA bar plot
############################

vpa_plot_df <- vpa_df %>%
  mutate(
    Component = factor(
      Component,
      levels = rev(
        c(
          "Host (pure)",
          "Environment (pure)",
          "Spatial (pure)",
          "Shared",
          "Residual"
        )
      )
    )
  )

p_vpa <- ggplot(
  vpa_plot_df,
  aes(
    x = Percent,
    y = Component,
    fill = Component
  )
) +

  geom_vline(
    xintercept = 0,
    colour = "grey60",
    linewidth = 0.4
  ) +

  geom_col(
    width = 0.65
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        Percent
      )
    ),
    hjust = ifelse(
      vpa_plot_df$Percent >= 0,
      -0.1,
      1.1
    ),
    size = 3.5
  ) +

  scale_fill_manual(
    values = vpa_colors
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.12,
        0.18
      )
    )
  ) +

  labs(
    x = "Adjusted explained variation (%)",
    y = NULL,
    title = "Variation partitioning of family-level microbiome"
  ) +

  theme_nature() +

  theme(
    legend.position = "none"
  )

save_fig(
  p_vpa,
  "04_VPA",
  width = 7,
  height = 4.8
)


############################################################
# PART VI. GEOGRAPHIC DISTANCE
############################################################

coords <- as.matrix(
  meta[
    ,
    c(
      "Longitude",
      "Latitude"
    )
  ]
)

if (anyNA(coords)) {
  stop(
    "Missing Longitude/Latitude values detected."
  )
}

geo_matrix_m <- geosphere::distm(
  coords,
  fun = geosphere::distHaversine
)

rownames(
  geo_matrix_m
) <- meta$SampleID

colnames(
  geo_matrix_m
) <- meta$SampleID

geo_dist <- as.dist(
  geo_matrix_m /
    1000
)

stopifnot(
  identical(
    labels(
      bc_dist
    ),
    labels(
      geo_dist
    )
  )
)


############################################################
# PART VII. OVERALL MANTEL DISTANCE-DECAY
############################################################

set.seed(20260813)

mantel_overall <- vegan::mantel(
  bc_dist,
  geo_dist,
  method = "spearman",
  permutations = 9999
)

mantel_overall_df <- data.frame(
  Group = "Overall",

  N_samples =
    nrow(
      meta
    ),

  Mantel_r =
    as.numeric(
      mantel_overall$statistic
    ),

  P_value =
    mantel_overall$signif
)

write.csv(
  mantel_overall_df,
  file.path(
    outdir,
    "05_Mantel_overall.csv"
  ),
  row.names = FALSE
)


############################################################
# PART VIII. CLADE-SPECIFIC MANTEL
############################################################

mantel_clade_list <- list()

for (
  cl in levels(
    droplevels(
      meta$Clade
    )
  )
) {

  idx <- which(
    meta$Clade == cl
  )

  # Do not run if too few samples.
  if (length(idx) < 4) {

    mantel_clade_list[[cl]] <-
      data.frame(
        Group = cl,
        N_samples = length(idx),
        Mantel_r = NA_real_,
        P_value = NA_real_
      )

    next
  }

  bc_sub <- as.dist(
    as.matrix(
      bc_dist
    )[
      idx,
      idx
    ]
  )

  geo_sub <- as.dist(
    as.matrix(
      geo_dist
    )[
      idx,
      idx
    ]
  )

  set.seed(
    20260813 +
      match(
        cl,
        levels(
          meta$Clade
        )
      )
  )

  m <- vegan::mantel(
    bc_sub,
    geo_sub,
    method = "spearman",
    permutations = 9999
  )

  mantel_clade_list[[cl]] <-
    data.frame(
      Group = cl,
      N_samples = length(idx),
      Mantel_r =
        as.numeric(
          m$statistic
        ),
      P_value =
        m$signif
    )
}

mantel_clade_df <- bind_rows(
  mantel_clade_list
)

all_mantel_results <- bind_rows(
  mantel_overall_df,
  mantel_clade_df
)

write.csv(
  all_mantel_results,
  file.path(
    outdir,
    "05_Mantel_overall_and_clade_specific.csv"
  ),
  row.names = FALSE
)


############################################################
# PART IX. PARTIAL MANTEL CONTROLLED FOR CLADE
############################################################

# Binary categorical distance:
# 0 = same clade
# 1 = different clade
#
# This avoids the incorrect ordinal coding:
# A=1, B=2, C=3.

clade_vector <- as.character(
  meta$Clade
)

clade_matrix <- outer(
  clade_vector,
  clade_vector,
  FUN = function(
      x,
      y
  ) {
    as.numeric(
      x != y
    )
  }
)

rownames(
  clade_matrix
) <- meta$SampleID

colnames(
  clade_matrix
) <- meta$SampleID

clade_dist <- as.dist(
  clade_matrix
)

set.seed(20260813)

partial_mantel <- vegan::mantel.partial(
  bc_dist,
  geo_dist,
  clade_dist,
  method = "spearman",
  permutations = 9999
)

partial_mantel_df <- data.frame(
  Test =
    "Geographic distance controlling for binary clade membership",

  Mantel_r =
    as.numeric(
      partial_mantel$statistic
    ),

  P_value =
    partial_mantel$signif
)

write.csv(
  partial_mantel_df,
  file.path(
    outdir,
    "06_Partial_Mantel_control_Clade.csv"
  ),
  row.names = FALSE
)


############################################################
# PART X. DISTANCE-DECAY VISUALIZATION DATA
############################################################

############################################################
# 1. Overall unique pairwise observations
############################################################

overall_pairs <- data.frame(
  
  Geographic_distance =
    as.vector(
      geo_dist
    ),
  
  Microbial_distance =
    as.vector(
      bc_dist
    ),
  
  Group =
    "Overall",
  
  stringsAsFactors =
    FALSE
)


############################################################
# 2. Check overall data
############################################################

cat("\nOverall pairwise comparisons:", nrow(overall_pairs), "\n")

stopifnot(
  nrow(overall_pairs) ==
    length(
      geo_dist
    )
)

stopifnot(
  nrow(overall_pairs) ==
    length(
      bc_dist
    )
)

############################################################
# 3. Clade-specific unique pairwise observations
############################################################

clade_pair_list <- list()


for (
  cl in levels(
    droplevels(
      meta$Clade
    )
  )
) {
  
  ##########################################################
  # Samples belonging to this clade
  ##########################################################
  
  idx <- which(
    meta$Clade == cl
  )
  
  n_clade <- length(
    idx
  )
  
  
  cat(
    "\nClade",
    cl,
    ":",
    n_clade,
    "samples\n"
  )
  
  
  ##########################################################
  # Skip groups with fewer than 2 samples
  ##########################################################
  
  if (
    n_clade < 2
  ) {
    
    warning(
      "Clade ",
      cl,
      " has fewer than 2 samples and was skipped."
    )
    
    next
  }
  
  
  ##########################################################
  # Extract Bray-Curtis submatrix
  ##########################################################
  
  bc_sub_mat <- as.matrix(
    bc_dist
  )[
    idx,
    idx,
    drop = FALSE
  ]
  
  
  ##########################################################
  # Extract geographic-distance submatrix
  ##########################################################
  
  geo_sub_mat <- as.matrix(
    geo_dist
  )[
    idx,
    idx,
    drop = FALSE
  ]
  
  
  ##########################################################
  # Keep lower triangle only
  #
  # This guarantees that every sample pair appears once:
  #
  # sample1-sample2
  # sample1-sample3
  # ...
  ##########################################################
  
  keep <- lower.tri(
    bc_sub_mat,
    diag = FALSE
  )
  
  
  ##########################################################
  # Construct clade-specific pairwise data
  ##########################################################
  
  temp_df <- data.frame(
    
    Geographic_distance =
      geo_sub_mat[
        keep
      ],
    
    Microbial_distance =
      bc_sub_mat[
        keep
      ],
    
    Group =
      cl,
    
    stringsAsFactors =
      FALSE
  )
  
  
  ##########################################################
  # Expected number of pairs
  ##########################################################
  
  expected_pairs <-
    choose(
      n_clade,
      2
    )
  
  
  cat(
    "Expected pairwise comparisons:",
    expected_pairs,
    "\n"
  )
  
  cat(
    "Observed pairwise comparisons:",
    nrow(temp_df),
    "\n"
  )
  
  
  stopifnot(
    nrow(temp_df) ==
      expected_pairs
  )
  
  
  ##########################################################
  # Save in list
  ##########################################################
  
  clade_pair_list[[cl]] <- temp_df
  
}


############################################################
# 4. Combine clade-specific data first
############################################################

clade_pairs <- dplyr::bind_rows(
  clade_pair_list
)


############################################################
# 5. Combine overall + clade-specific data
############################################################

distance_decay_df <- dplyr::bind_rows(
  overall_pairs,
  clade_pairs
)


############################################################
# 6. Set group order
############################################################

distance_decay_df$Group <- factor(
  distance_decay_df$Group,
  levels = c(
    "Overall",
    "A",
    "B",
    "C"
  )
)


############################################################
# 7. Check final data
############################################################

cat("\n==========================================\n")
cat(  "Distance-decay visualization dataset\n")
cat("==========================================\n")
print(
  table(
    distance_decay_df$Group
  )
)


############################################################
# 8. Check missing/infinite values
############################################################

qc_decay <- distance_decay_df %>%
  summarise(
    N =
      n(),
    Missing_geo =
      sum(
        is.na(
          Geographic_distance
        )
      ),
    Missing_microbiome = sum(
        is.na(
          Microbial_distance
        )
      ),
    Infinite_geo =
      sum(
        !is.finite(
          Geographic_distance
        )
      ),
    Infinite_microbiome =
      sum(
        !is.finite(
          Microbial_distance
        )
      )
  )

print(
  qc_decay
)


############################################################
# 9. Remove invalid values if present
############################################################
distance_decay_df <- distance_decay_df %>%
  filter(
    is.finite(
      Geographic_distance
    ),
    is.finite(
      Microbial_distance
    )
  )


############################################################
# 10. Export pairwise plotting data
############################################################
write.csv(
  distance_decay_df,
  file.path(
    outdir,
    "07_distance_decay_pairwise_data.csv"
  ),
  row.names =
    FALSE
)

############################################################
# PART XI. AUTOMATIC MANTEL LABELS
############################################################

mantel_labels <- all_mantel_results %>%
  mutate(
    label = paste0(
      Group,
      ": Mantel r = ",
      sprintf(
        "%.3f",
        Mantel_r
      ),
      ", P = ",
      ifelse(
        is.na(
          P_value
        ),
        "NA",
        ifelse(
          P_value < 0.0001,
          "<0.0001",
          sprintf(
            "%.4f",
            P_value
          )
        )
      )
    )
  )

write.csv(
  mantel_labels,
  file.path(
    outdir,
    "07_Mantel_plot_labels.csv"
  ),
  row.names = FALSE
)


############################################################
# PART XII. FACETED DISTANCE-DECAY FIGURE
############################################################

# Regression lines are descriptive visualization only.
# Mantel permutation P values provide inference.

p_decay_facet <- ggplot(
  distance_decay_df,
  aes(
    x = Geographic_distance,
    y = Microbial_distance,
    colour = Group
  )
) +

  geom_point(
    alpha = 0.25,
    size = 1.4
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.9
  ) +

  facet_wrap(
    ~ Group,
    ncol = 2,
    scales = "free_x"
  ) +

  scale_colour_manual(
    values = group_colors,
    drop = FALSE
  ) +

  labs(
    x = "Geographic distance (km)",
    y = "Bray-Curtis dissimilarity",
    colour = NULL,
    title = "Family-level microbiome distance-decay"
  ) +

  theme_nature() +

  theme(
    legend.position = "none"
  )

# Add group-specific labels in each facet.
label_positions <- distance_decay_df %>%
  group_by(
    Group
  ) %>%
  summarise(
    x = min(
      Geographic_distance,
      na.rm = TRUE
    ),

    y = max(
      Microbial_distance,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  left_join(
    mantel_labels %>%
      select(
        Group,
        label
      ),
    by = "Group"
  )

p_decay_facet <- p_decay_facet +

  geom_text(
    data = label_positions,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.3
  )

save_fig(
  p_decay_facet,
  "07_distance_decay_faceted",
  width = 8.5,
  height = 7
)


############################################################
# PART XIII. COMBINED DISTANCE-DECAY FIGURE
############################################################

p_decay_combined <- ggplot(
  distance_decay_df,
  aes(
    x = Geographic_distance,
    y = Microbial_distance,
    colour = Group
  )
) +

  geom_point(
    alpha = 0.15,
    size = 1.2
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 1
  ) +

  scale_colour_manual(
    values = group_colors
  ) +

  labs(
    x = "Geographic distance (km)",
    y = "Bray-Curtis dissimilarity",
    colour = NULL,
    title = "Overall and clade-specific distance-decay relationships"
  ) +

  theme_nature()

save_fig(
  p_decay_combined,
  "07_distance_decay_combined",
  width = 7.5,
  height = 5.8
)


############################################################
# PART XIV. PARTIAL MANTEL SUMMARY FIGURE
############################################################

mantel_summary_plot <- bind_rows(

  all_mantel_results %>%
    transmute(
      Analysis =
        paste0(
          "Mantel: ",
          Group
        ),

      Mantel_r,

      P_value
    ),

  partial_mantel_df %>%
    transmute(
      Analysis =
        "Partial Mantel:\ngeography | clade",

      Mantel_r,

      P_value
    )
) %>%
  mutate(
    Significant =
      !is.na(
        P_value
      ) &
      P_value < 0.05,

    Analysis = factor(
      Analysis,
      levels = rev(
        Analysis
      )
    )
  )

p_mantel_summary <- ggplot(
  mantel_summary_plot,
  aes(
    x = Mantel_r,
    y = Analysis,
    fill = Significant
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey60"
  ) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label = ifelse(
        is.na(
          P_value
        ),
        "P = NA",
        ifelse(
          P_value < 0.0001,
          "P < 0.0001",
          paste0(
            "P = ",
            sprintf(
              "%.4f",
              P_value
            )
          )
        )
      )
    ),
    hjust = ifelse(
      mantel_summary_plot$Mantel_r >= 0,
      -0.08,
      1.08
    ),
    size = 3.4
  ) +

  scale_fill_manual(
    values = c(
      "TRUE" = "#E64B35",
      "FALSE" = "grey75"
    )
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.18,
        0.25
      )
    )
  ) +

  labs(
    x = "Spearman Mantel r",
    y = NULL,
    fill = "P < 0.05",
    title = "Spatial distance-decay tests"
  ) +

  theme_nature()

save_fig(
  p_mantel_summary,
  "08_Mantel_summary",
  width = 7,
  height = 4.8
)


############################################################
# PART XV. COMBINED dbRDA + VPA FIGURE
############################################################

final_db_vpa <-
  p_db /
  p_vpa +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) +
  patchwork::plot_layout(
    heights = c(
      1.35,
      1
    )
  )

ggsave(
  file.path(
    outdir,
    "09_dbRDA_VPA_combined.pdf"
  ),
  final_db_vpa,
  width = 8,
  height = 10,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "09_dbRDA_VPA_combined.png"
  ),
  final_db_vpa,
  width = 8,
  height = 10,
  dpi = 600
)


############################################################
# PART XVI. REPRODUCIBILITY
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)


############################################################
# PART XVII. FINAL CONSOLE SUMMARY
############################################################

cat("\n======================================================\n")
cat("Family-level ecological analyses completed.\n")
cat("Samples analysed:", nrow(meta), "\n")
cat("Families analysed:", ncol(family_rel), "\n")
cat("Output directory:", normalizePath(outdir), "\n")
cat("\nMantel results:\n")

print(all_mantel_results)
cat("\nPartial Mantel:\n")
print(partial_mantel_df)
cat("\nVPA summary:\n")
print(vpa_df)
cat("\nImportant: inspect 03_dbRDA_VIF.csv before interpreting ", "individual dbRDA coefficients.\n")
cat("======================================================\n")
