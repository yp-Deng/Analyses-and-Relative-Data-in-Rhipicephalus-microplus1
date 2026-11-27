############################################################
# Integrated genus evidence matrix
# Actual uploaded result-column version
#
# Sources used:
#   ANCOM-BC2 pairwise:
#     taxon, lfc_CladeB, lfc_CladeC, lfc_CladeC_CladeB,
#     q_*, diff_robust_*
#
#   MaAsLin2:
#     feature, metadata, value, coef, qval, Significant, Contrast
#
#   LEfSe:
#     Genus, Comparison, scores, EnrichedClade
#
#   Random forest stability:
#     Genus, SelectionCount, SelectionFrequency,
#     MeanImportance, MedianImportance
#
#   SPIEC-EASI node metrics:
#     Genus, Degree, Betweenness, MeanRelativeAbundance,
#     EnrichedClade, Mean
#
# Evidence classes:
#   1. Robust association
#   2. Stable discriminator
#   3. Network-central
#   4. Multi-layer evidence
#
# Classification logic:
# Robust association =
#   ANCOM-BC2 robust AND MaAsLin2 q < 0.05
#   AND same direction in >=1 matched contrast
#
# Stable discriminator =
#   RF SelectionFrequency >= 0.80
#
# Network-central =
#   Degree >= 90th percentile OR
#   Betweenness >= 90th percentile
#
# Multi-layer evidence =
#   >=2 supported evidence layers among:
#     robust association / LEfSe / stable RF / network central
#
# IMPORTANT:
# LEfSe alone is treated as supporting biomarker evidence and does not
# define a standalone final class.
############################################################

rm(list = ls())

############################
# 1. Packages
############################

pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "scales",
  "stringr", "forcats", "patchwork"
)

need <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(need) > 0) {
  install.packages(need)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(forcats)
  library(patchwork)
})

############################
# 2. Locate files
############################
ancom_file <- "03_ANCOMBC2_pairwise_Clade.csv"
maaslin_file <- "03_MaAsLin2_Clade_results.csv"
lefse_file <- "LEfSe_all_pairwise_results.csv"
rf_file <- "04E_RF_feature_stability.csv"
network_file <- "02D_SPIECEASI_node_metrics.csv"

outdir <- "Consensus_evidence_matrix"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

############################
# 3. Read files
############################

ancom <- read.csv(
  ancom_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

maas <- read.csv(
  maaslin_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

lefse <- read.csv(
  lefse_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rf <- read.csv(
  rf_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

network <- read.csv(
  network_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################
# 4. Verify actual columns
############################

required_ancom <- c(
  "taxon",
  "lfc_CladeB",
  "lfc_CladeC",
  "lfc_CladeC_CladeB",
  "q_CladeB",
  "q_CladeC",
  "q_CladeC_CladeB",
  "diff_robust_CladeB",
  "diff_robust_CladeC",
  "diff_robust_CladeC_CladeB"
)

required_maas <- c(
  "feature",
  "metadata",
  "value",
  "coef",
  "qval",
  "Significant",
  "Contrast"
)

required_lefse <- c(
  "Genus",
  "scores",
  "Comparison",
  "EnrichedClade"
)

required_rf <- c(
  "Genus",
  "SelectionCount",
  "SelectionFrequency",
  "MeanImportance",
  "MedianImportance"
)

required_network <- c(
  "Genus",
  "Degree",
  "Betweenness",
  "MeanRelativeAbundance",
  "EnrichedClade"
)


check_cols <- function(df, needed, label) {

  missing <- setdiff(needed, names(df))

  if (length(missing) > 0) {

    stop(
      label,
      " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

check_cols(ancom, required_ancom, "ANCOM-BC2")
check_cols(maas, required_maas, "MaAsLin2")
check_cols(lefse, required_lefse, "LEfSe")
check_cols(rf, required_rf, "Random forest")
check_cols(network, required_network, "Network")

############################
# 5. Nature/NPG palette
############################

clade_colours <- c(
  "A enriched" = "#E64B35",
  "B enriched" = "#4DBBD5",
  "C enriched" = "#00A087"
)

neutral_colours <- c(
  "RF stability" = "#3C5488",
  "Network degree" = "#7E6148",
  "Network betweenness" = "#B09C85"
)

all_evidence_colours <- c(
  clade_colours,
  neutral_colours
)

class_colours <- c(
  "Robust association" = "#E64B35",
  "Stable discriminator" = "#3C5488",
  "Network-central" = "#7E6148",
  "Multi-layer evidence" = "#00A087"
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

      panel.grid.major = element_line(
        colour = "grey92",
        linewidth = 0.30
      ),

      panel.grid.minor = element_blank()
    )
}

############################
# 6. ANCOM-BC2 pairwise -> long format
############################

ancom_long <- bind_rows(

  ancom %>%
    transmute(
      Genus = taxon,
      Contrast = "B vs A",
      Effect = lfc_CladeB,
      Q = q_CladeB,
      Robust = as.logical(diff_robust_CladeB),
      Direction = case_when(
        lfc_CladeB < 0 ~ "A enriched",
        lfc_CladeB > 0 ~ "B enriched",
        TRUE ~ NA_character_
      )
    ),

  ancom %>%
    transmute(
      Genus = taxon,
      Contrast = "C vs A",
      Effect = lfc_CladeC,
      Q = q_CladeC,
      Robust = as.logical(diff_robust_CladeC),
      Direction = case_when(
        lfc_CladeC < 0 ~ "A enriched",
        lfc_CladeC > 0 ~ "C enriched",
        TRUE ~ NA_character_
      )
    ),

  ancom %>%
    transmute(
      Genus = taxon,
      Contrast = "C vs B",
      Effect = lfc_CladeC_CladeB,
      Q = q_CladeC_CladeB,
      Robust = as.logical(diff_robust_CladeC_CladeB),
      Direction = case_when(
        lfc_CladeC_CladeB < 0 ~ "B enriched",
        lfc_CladeC_CladeB > 0 ~ "C enriched",
        TRUE ~ NA_character_
      )
    )

) %>%
  mutate(
    ANCOM_strength = pmin(
      -log10(pmax(Q, 1e-300)),
      20
    )
  )

############################
# 7. MaAsLin2 -> matched contrasts
############################

# Actual file uses value B/C and Contrast "Clade B vs A", etc.
maas_clade <- maas %>%
  filter(
    metadata == "Clade",
    value %in% c("B", "C")
  ) %>%
  transmute(
    Genus = feature,
    Contrast = case_when(
      value == "B" ~ "B vs A",
      value == "C" ~ "C vs A",
      TRUE ~ NA_character_
    ),
    Effect = coef,
    Q = qval,
    Significant = as.logical(Significant) | qval < 0.05,
    Direction = case_when(
      value == "B" & coef < 0 ~ "A enriched",
      value == "B" & coef > 0 ~ "B enriched",
      value == "C" & coef < 0 ~ "A enriched",
      value == "C" & coef > 0 ~ "C enriched",
      TRUE ~ NA_character_
    ),
    Maaslin_strength = pmin(
      -log10(pmax(qval, 1e-300)),
      20
    )
  )

############################
# 8. ANCOM + MaAsLin2 direction concordance
############################

da_match <- ancom_long %>%
  filter(
    Contrast %in% c("B vs A", "C vs A")
  ) %>%
  select(
    Genus,
    Contrast,
    ANCOM_effect = Effect,
    ANCOM_Q = Q,
    ANCOM_robust = Robust,
    ANCOM_direction = Direction
  ) %>%
  full_join(
    maas_clade %>%
      select(
        Genus,
        Contrast,
        Maaslin_effect = Effect,
        Maaslin_Q = Q,
        Maaslin_sig = Significant,
        Maaslin_direction = Direction
      ),
    by = c(
      "Genus",
      "Contrast"
    )
  ) %>%
  mutate(
    DirectionConcordant =
      !is.na(ANCOM_direction) &
      !is.na(Maaslin_direction) &
      ANCOM_direction == Maaslin_direction,

    RobustAssociationContrast =
      replace_na(ANCOM_robust, FALSE) &
      replace_na(Maaslin_sig, FALSE) &
      replace_na(DirectionConcordant, FALSE),

    DirectionConflict =
      replace_na(ANCOM_robust, FALSE) &
      replace_na(Maaslin_sig, FALSE) &
      !replace_na(DirectionConcordant, FALSE)
  )

write.csv(
  da_match,
  file.path(
    outdir,
    "01_ANCOM_MaAsLin2_direction_check.csv"
  ),
  row.names = FALSE
)

da_summary <- da_match %>%
  group_by(Genus) %>%
  summarise(
    RobustAssociation =
      any(RobustAssociationContrast, na.rm = TRUE),

    DA_DirectionConflict =
      any(DirectionConflict, na.rm = TRUE),

    ConcordantContrasts =
      paste(
        Contrast[RobustAssociationContrast %in% TRUE],
        collapse = "; "
      ),

    ConflictContrasts =
      paste(
        Contrast[DirectionConflict %in% TRUE],
        collapse = "; "
      ),

    .groups = "drop"
  )

############################
# 9. LEfSe support
############################

lefse_summary <- lefse %>%
  mutate(
    LEfSeDirection = paste0(
      EnrichedClade,
      " enriched"
    )
  ) %>%
  group_by(Genus) %>%
  summarise(
    LEfSeSupport = n() > 0,
    LEfSeNComparisons = n_distinct(Comparison),
    MaxLEfSeLDA = max(abs(scores), na.rm = TRUE),
    LEfSeComparisons = paste(
      unique(
        paste0(
          Comparison,
          " (",
          EnrichedClade,
          ")"
        )
      ),
      collapse = "; "
    ),
    .groups = "drop"
  )

############################
# 10. RF stability
############################

rf_threshold <- 0.80

rf_summary <- rf %>%
  transmute(
    Genus,
    SelectionCount,
    SelectionFrequency,
    MeanImportance,
    MedianImportance,
    RFStable =
      SelectionFrequency >= rf_threshold
  )

############################
# 11. Network centrality
############################

# Use upper 10% of the observed network as a stringent "central" threshold.
degree_cut <- as.numeric(
  quantile(
    network$Degree,
    probs = 0.90,
    na.rm = TRUE
  )
)

betweenness_cut <- as.numeric(
  quantile(
    network$Betweenness,
    probs = 0.90,
    na.rm = TRUE
  )
)

network_summary <- network %>%
  transmute(
    Genus,
    Degree,
    Betweenness,
    MeanRelativeAbundance,
    NetworkEnrichedClade = EnrichedClade,
    DegreeTop10 =
      Degree >= degree_cut,
    BetweennessTop10 =
      Betweenness >= betweenness_cut,
    NetworkCentral =
      Degree >= degree_cut |
      Betweenness >= betweenness_cut
  )

writeLines(
  c(
    paste0(
      "RF stable threshold: SelectionFrequency >= ",
      rf_threshold
    ),
    paste0(
      "Network Degree 90th percentile = ",
      degree_cut
    ),
    paste0(
      "Network Betweenness 90th percentile = ",
      signif(betweenness_cut, 5)
    )
  ),
  file.path(
    outdir,
    "00_thresholds.txt"
  )
)

############################
# 12. Build master genus table
############################

all_genera <- unique(
  c(
    ancom$taxon,
    maas$feature,
    lefse$Genus,
    rf$Genus,
    network$Genus
  )
)

master <- data.frame(
  Genus = all_genera,
  stringsAsFactors = FALSE
) %>%
  left_join(
    da_summary,
    by = "Genus"
  ) %>%
  left_join(
    lefse_summary,
    by = "Genus"
  ) %>%
  left_join(
    rf_summary,
    by = "Genus"
  ) %>%
  left_join(
    network_summary,
    by = "Genus"
  ) %>%
  mutate(
    RobustAssociation =
      replace_na(
        RobustAssociation,
        FALSE
      ),

    DA_DirectionConflict =
      replace_na(
        DA_DirectionConflict,
        FALSE
      ),

    LEfSeSupport =
      replace_na(
        LEfSeSupport,
        FALSE
      ),

    RFStable =
      replace_na(
        RFStable,
        FALSE
      ),

    NetworkCentral =
      replace_na(
        NetworkCentral,
        FALSE
      ),

    # Number of independent supporting layers.
    # LEfSe counts as a support layer, but LEfSe-only taxa are
    # not retained as a final core class.
    EvidenceLayers =
      as.integer(RobustAssociation) +
      as.integer(LEfSeSupport) +
      as.integer(RFStable) +
      as.integer(NetworkCentral),

    CoreEvidence =
      RobustAssociation |
      RFStable |
      NetworkCentral,

    EvidenceClass =
      case_when(

        CoreEvidence &
          EvidenceLayers >= 2 ~
          "Multi-layer evidence",

        RobustAssociation ~
          "Robust association",

        RFStable ~
          "Stable discriminator",

        NetworkCentral ~
          "Network-central",

        TRUE ~
          NA_character_
      )
  )

############################
# 13. Candidate score / ranking
############################

# Continuous support score only for row ordering.
# It does NOT replace the evidence classes.
master <- master %>%
  mutate(
    RF_component =
      replace_na(
        SelectionFrequency,
        0
      ),

    Degree_component =
      ifelse(
        is.na(Degree),
        0,
        Degree / max(
          network$Degree,
          na.rm = TRUE
        )
      ),

    Betweenness_component =
      ifelse(
        is.na(Betweenness),
        0,
        Betweenness / max(
          network$Betweenness,
          na.rm = TRUE
        )
      ),

    EvidenceScore =
      2 * as.integer(RobustAssociation) +
      1 * as.integer(LEfSeSupport) +
      2 * as.integer(RFStable) +
      2 * as.integer(NetworkCentral) +
      RF_component +
      0.5 * Degree_component +
      0.5 * Betweenness_component
  )

core_table <- master %>%
  filter(
    !is.na(EvidenceClass)
  ) %>%
  arrange(
    desc(EvidenceLayers),
    desc(EvidenceScore)
  )

write.csv(
  master,
  file.path(
    outdir,
    "02_all_genera_integrated_evidence.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_table,
  file.path(
    outdir,
    "03_core_genera_evidence_classes.csv"
  ),
  row.names = FALSE
)

############################
# 14. Explicitly report key genera
############################

key_genera <- c(
  "Rhodococcus",
  "Leclercia",
  "Psychrobacter",
  "Campylobacter",
  "Mycobacterium"
)

key_report <- master %>%
  filter(
    Genus %in% key_genera
  ) %>%
  select(
    Genus,
    EvidenceClass,
    EvidenceLayers,
    RobustAssociation,
    DA_DirectionConflict,
    ConcordantContrasts,
    ConflictContrasts,
    LEfSeSupport,
    LEfSeComparisons,
    SelectionFrequency,
    RFStable,
    Degree,
    Betweenness,
    NetworkCentral
  ) %>%
  arrange(
    match(
      Genus,
      key_genera
    )
  )

write.csv(
  key_report,
  file.path(
    outdir,
    "04_key_genera_classification.csv"
  ),
  row.names = FALSE
)

cat("\n==============================================\n")
cat("KEY GENERA CLASSIFICATION\n")
cat("==============================================\n")
print(
  key_report,
  row.names = FALSE
)

############################
# 15. Build detailed bubble matrix data
############################

# ANCOM evidence: only robust contrasts are displayed.
bubble_ancom <- ancom_long %>%
  filter(
    Robust %in% TRUE
  ) %>%
  transmute(
    Genus,
    EvidenceColumn = paste0(
      "ANCOM-BC2\n",
      Contrast
    ),
    RawStrength = ANCOM_strength,
    EvidenceType = Direction,
    Source = "ANCOM-BC2"
  )

# MaAsLin2: q < 0.05 only.
bubble_maas <- maas_clade %>%
  filter(
    Significant %in% TRUE
  ) %>%
  transmute(
    Genus,
    EvidenceColumn = paste0(
      "MaAsLin2\n",
      Contrast
    ),
    RawStrength = Maaslin_strength,
    EvidenceType = Direction,
    Source = "MaAsLin2"
  )

# LEfSe.
bubble_lefse <- lefse %>%
  transmute(
    Genus,
    EvidenceColumn = paste0(
      "LEfSe\n",
      Comparison
    ),
    RawStrength = abs(scores),
    EvidenceType = paste0(
      EnrichedClade,
      " enriched"
    ),
    Source = "LEfSe"
  )

# RF.
bubble_rf <- rf %>%
  filter(
    SelectionFrequency >= rf_threshold
  ) %>%
  transmute(
    Genus,
    EvidenceColumn =
      "RF\nstability",
    RawStrength =
      SelectionFrequency,
    EvidenceType =
      "RF stability",
    Source =
      "Random forest"
  )

# Degree.
bubble_degree <- network %>%
  filter(
    Degree > 0
  ) %>%
  transmute(
    Genus,
    EvidenceColumn =
      "Network\ndegree",
    RawStrength =
      Degree,
    EvidenceType =
      "Network degree",
    Source =
      "SPIEC-EASI"
  )

# Betweenness.
bubble_between <- network %>%
  filter(
    Betweenness > 0
  ) %>%
  transmute(
    Genus,
    EvidenceColumn =
      "Network\nbetweenness",
    RawStrength =
      Betweenness,
    EvidenceType =
      "Network betweenness",
    Source =
      "SPIEC-EASI"
  )

bubble_all <- bind_rows(
  bubble_ancom,
  bubble_maas,
  bubble_lefse,
  bubble_rf,
  bubble_degree,
  bubble_between
)

############################
# 16. Select rows for main evidence matrix
############################

# Main figure:
# show the top 25 final core genera.
# You can change 25 to 15/20/30.
n_show <- 25

selected_genera <- core_table %>%
  slice_head(
    n = n_show
  ) %>%
  pull(
    Genus
  )

# Ensure the five user-highlighted genera are included whenever available.
selected_genera <- unique(
  c(
    selected_genera,
    intersect(
      key_genera,
      master$Genus
    )
  )
)

selected_master <- master %>%
  filter(
    Genus %in% selected_genera
  ) %>%
  arrange(
    EvidenceClass,
    EvidenceScore
  )

genus_levels <- selected_master$Genus

bubble_plot <- bubble_all %>%
  filter(
    Genus %in% selected_genera
  ) %>%
  group_by(
    EvidenceColumn
  ) %>%
  mutate(
    SizeScaled =
      if (
        length(unique(RawStrength)) > 1
      ) {
        scales::rescale(
          RawStrength,
          to = c(2.5, 9)
        )
      } else {
        rep(5.5, n())
      }
  ) %>%
  ungroup() %>%
  mutate(
    Genus = factor(
      Genus,
      levels = genus_levels
    )
  )

############################
# 17. Fixed evidence-column ordering
############################

column_order <- c(
  "ANCOM-BC2\nB vs A",
  "ANCOM-BC2\nC vs A",
  "ANCOM-BC2\nC vs B",
  "MaAsLin2\nB vs A",
  "MaAsLin2\nC vs A",
  "LEfSe\nA vs B",
  "LEfSe\nA vs C",
  "LEfSe\nB vs C",
  "RF\nstability",
  "Network\ndegree",
  "Network\nbetweenness"
)

bubble_plot$EvidenceColumn <- factor(
  bubble_plot$EvidenceColumn,
  levels = column_order
)

############################
# 18. Bubble evidence matrix
############################

p_matrix <- ggplot(
  bubble_plot,
  aes(
    x = EvidenceColumn,
    y = Genus
  )
) +

  geom_point(
    aes(
      size = SizeScaled,
      fill = EvidenceType
    ),
    shape = 21,
    colour = "black",
    stroke = 0.30,
    alpha = 0.92
  ) +

  scale_fill_manual(
    values = all_evidence_colours,
    drop = FALSE
  ) +

  scale_size_identity(
    guide = "none"
  ) +

  scale_x_discrete(
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    fill = "Evidence",
    title =
      "Integrated evidence for clade-associated microbial genera",
    subtitle =
      "Bubble size denotes within-analysis evidence strength"
  ) +

  theme_nature(
    base_size = 10.5
  ) +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    axis.text.y = element_text(
      face = "italic"
    ),

    panel.grid.major = element_line(
      colour = "grey92",
      linewidth = 0.30
    ),

    legend.position = "right"
  )

############################
# 19. Evidence-class annotation strip
############################

class_plot_data <- selected_master %>%
  mutate(
    Genus = factor(
      Genus,
      levels = genus_levels
    ),
    EvidenceClass = factor(
      EvidenceClass,
      levels = c(
        "Robust association",
        "Stable discriminator",
        "Network-central",
        "Multi-layer evidence"
      )
    )
  )

p_class <- ggplot(
  class_plot_data,
  aes(
    x = "Class",
    y = Genus,
    fill = EvidenceClass
  )
) +

  geom_tile(
    width = 0.85,
    height = 0.85,
    colour = "white",
    linewidth = 0.35
  ) +

  scale_fill_manual(
    values = class_colours,
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    fill = "Evidence class"
  ) +

  theme_void(
    base_size = 10.5
  ) +

  theme(
    axis.text.y = element_blank(),
    legend.position = "right"
  )

############################
# 20. Combined Nature-style figure
############################

p_combined <-
  p_class +
  p_matrix +
  plot_layout(
    widths = c(0.8, 9),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(
    outdir,
    "05_consensus_evidence_matrix.pdf"
  ),
  p_combined,
  width = 13,
  height = 8.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "05_consensus_evidence_matrix.png"
  ),
  p_combined,
  width = 13,
  height = 8.5,
  dpi = 600
)

############################
# 21. Simpler heatmap-style binary support matrix
############################

binary <- master %>%
  filter(
    Genus %in% selected_genera
  ) %>%
  transmute(
    Genus,
    `Robust association` =
      RobustAssociation,
    `LEfSe support` =
      LEfSeSupport,
    `RF stable` =
      RFStable,
    `Network central` =
      NetworkCentral
  ) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "Evidence",
    values_to = "Supported"
  ) %>%
  mutate(
    Genus = factor(
      Genus,
      levels = genus_levels
    )
  )

p_binary <- ggplot(
  binary,
  aes(
    Evidence,
    Genus,
    fill = Supported
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +

  scale_fill_manual(
    values = c(
      "FALSE" = "grey94",
      "TRUE" = "#3C5488"
    ),
    labels = c(
      "FALSE" = "No",
      "TRUE" = "Yes"
    )
  ) +

  labs(
    x = NULL,
    y = NULL,
    fill = "Supported",
    title =
      "Cross-method support for core microbial genera"
  ) +

  theme_nature() +

  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "italic"
    )
  )

ggsave(
  file.path(
    outdir,
    "06_binary_evidence_heatmap.pdf"
  ),
  p_binary,
  width = 7.5,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "06_binary_evidence_heatmap.png"
  ),
  p_binary,
  width = 7.5,
  height = 8,
  dpi = 600
)

############################
# 22. Evidence-class count figure
############################

class_counts <- core_table %>%
  count(
    EvidenceClass
  ) %>%
  mutate(
    EvidenceClass = factor(
      EvidenceClass,
      levels = c(
        "Robust association",
        "Stable discriminator",
        "Network-central",
        "Multi-layer evidence"
      )
    )
  )

p_class_count <- ggplot(
  class_counts,
  aes(
    EvidenceClass,
    n,
    fill = EvidenceClass
  )
) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.30
  ) +

  geom_text(
    aes(
      label = n
    ),
    vjust = -0.35,
    size = 4
  ) +

  scale_fill_manual(
    values = class_colours
  ) +

  labs(
    x = NULL,
    y = "Number of genera",
    title = "Integrated evidence classes"
  ) +

  theme_nature() +

  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

ggsave(
  file.path(
    outdir,
    "07_evidence_class_counts.pdf"
  ),
  p_class_count,
  width = 7,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "07_evidence_class_counts.png"
  ),
  p_class_count,
  width = 7,
  height = 4.8,
  dpi = 600
)

############################
# 23. Direction-conflict table
############################

conflicts <- master %>%
  filter(
    DA_DirectionConflict
  ) %>%
  arrange(
    desc(EvidenceScore)
  )

write.csv(
  conflicts,
  file.path(
    outdir,
    "08_ANCOM_MaAsLin2_direction_conflicts.csv"
  ),
  row.names = FALSE
)

############################
# 24. Session
############################

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)

cat("\n==============================================\n")
cat("Finished.\n")
cat("Output directory:", outdir, "\n")
cat("Degree centrality cutoff:", degree_cut, "\n")
cat("Betweenness cutoff:", betweenness_cut, "\n")
cat("RF stability cutoff:", rf_threshold, "\n")
cat("==============================================\n")

