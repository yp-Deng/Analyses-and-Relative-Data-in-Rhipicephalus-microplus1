############################################################
# Chromosomal virulence-factor analysis
# VFDB-style annotation workflow
#
# Designed for tables containing:
# Gene_id, Identity, E_value, VF_internal_id, VF_id, VF_name,
# Characteristics, Structure_features, Functions, Mechanism
#
# ANALYSES
#  1. VF annotation QC and best-hit filtering
#  2. Chromosomal virulence-gene inventory
#  3. Functional virulence-category classification
#  4. VF gene comparison with non-pathogenic relatives
#  5. Candidate-specific / shared / relative-specific VF analysis
#  6. VF prevalence comparison across genomes
#  7. Virulence-category enrichment analysis
#  8. Optional CDS-normalized enrichment models
#  9. Publication-ready figures and tables
#
# IMPORTANT
#  - All genomes used for comparison should be annotated with
#    the SAME VFDB version, search program, identity cutoff,
#    E-value cutoff and coverage criteria.
#  - "Virulence-associated homolog" is safer wording than
#    "virulence gene" when annotation is based only on homology.
#  - Enrichment does not establish pathogenicity.
############################################################

rm(list = ls())
set.seed(20260813)

############################
# 0. USER SETTINGS
############################

# Put all VF annotation files here.
# Supported extensions: .txt .tsv .csv
input_dir <- "VF_annotations"

# Genome metadata.
# Required:
#   Genome
#   Group
#
# Recommended:
#   Species
#   Total_CDS
#
# Example:
# Genome,Group,Species,Total_CDS
# Ano5,Candidate,Your_species,3150
# Relative1,Non_pathogenic,Closest_relative_1,3020
# Relative2,Non_pathogenic,Closest_relative_2,2980

metadata_file <- "genome_metadata.csv"

# Group labels used in the comparison.
candidate_label <- "Candidate"
relative_label  <- "Non_pathogenic"

# VF homology thresholds.
# Your uploaded table appears already filtered around >=40% identity;
# these are kept configurable.
min_identity <- 40
max_evalue   <- 1e-5

# Optional: require a named VF?
# FALSE is recommended because many valid VFDB hits have VF_name = "-".
require_named_vf <- FALSE

# Prevalence threshold for describing a VF as "candidate-core".
candidate_core_prevalence <- 0.50

# Stable comparison unit.
# VF_id is recommended because VF_name is missing for many hits.
comparison_feature <- "VF_id"

# Number of VFs shown in heatmaps/bubble plots.
top_n_vf <- 50

# Output directory.
outdir <- "Virulence_gene_analysis_results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


############################
# 1. PACKAGES
############################

pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "readr",
  "stringr", "purrr", "tibble", "scales",
  "broom", "forcats", "patchwork"
)

need <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(need) > 0) {
  install.packages(need)
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
  library(broom)
  library(forcats)
  library(patchwork)
})


############################
# 2. FIGURE STYLE
############################

nature_group <- c(
  "Candidate" = "#E64B35",
  "Non_pathogenic" = "#4DBBD5"
)

vf_class_colours <- c(
  "Candidate-specific" = "#E64B35",
  "Shared" = "#8491B4",
  "Relative-specific" = "#4DBBD5"
)

category_colours <- c(
  "Adhesion/Colonization" = "#E64B35",
  "Toxin/Cytolysis" = "#4DBBD5",
  "Secretion system/Effector" = "#00A087",
  "Iron/Nutrient acquisition" = "#3C5488",
  "Capsule/Surface polysaccharide" = "#F39B7F",
  "Motility/Flagella" = "#91D1C2",
  "Biofilm" = "#7E6148",
  "Immune evasion/Intracellular survival" = "#B09C85",
  "Stress resistance/Regulation" = "#8491B4",
  "Metabolism/Persistence" = "#DC0000",
  "Other/Unclassified" = "grey70"
)

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      plot.subtitle = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.45
      ),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_fig <- function(
    p,
    name,
    width = 7,
    height = 5
) {
  ggsave(
    file.path(outdir, paste0(name, ".pdf")),
    p,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    file.path(outdir, paste0(name, ".png")),
    p,
    width = width,
    height = height,
    dpi = 600
  )
}


############################################################
# PART I. INPUT AND QC
############################################################

############################
# 3. CREATE INPUT DIRECTORY
############################

if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE)

  stop(
    "\nCreated folder '", input_dir, "'.\n",
    "Place VF annotation .txt/.tsv/.csv files in this folder and rerun.\n"
  )
}


############################
# 4. FIND VF FILES
############################

vf_files <- list.files(
  input_dir,
  pattern = "\\.(txt|tsv|csv)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(vf_files) == 0) {
  stop(
    "No VF annotation files found in ",
    input_dir
  )
}


############################
# 5. FLEXIBLE TABLE READER
############################

read_vf_file <- function(path) {

  ext <- tolower(
    tools::file_ext(path)
  )

  if (ext == "csv") {

    x <- readr::read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    )

  } else {

    x <- readr::read_tsv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    )
  }

  x <- as.data.frame(
    x,
    check.names = FALSE
  )

  required <- c(
    "Gene_id",
    "Identity",
    "E_value",
    "VF_internal_id",
    "VF_id",
    "VF_name",
    "Characteristics",
    "Structure_features",
    "Functions",
    "Mechanism"
  )

  missing <- setdiff(
    required,
    names(x)
  )

  if (length(missing) > 0) {
    stop(
      basename(path),
      " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  # Infer genome from Gene_id:
  # Ano5_GM000050 -> Ano5
  inferred <- unique(
    sub(
      "_GM.*$",
      "",
      as.character(x$Gene_id)
    )
  )

  # If inference is ambiguous, use file basename.
  if (length(inferred) == 1) {
    genome <- inferred
  } else {
    genome <- tools::file_path_sans_ext(
      basename(path)
    )
  }

  x$Genome <- genome
  x$SourceFile <- basename(path)

  x
}


############################
# 6. READ ALL GENOMES
############################

vf_raw <- purrr::map_dfr(
  vf_files,
  read_vf_file
)

cat(
  "\nGenomes inferred from VF files:\n"
)

print(
  table(vf_raw$Genome)
)


############################
# 7. CLEAN DATA
############################

vf <- vf_raw %>%
  mutate(
    Gene_id = as.character(Gene_id),

    Identity = suppressWarnings(
      as.numeric(Identity)
    ),

    E_value = suppressWarnings(
      as.numeric(E_value)
    ),

    across(
      c(
        VF_internal_id,
        VF_id,
        VF_name,
        Characteristics,
        Structure_features,
        Functions,
        Mechanism
      ),
      ~ na_if(
        trimws(as.character(.x)),
        "-"
      )
    ),

    VF_label = ifelse(
      is.na(VF_name) |
        VF_name == "",
      VF_id,
      VF_name
    )
  )


############################
# 8. QC FILTERS
############################

vf_filtered <- vf %>%
  filter(
    !is.na(Identity),
    Identity >= min_identity,

    # E-value = 0 is valid and extremely significant.
    !is.na(E_value),
    E_value <= max_evalue,

    !is.na(VF_id),
    VF_id != ""
  )

if (require_named_vf) {
  vf_filtered <- vf_filtered %>%
    filter(
      !is.na(VF_name),
      VF_name != ""
    )
}


############################
# 9. BEST HIT PER CHROMOSOMAL GENE
############################
# If one chromosomal gene matches multiple VF entries:
# 1) lowest E-value
# 2) highest identity
#
# This avoids counting one CDS multiple times in total VF burden.
#
# A separate all-hits table is retained.

vf_best <- vf_filtered %>%
  arrange(
    Genome,
    Gene_id,
    E_value,
    desc(Identity)
  ) %>%
  group_by(
    Genome,
    Gene_id
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup()


############################
# 10. QC OUTPUT
############################

qc_summary <- vf_raw %>%
  group_by(Genome) %>%
  summarise(
    Raw_hits = n(),
    Raw_genes = n_distinct(Gene_id),
    .groups = "drop"
  ) %>%
  full_join(
    vf_filtered %>%
      group_by(Genome) %>%
      summarise(
        Filtered_hits = n(),
        Filtered_genes = n_distinct(Gene_id),
        Unique_VF_IDs = n_distinct(VF_id),
        Median_identity = median(
          Identity,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "Genome"
  )

write.csv(
  qc_summary,
  file.path(
    outdir,
    "01_VF_QC_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  vf_filtered,
  file.path(
    outdir,
    "01_all_filtered_VF_hits.csv"
  ),
  row.names = FALSE
)

write.csv(
  vf_best,
  file.path(
    outdir,
    "01_best_VF_hit_per_gene.csv"
  ),
  row.names = FALSE
)


############################################################
# PART II. FUNCTIONAL VIRULENCE CATEGORIES
############################################################

############################
# 11. CATEGORY CLASSIFIER
############################
# This is a reproducible keyword-based first-pass classifier.
# Always inspect the exported assignments.
#
# For a final manuscript, a manually curated mapping file is preferable
# for biologically important genes.

classify_vf_category <- function(
    vf_name,
    functions,
    characteristics,
    mechanism
) {

  txt <- tolower(
    paste(
      vf_name,
      functions,
      characteristics,
      mechanism,
      sep = " "
    )
  )

  case_when(

    str_detect(
      txt,
      "toxin|cytolysin|hemolysin|haemolysin|cytotox|phospholipase|protease.*toxin|rtx"
    ) ~
      "Toxin/Cytolysis",

    str_detect(
      txt,
      "type iv|type vi|type iii|t4ss|t6ss|t3ss|dot/icm|secreted effector|secretion system|effector"
    ) ~
      "Secretion system/Effector",

    str_detect(
      txt,
      "adhesi|attachment|colonization|colonisation|fibronectin|collagen binding|mscramm|invasion"
    ) ~
      "Adhesion/Colonization",

    str_detect(
      txt,
      "iron|feo|fbp|hitabc|siderophore|heme|haem|transferrin|nutrient acquisition|magnesium transporter|mgt"
    ) ~
      "Iron/Nutrient acquisition",

    str_detect(
      txt,
      "capsule|capsular|lps|lipopolysaccharide|los|lipooligosaccharide|alginate|surface polysaccharide"
    ) ~
      "Capsule/Surface polysaccharide",

    str_detect(
      txt,
      "flagell|motility|chemotaxis|swarming|swimming"
    ) ~
      "Motility/Flagella",

    str_detect(
      txt,
      "biofilm|pili|pilus|fimbr"
    ) ~
      "Biofilm",

    str_detect(
      txt,
      "immune evasion|macrophage|intracellular|phagosome|phagosomal|serum killing|complement|autophagy|persistent infection|persistence"
    ) ~
      "Immune evasion/Intracellular survival",

    str_detect(
      txt,
      "stress|oxidative|peroxide|katg|clpc|clpp|sigma|mprab|two-component|regulator|regulation|msrab"
    ) ~
      "Stress resistance/Regulation",

    str_detect(
      txt,
      "metabolism|biosynthesis|pantothen|pyrimidine|lipoate|lipid biosynthesis"
    ) ~
      "Metabolism/Persistence",

    TRUE ~
      "Other/Unclassified"
  )
}


vf_best <- vf_best %>%
  mutate(
    VF_Category =
      classify_vf_category(
        VF_name,
        Functions,
        Characteristics,
        Mechanism
      )
  )

write.csv(
  vf_best,
  file.path(
    outdir,
    "02_VF_functional_category_assignments.csv"
  ),
  row.names = FALSE
)


############################
# 12. MANUAL CATEGORY OVERRIDE
############################
# Optional file:
# VF_category_manual.csv
#
# VF_id,VF_Category_manual
# VF0356,Toxin/Cytolysis
# ...
#
# If present, manual assignments override keyword assignments.

manual_file <- "VF_category_manual.csv"

if (file.exists(manual_file)) {

  manual_cat <- read.csv(
    manual_file,
    stringsAsFactors = FALSE
  )

  if (
    all(
      c(
        "VF_id",
        "VF_Category_manual"
      ) %in%
        names(manual_cat)
    )
  ) {

    vf_best <- vf_best %>%
      left_join(
        manual_cat,
        by = "VF_id"
      ) %>%
      mutate(
        VF_Category =
          ifelse(
            !is.na(VF_Category_manual) &
              VF_Category_manual != "",
            VF_Category_manual,
            VF_Category
          )
      ) %>%
      select(
        -VF_Category_manual
      )
  }
}


############################################################
# PART III. SINGLE-GENOME VF INVENTORY
############################################################

############################
# 13. GENOME-LEVEL SUMMARY
############################

genome_vf_summary <- vf_best %>%
  group_by(Genome) %>%
  summarise(
    VF_associated_CDS =
      n_distinct(Gene_id),

    Unique_VF_IDs =
      n_distinct(VF_id),

    Named_VF_classes =
      n_distinct(
        VF_name[
          !is.na(VF_name)
        ]
      ),

    Median_identity =
      median(
        Identity,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write.csv(
  genome_vf_summary,
  file.path(
    outdir,
    "03_genome_VF_inventory_summary.csv"
  ),
  row.names = FALSE
)


############################
# 14. CATEGORY COUNTS
############################

category_counts <- vf_best %>%
  count(
    Genome,
    VF_Category,
    name = "VF_gene_count"
  )

write.csv(
  category_counts,
  file.path(
    outdir,
    "03_VF_category_counts.csv"
  ),
  row.names = FALSE
)


############################
# 15. VF NAME COUNTS
############################

vf_name_counts <- vf_best %>%
  count(
    Genome,
    VF_label,
    VF_id,
    VF_Category,
    name = "Gene_count"
  ) %>%
  arrange(
    Genome,
    desc(Gene_count)
  )

write.csv(
  vf_name_counts,
  file.path(
    outdir,
    "03_VF_gene_counts_by_name.csv"
  ),
  row.names = FALSE
)


############################
# 16. CATEGORY BAR PLOT
############################

p_cat <- ggplot(
  category_counts,
  aes(
    x = Genome,
    y = VF_gene_count,
    fill = VF_Category
  )
) +

  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.18
  ) +

  scale_fill_manual(
    values = category_colours,
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = "Virulence-associated CDS",
    fill = "VF category",
    title = "Chromosomal virulence-factor repertoire"
  ) +

  theme_nature() +

  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

save_fig(
  p_cat,
  "03_VF_category_composition",
  width = max(
    7,
    length(
      unique(
        category_counts$Genome
      )
    ) * 0.65
  ),
  height = 5.5
)


############################################################
# PART IV. GENOME METADATA
############################################################

############################
# 17. CREATE METADATA TEMPLATE IF ABSENT
############################

genome_names <- sort(
  unique(
    vf_best$Genome
  )
)

if (!file.exists(metadata_file)) {

  template <- data.frame(
    Genome = genome_names,
    Group = ifelse(
      seq_along(genome_names) == 1,
      candidate_label,
      relative_label
    ),
    Species = "",
    Total_CDS = NA_real_
  )

  write.csv(
    template,
    metadata_file,
    row.names = FALSE
  )

  warning(
    "\nCreated ", metadata_file, ". ",
    "Please verify Group labels and add Total_CDS values ",
    "before interpreting comparative enrichment analyses.\n"
  )
}

genome_meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

if (
  !"Genome" %in% names(genome_meta) ||
  !"Group" %in% names(genome_meta)
) {
  stop(
    "genome_metadata.csv must contain Genome and Group columns."
  )
}

missing_genomes <- setdiff(
  genome_names,
  genome_meta$Genome
)

if (length(missing_genomes) > 0) {

  stop(
    "These VF genomes are missing from genome_metadata.csv: ",
    paste(
      missing_genomes,
      collapse = ", "
    )
  )
}

genome_meta <- genome_meta %>%
  filter(
    Genome %in% genome_names
  )

if (
  !candidate_label %in% genome_meta$Group
) {
  stop(
    "No Candidate genome found in genome_metadata.csv."
  )
}

if (
  !relative_label %in% genome_meta$Group
) {

  warning(
    "No Non_pathogenic genomes are currently available. ",
    "Single-genome inventory will run, but comparative analyses will be skipped."
  )
}


############################################################
# PART V. VF GENE COMPARISON WITH NON-PATHOGENIC RELATIVES
############################################################

if (
  candidate_label %in% genome_meta$Group &&
  relative_label %in% genome_meta$Group
) {

  ############################
  # 18. FEATURE PRESENCE MATRIX
  ############################

  vf_presence <- vf_best %>%
    distinct(
      Genome,
      .data[[comparison_feature]]
    ) %>%
    mutate(
      Present = 1L
    ) %>%
    tidyr::complete(
      Genome =
        genome_meta$Genome,
      !!rlang::sym(comparison_feature),
      fill = list(
        Present = 0L
      )
    )

  presence_wide <- vf_presence %>%
    tidyr::pivot_wider(
      names_from = Genome,
      values_from = Present,
      values_fill = 0
    )

  write.csv(
    presence_wide,
    file.path(
      outdir,
      "04_VF_presence_absence_matrix.csv"
    ),
    row.names = FALSE
  )


  ############################
  # 19. GROUP PREVALENCE
  ############################

  vf_prevalence <- vf_presence %>%
    left_join(
      genome_meta %>%
        select(
          Genome,
          Group
        ),
      by = "Genome"
    ) %>%
    group_by(
      .data[[comparison_feature]],
      Group
    ) %>%
    summarise(
      N_genomes = n(),
      N_present = sum(Present),
      Prevalence = mean(Present),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = c(
        N_genomes,
        N_present,
        Prevalence
      ),
      values_fill = 0
    )

  # Standard column names.
  cand_prev_col <- paste0(
    "Prevalence_",
    candidate_label
  )

  rel_prev_col <- paste0(
    "Prevalence_",
    relative_label
  )

  vf_prevalence <- vf_prevalence %>%
    mutate(
      VF_comparison_class = case_when(

        .data[[cand_prev_col]] >=
          candidate_core_prevalence &
          .data[[rel_prev_col]] == 0 ~
          "Candidate-specific",

        .data[[cand_prev_col]] > 0 &
          .data[[rel_prev_col]] > 0 ~
          "Shared",

        .data[[cand_prev_col]] == 0 &
          .data[[rel_prev_col]] > 0 ~
          "Relative-specific",

        TRUE ~
          "Other"
      )
    )

  write.csv(
    vf_prevalence,
    file.path(
      outdir,
      "04_VF_candidate_vs_nonpathogenic_prevalence.csv"
    ),
    row.names = FALSE
  )


  ############################
  # 20. CANDIDATE-SPECIFIC VF TABLE
  ############################

  candidate_specific_ids <- vf_prevalence %>%
    filter(
      VF_comparison_class ==
        "Candidate-specific"
    ) %>%
    pull(
      .data[[comparison_feature]]
    )

  candidate_specific <- vf_best %>%
    filter(
      .data[[comparison_feature]] %in%
        candidate_specific_ids
    ) %>%
    left_join(
      genome_meta %>%
        select(
          Genome,
          Group
        ),
      by = "Genome"
    ) %>%
    filter(
      Group ==
        candidate_label
    )

  write.csv(
    candidate_specific,
    file.path(
      outdir,
      "04_candidate_specific_VF_homologs.csv"
    ),
    row.names = FALSE
  )


  ############################
  # 21. FISHER TEST PER VF ID
  ############################
  # Genome-level presence/absence test.
  #
  # This is preferable to treating each gene as an independent replicate.
  # With only 1 candidate genome, inferential power will be extremely low;
  # results should then be treated descriptively.

  group_sizes <- table(
    genome_meta$Group
  )

  candidate_n <- as.numeric(
    group_sizes[
      candidate_label
    ]
  )

  relative_n <- as.numeric(
    group_sizes[
      relative_label
    ]
  )

  fisher_vf <- vf_presence %>%
    left_join(
      genome_meta %>%
        select(
          Genome,
          Group
        ),
      by = "Genome"
    ) %>%
    filter(
      Group %in%
        c(
          candidate_label,
          relative_label
        )
    ) %>%
    group_by(
      .data[[comparison_feature]]
    ) %>%
    summarise(
      Candidate_present =
        sum(
          Present[
            Group == candidate_label
          ]
        ),

      Relative_present =
        sum(
          Present[
            Group == relative_label
          ]
        ),

      .groups = "drop"
    )

  fisher_results <- purrr::pmap_dfr(
    fisher_vf,
    function(
        ...,
        Candidate_present,
        Relative_present
    ) {

      vals <- list(...)

      a <- Candidate_present
      b <- candidate_n - a
      c <- Relative_present
      d <- relative_n - c

      ft <- fisher.test(
        matrix(
          c(
            a,
            b,
            c,
            d
          ),
          nrow = 2,
          byrow = TRUE
        )
      )

      data.frame(
        Candidate_present = a,
        Candidate_absent = b,
        Relative_present = c,
        Relative_absent = d,
        OddsRatio =
          unname(
            ft$estimate
          ),
        P_value =
          ft$p.value
      )
    }
  )

  # Reattach VF feature column explicitly.
  fisher_results[
    [
      comparison_feature
    ]
  ] <- fisher_vf[
    [
      comparison_feature
    ]
  ]

  fisher_results <- fisher_results %>%
    relocate(
      all_of(
        comparison_feature
      ),
      .before = 1
    ) %>%
    mutate(
      Q_value =
        p.adjust(
          P_value,
          method = "BH"
        ),

      Enriched_in_candidate =
        Q_value < 0.05 &
        OddsRatio > 1
    ) %>%
    arrange(
      Q_value,
      desc(
        OddsRatio
      )
    )

  write.csv(
    fisher_results,
    file.path(
      outdir,
      "05_VF_ID_genome_prevalence_Fisher_enrichment.csv"
    ),
    row.names = FALSE
  )


  ############################
  # 22. VF PREVALENCE BUBBLE PLOT
  ############################

  plot_prev <- vf_prevalence %>%
    filter(
      VF_comparison_class %in%
        c(
          "Candidate-specific",
          "Shared",
          "Relative-specific"
        )
    ) %>%
    arrange(
      desc(
        .data[[cand_prev_col]] -
          .data[[rel_prev_col]]
      )
    ) %>%
    slice_head(
      n = top_n_vf
    ) %>%
    select(
      all_of(
        comparison_feature
      ),
      all_of(
        cand_prev_col
      ),
      all_of(
        rel_prev_col
      ),
      VF_comparison_class
    ) %>%
    pivot_longer(
      cols =
        all_of(
          c(
            cand_prev_col,
            rel_prev_col
          )
        ),
      names_to = "Group",
      values_to = "Prevalence"
    ) %>%
    mutate(
      Group = sub(
        "^Prevalence_",
        "",
        Group
      )
    )

  vf_order <- plot_prev %>%
    filter(
      Group ==
        candidate_label
    ) %>%
    arrange(
      Prevalence
    ) %>%
    pull(
      .data[[comparison_feature]]
    )

  plot_prev[
    [
      comparison_feature
    ]
  ] <- factor(
    plot_prev[
      [
        comparison_feature
      ]
    ],
    levels = vf_order
  )

  p_prev <- ggplot(
    plot_prev,
    aes(
      x = Group,
      y = .data[[comparison_feature]]
    )
  ) +

    geom_point(
      aes(
        size = Prevalence,
        fill = Group
      ),
      shape = 21,
      colour = "black",
      stroke = 0.25
    ) +

    scale_fill_manual(
      values = nature_group
    ) +

    scale_size_continuous(
      range = c(
        1.5,
        8
      ),
      limits = c(
        0,
        1
      )
    ) +

    labs(
      x = NULL,
      y = NULL,
      size = "Genome\nprevalence",
      fill = "Genome group",
      title =
        "Virulence-factor comparison with non-pathogenic relatives"
    ) +

    theme_nature() +

    theme(
      axis.text.y =
        element_text(
          size = 8
        )
    )

  save_fig(
    p_prev,
    "05_VF_candidate_vs_nonpathogenic_bubble",
    width = 7,
    height = max(
      6,
      length(
        unique(
          plot_prev[
            [
              comparison_feature
            ]
          ]
        )
      ) * 0.16
    )
  )


  ############################################################
  # PART VI. VIRULENCE CATEGORY ENRICHMENT
  ############################################################

  ############################
  # 23. CATEGORY COUNTS PER GENOME
  ############################

  cat_genome <- vf_best %>%
    count(
      Genome,
      VF_Category,
      name = "Category_VF_CDS"
    ) %>%
    tidyr::complete(
      Genome =
        genome_meta$Genome,
      VF_Category,
      fill = list(
        Category_VF_CDS = 0
      )
    ) %>%
    left_join(
      genome_meta,
      by = "Genome"
    )

  write.csv(
    cat_genome,
    file.path(
      outdir,
      "06_VF_category_counts_per_genome.csv"
    ),
    row.names = FALSE
  )


  ############################
  # 24. AGGREGATED CATEGORY FISHER TEST
  ############################
  # Tests whether the composition of VF hits differs between groups.
  #
  # This is an exploratory compositional test.
  # It should NOT be interpreted as biological replication when only
  # one genome exists in a group.

  cat_group <- cat_genome %>%
    filter(
      Group %in%
        c(
          candidate_label,
          relative_label
        )
    ) %>%
    group_by(
      Group,
      VF_Category
    ) %>%
    summarise(
      Category_count =
        sum(
          Category_VF_CDS
        ),
      .groups = "drop"
    )

  group_total_vf <- cat_group %>%
    group_by(
      Group
    ) %>%
    summarise(
      Total_VF =
        sum(
          Category_count
        ),
      .groups = "drop"
    )

  category_levels <- sort(
    unique(
      cat_group$VF_Category
    )
  )

  category_fisher <- purrr::map_dfr(
    category_levels,
    function(cat) {

      cand_cat <- cat_group %>%
        filter(
          Group ==
            candidate_label,
          VF_Category ==
            cat
        ) %>%
        pull(
          Category_count
        )

      rel_cat <- cat_group %>%
        filter(
          Group ==
            relative_label,
          VF_Category ==
            cat
        ) %>%
        pull(
          Category_count
        )

      if (length(cand_cat) == 0) {
        cand_cat <- 0
      }

      if (length(rel_cat) == 0) {
        rel_cat <- 0
      }

      cand_total <- group_total_vf %>%
        filter(
          Group ==
            candidate_label
        ) %>%
        pull(
          Total_VF
        )

      rel_total <- group_total_vf %>%
        filter(
          Group ==
            relative_label
        ) %>%
        pull(
          Total_VF
        )

      tab <- matrix(
        c(
          cand_cat,
          cand_total - cand_cat,
          rel_cat,
          rel_total - rel_cat
        ),
        nrow = 2,
        byrow = TRUE
      )

      ft <- fisher.test(
        tab
      )

      data.frame(
        VF_Category = cat,

        Candidate_category_count =
          cand_cat,

        Candidate_other_VF =
          cand_total -
          cand_cat,

        Relative_category_count =
          rel_cat,

        Relative_other_VF =
          rel_total -
          rel_cat,

        OddsRatio =
          unname(
            ft$estimate
          ),

        P_value =
          ft$p.value
      )
    }
  ) %>%
    mutate(
      Q_value =
        p.adjust(
          P_value,
          method = "BH"
        ),

      Candidate_enriched =
        Q_value < 0.05 &
        OddsRatio > 1
    ) %>%
    arrange(
      Q_value
    )

  write.csv(
    category_fisher,
    file.path(
      outdir,
      "06_VF_category_compositional_enrichment_Fisher.csv"
    ),
    row.names = FALSE
  )


  ############################
  # 25. CDS-NORMALIZED VF BURDEN
  ############################

  if (
    "Total_CDS" %in%
      names(
        genome_meta
      )
  ) {

    genome_meta$Total_CDS <-
      suppressWarnings(
        as.numeric(
          genome_meta$Total_CDS
        )
      )

    if (
      all(
        !is.na(
          genome_meta$Total_CDS
        )
      ) &&
      all(
        genome_meta$Total_CDS > 0
      )
    ) {

      burden <- genome_vf_summary %>%
        left_join(
          genome_meta,
          by = "Genome"
        ) %>%
        mutate(
          VF_CDS_per_1000_CDS =
            1000 *
            VF_associated_CDS /
            Total_CDS
        )

      write.csv(
        burden,
        file.path(
          outdir,
          "07_VF_burden_normalized_by_total_CDS.csv"
        ),
        row.names = FALSE
      )


      p_burden <- ggplot(
        burden,
        aes(
          x = Group,
          y = VF_CDS_per_1000_CDS,
          fill = Group
        )
      ) +

        geom_boxplot(
          width = 0.50,
          outlier.shape = NA
        ) +

        geom_jitter(
          width = 0.08,
          size = 2.4,
          alpha = 0.85
        ) +

        scale_fill_manual(
          values =
            nature_group
        ) +

        labs(
          x = NULL,
          y =
            "Virulence-associated CDS\nper 1,000 CDS",
          fill = "Genome group",
          title =
            "Genome-size-normalized virulence-factor burden"
        ) +

        theme_nature()

      save_fig(
        p_burden,
        "07_VF_burden_normalized_by_CDS",
        width = 5.5,
        height = 4.8
      )


      ############################
      # 26. CATEGORY RATE MODELS
      ############################
      # Only run inferential GLMs when both groups contain >=3 genomes.
      #
      # Response: category-associated CDS count
      # Offset: total CDS
      #
      # Quasipoisson helps accommodate overdispersion.

      n_by_group <- table(
        genome_meta$Group
      )

      enough_replication <-
        all(
          n_by_group[
            c(
              candidate_label,
              relative_label
            )
          ] >= 3
        )

      if (enough_replication) {

        cat_glm_data <- cat_genome %>%
          left_join(
            genome_meta %>%
              select(
                Genome,
                Total_CDS
              ),
            by = "Genome",
            suffix = c(
              "",
              ".meta"
            )
          )

        if (
          "Total_CDS.meta" %in%
            names(
              cat_glm_data
            )
        ) {

          cat_glm_data$Total_CDS <-
            ifelse(
              is.na(
                cat_glm_data$Total_CDS
              ),
              cat_glm_data$Total_CDS.meta,
              cat_glm_data$Total_CDS
            )
        }

        glm_categories <- unique(
          cat_glm_data$VF_Category
        )

        category_glm <- purrr::map_dfr(
          glm_categories,
          function(cat) {

            dat <- cat_glm_data %>%
              filter(
                VF_Category ==
                  cat,
                Group %in%
                  c(
                    relative_label,
                    candidate_label
                  )
              )

            dat$Group <- factor(
              dat$Group,
              levels = c(
                relative_label,
                candidate_label
              )
            )

            fit <- glm(
              Category_VF_CDS ~
                Group +
                offset(
                  log(
                    Total_CDS
                  )
                ),
              data = dat,
              family =
                quasipoisson(
                  link = "log"
                )
            )

            tt <- broom::tidy(
              fit
            )

            target_term <-
              paste0(
                "Group",
                candidate_label
              )

            rr <- tt %>%
              filter(
                term ==
                  target_term
              )

            if (nrow(rr) == 0) {
              return(
                NULL
              )
            }

            data.frame(
              VF_Category = cat,
              Log_rate_ratio =
                rr$estimate,
              Rate_ratio =
                exp(
                  rr$estimate
                ),
              SE =
                rr$std.error,
              P_value =
                rr$p.value
            )
          }
        ) %>%
          mutate(
            Q_value =
              p.adjust(
                P_value,
                method = "BH"
              ),

            Candidate_enriched =
              Q_value < 0.05 &
              Rate_ratio > 1
          ) %>%
          arrange(
            Q_value
          )

        write.csv(
          category_glm,
          file.path(
            outdir,
            "07_VF_category_CDS_normalized_quasiPoisson_enrichment.csv"
          ),
          row.names = FALSE
        )

      } else {

        writeLines(
          paste0(
            "Category-level CDS-normalized GLM was skipped because ",
            "at least 3 genomes per group are recommended. ",
            "Current counts: ",
            paste(
              names(
                n_by_group
              ),
              n_by_group,
              sep = "=",
              collapse = "; "
            )
          ),
          file.path(
            outdir,
            "07_GLM_SKIPPED_insufficient_genome_replication.txt"
          )
        )
      }
    }
  }


  ############################################################
  # PART VII. FUNCTIONAL CATEGORY COMPARISON FIGURE
  ############################################################

  cat_plot <- cat_genome %>%
    filter(
      Group %in%
        c(
          candidate_label,
          relative_label
        )
    ) %>%
    group_by(
      Group,
      VF_Category
    ) %>%
    summarise(
      Mean_count =
        mean(
          Category_VF_CDS
        ),
      SE =
        ifelse(
          n() > 1,
          sd(
            Category_VF_CDS
          ) /
            sqrt(
              n()
            ),
          NA_real_
        ),
      .groups = "drop"
    )

  p_cat_compare <- ggplot(
    cat_plot,
    aes(
      x = VF_Category,
      y = Mean_count,
      fill = Group
    )
  ) +

    geom_col(
      position =
        position_dodge(
          width = 0.78
        ),
      width = 0.70,
      colour = "black",
      linewidth = 0.20
    ) +

    geom_errorbar(
      aes(
        ymin = pmax(
          0,
          Mean_count - SE
        ),
        ymax =
          Mean_count + SE
      ),
      position =
        position_dodge(
          width = 0.78
        ),
      width = 0.20,
      na.rm = TRUE
    ) +

    scale_fill_manual(
      values =
        nature_group
    ) +

    labs(
      x = NULL,
      y =
        "Mean number of virulence-associated CDS",
      fill = "Genome group",
      title =
        "Virulence-category comparison with non-pathogenic relatives"
    ) +

    theme_nature() +

    theme(
      axis.text.x =
        element_text(
          angle = 50,
          hjust = 1
        )
    )

  save_fig(
    p_cat_compare,
    "08_VF_category_candidate_vs_nonpathogenic",
    width = 9,
    height = 5.8
  )


  ############################################################
  # PART VIII. COMPARISON CLASS FIGURE
  ############################################################

  comparison_counts <- vf_prevalence %>%
    count(
      VF_comparison_class,
      name = "N_VF_IDs"
    ) %>%
    filter(
      VF_comparison_class %in%
        names(
          vf_class_colours
        )
    )

  p_compare_class <- ggplot(
    comparison_counts,
    aes(
      x = VF_comparison_class,
      y = N_VF_IDs,
      fill = VF_comparison_class
    )
  ) +

    geom_col(
      width = 0.62,
      colour = "black",
      linewidth = 0.25
    ) +

    geom_text(
      aes(
        label = N_VF_IDs
      ),
      vjust = -0.35,
      size = 4
    ) +

    scale_fill_manual(
      values =
        vf_class_colours
    ) +

    labs(
      x = NULL,
      y = "Number of VF IDs",
      title =
        "Candidate-specific and shared virulence-associated homologs"
    ) +

    theme_nature() +

    theme(
      legend.position =
        "none",
      axis.text.x =
        element_text(
          angle = 20,
          hjust = 1
        )
    )

  save_fig(
    p_compare_class,
    "08_VF_specific_shared_counts",
    width = 6.2,
    height = 4.8
  )
}


############################################################
# PART IX. HIGH-CONFIDENCE / BIOLOGICALLY INTERPRETABLE HITS
############################################################

############################
# 27. PRIORITIZATION TABLE
############################
# This does not claim pathogenicity.
#
# Priority is based on:
# - strong homology
# - interpretable VF name
# - major virulence-function category

major_categories <- c(
  "Adhesion/Colonization",
  "Toxin/Cytolysis",
  "Secretion system/Effector",
  "Iron/Nutrient acquisition",
  "Capsule/Surface polysaccharide",
  "Motility/Flagella",
  "Biofilm",
  "Immune evasion/Intracellular survival"
)

vf_priority <- vf_best %>%
  mutate(
    Strong_identity =
      Identity >= 60,

    Named =
      !is.na(VF_name) &
      VF_name != "",

    Major_category =
      VF_Category %in%
        major_categories,

    Priority_score =
      as.integer(
        Strong_identity
      ) +
      as.integer(
        Named
      ) +
      as.integer(
        Major_category
      )
  ) %>%
  arrange(
    Genome,
    desc(
      Priority_score
    ),
    desc(
      Identity
    ),
    E_value
  )

write.csv(
  vf_priority,
  file.path(
    outdir,
    "09_prioritized_virulence_associated_homologs.csv"
  ),
  row.names = FALSE
)


############################################################
# PART X. METHOD-SAFE SUMMARY
############################################################

############################
# 28. AUTOMATIC SUMMARY
############################

summary_lines <- c(
  paste0(
    "VF homology filters: identity >= ",
    min_identity,
    "%; E-value <= ",
    format(
      max_evalue,
      scientific = TRUE
    ),
    "."
  ),

  paste0(
    "Comparison unit: ",
    comparison_feature,
    "."
  ),

  "One best VFDB hit per chromosomal CDS was used for VF-burden counting.",

  paste0(
    "Candidate-core prevalence threshold: ",
    candidate_core_prevalence,
    "."
  ),

  "Candidate-specific means present in candidate genomes at/above the configured prevalence threshold and absent from all non-pathogenic relatives.",

  "Functional categories are keyword-derived unless overridden by VF_category_manual.csv.",

  "Genome-level VF presence/absence Fisher tests compare prevalence between genome groups.",

  "Category Fisher tests compare the composition of VF hits and should be treated as exploratory when genome replication is low.",

  "CDS-normalized quasi-Poisson enrichment models are run only when both groups contain at least three genomes.",

  "Virulence-associated homology alone does not establish a pathogenic phenotype."
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "00_analysis_notes.txt"
  )
)


############################################################
# PART XI. SESSION INFO
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)


############################################################
# PART XII. CONSOLE SUMMARY
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "Chromosomal virulence-factor analysis completed.\n"
)

cat(
  "Genomes:",
  paste(
    genome_names,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Filtered VF-associated CDS:",
  nrow(
    vf_best
  ),
  "\n"
)

cat(
  "Unique VF IDs:",
  n_distinct(
    vf_best$VF_id
  ),
  "\n"
)

cat(
  "Output directory:",
  normalizePath(
    outdir
  ),
  "\n"
)

cat(
  "\nImportant interpretation:\n",
  "- Prefer 'virulence-associated homolog' over 'virulence gene' ",
  "unless orthology/function is independently validated.\n",
  "- Use non-pathogenic relatives annotated with the same pipeline.\n",
  "- Treat enrichment as comparative genomic evidence, not proof of pathogenicity.\n"
)

cat(
  "====================================================\n"
)
