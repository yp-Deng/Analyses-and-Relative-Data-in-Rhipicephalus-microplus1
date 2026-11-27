############################################################
# Integrated taxonomic-functional microbiome analysis
#
# DATA LAYERS
#   Genus relative abundance
#   KO abundance / relative abundance
#   COG / eggNOG ortholog abundance
#   CAZy level 1 relative abundance
#   CAZy level 2 absolute abundance -> converted to relative abundance
#   Metadata
#
# CORE ANALYSES
#   A. Data QC / orientation / normalization
#   B. PCoA
#   C. Sequential + marginal PERMANOVA
#   D. PERMDISP
#   E. dbRDA + marginal / sequential tests + VIF
#   F. Variation partitioning:
#          Host = Clade + Stage
#          Environment = Temperature + Humidity + Altitude
#          Space = Longitude + Latitude
#   G. MaAsLin2 differential functional analysis
#          KO / COG / CAZy L1 / CAZy L2
#   H. KO -> KEGG pathway analysis and enrichment
#   I. Taxonomy-function concordance:
#          pairwise Mantel
#          Procrustes / PROTEST
#   J. Cross-layer synthesis:
#          marginal Clade R2 comparison
#          VPA pure-fraction comparison
#          Mantel correlation heatmap
#
# IMPORTANT
#   - Genus is used mainly as the taxonomic reference layer because
#     genus-level DA has already been analysed separately.
#   - CAZy level 2 is explicitly converted from absolute abundance
#     to sample-wise relative abundance.
#   - The uploaded COG table is actually a large count/abundance
#     matrix with features in rows and samples in columns; it is
#     transposed and normalized here.
#   - "Others" is excluded from COG composition by default.
#
############################################################

rm(list = ls())
set.seed(20260813)

############################################################
# 0. USER SETTINGS
############################################################

# The resolver accepts either user-friendly filenames or the current
# uploaded filenames. Rename your local files if desired.

resolve_file <- function(candidates) {

  hit <- candidates[
    file.exists(
      candidates
    )
  ]

  if (
    length(
      hit
    ) == 0
  ) {

    stop(
      "Cannot find any of these files:\n",
      paste(
        candidates,
        collapse = "\n"
      )
    )
  }

  hit[1]
}


genus_file <- "Genus.xlsx"
cazy_l1_file <- "CAZy_level1.xlsx"
cazy_l2_file <- "CAZy_level2_Absolute.xlsx"
cog_file <- "COG_table.txt"
metadata_file <- "metadata.xlsx"
ko_file <- "KO.xlsx"
outdir <- "Integrated_functional_microbiome_results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ----------------------------------------------------------
# Filtering thresholds
# ----------------------------------------------------------

# Multivariate beta diversity / dbRDA
prevalence_beta <- c(
  Genus = 0.10,
  KO = 0.10,
  COG = 0.10,
  CAZy_L1 = 0,
  CAZy_L2 = 0.10
)

# Differential functional analysis
prevalence_da <- c(
  KO = 0.20,
  COG = 0.20,
  CAZy_L1 = 0,
  CAZy_L2 = 0.20
)

mean_abundance_da <- c(
  KO = 1e-7,
  COG = 1e-6,
  CAZy_L1 = 0,
  CAZy_L2 = 1e-6
)


# Use only conventional COG IDs by default.
#
# "COG_only":
#    analyse COG0001 / COGxxxx entries only
#
# "all_orthologs":
#    include ENOG entries as well
#
# COG_only is recommended for a tractable and interpretable main analysis.
cog_feature_mode <- "COG_only"


# Remove "Others" from COG before relative normalization.
drop_cog_others <- TRUE


# Permutations
nperm <- 9999


# MaAsLin2 significance
fdr_cut <- 0.05


# Top features shown in plots
top_n_features <- 25


# KEGG mapping / pathway enrichment
run_kegg_mapping <- TRUE

remove_broad_kegg_pathways <- TRUE

broad_kegg_pathways <- c(
  "ko01100",
  "ko01110",
  "ko01120",
  "ko01200"
)


############################################################
# 1. PACKAGES
############################################################

cran_pkgs <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "vegan",
  "scales",
  "patchwork",
  "pheatmap",
  "purrr",
  "stringr"
)


need_cran <- cran_pkgs[
  !vapply(
    cran_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (
  length(
    need_cran
  ) > 0
) {

  install.packages(
    need_cran
  )
}


if (
  !requireNamespace(
    "BiocManager",
    quietly = TRUE
  )
) {

  install.packages(
    "BiocManager"
  )
}


bioc_pkgs <- c(
  "Maaslin2",
  "KEGGREST"
)


need_bioc <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (
  length(
    need_bioc
  ) > 0
) {

  BiocManager::install(
    need_bioc,
    ask = FALSE,
    update = FALSE
  )
}


suppressPackageStartupMessages({

  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(scales)
  library(patchwork)
  library(pheatmap)
  library(purrr)
  library(stringr)

})


############################################################
# 2. FIGURE STYLE
############################################################

nature_clade <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)


layer_colours <- c(
  "Genus" = "#E64B35",
  "KO" = "#4DBBD5",
  "COG" = "#00A087",
  "CAZy_L1" = "#3C5488",
  "CAZy_L2" = "#F39B7F"
)


vpa_colours <- c(
  "Host" = "#3C5488",
  "Environment" = "#E64B35",
  "Space" = "#00A087"
)


theme_nature <- function(
    base_size = 11
) {

  theme_classic(
    base_size = base_size
  ) +

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
          face = "bold",
          colour = "black"
        ),

      plot.subtitle =
        element_text(
          colour = "black"
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
        ),

      panel.grid.minor =
        element_blank()
    )
}


save_fig <- function(
    plot,
    name,
    width = 7,
    height = 5
) {

  ggsave(
    file.path(
      outdir,
      paste0(
        name,
        ".pdf"
      )
    ),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    file.path(
      outdir,
      paste0(
        name,
        ".png"
      )
    ),
    plot,
    width = width,
    height = height,
    dpi = 600
  )
}


############################################################
# 3. GENERIC INPUT FUNCTIONS
############################################################

read_table_flexible <- function(
    path
) {

  # Some files have an .xlsx suffix but are actually tab-delimited text.
  # Test whether readxl can open it; if not, use read_tsv.

  x <- tryCatch(

    readxl::read_excel(
      path
    ),

    error = function(e) {

      readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE
      )
    }
  )

  as.data.frame(
    x,
    check.names = FALSE
  )
}


to_numeric_matrix <- function(
    df
) {

  x <- df %>%

    mutate(
      across(
        everything(),
        ~ suppressWarnings(
          as.numeric(
            .x
          )
        )
      )
    ) %>%

    as.data.frame()

  x
}


normalize_rows <- function(
    mat
) {

  mat <- as.matrix(
    mat
  )

  totals <- rowSums(
    mat,
    na.rm = TRUE
  )

  if (
    any(
      totals <= 0
    )
  ) {

    stop(
      "At least one sample has zero total abundance."
    )
  }

  sweep(
    mat,
    1,
    totals,
    "/"
  )
}


remove_zero_features <- function(
    mat
) {

  mat[
    ,
    colSums(
      mat,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]
}


############################################################
# 4. READ METADATA
############################################################

meta <- read_table_flexible(
  metadata_file
)


rename_if_present <- function(
    df,
    old,
    new
) {

  if (
    old %in%
      names(
        df
      ) &&
      !new %in%
        names(
          df
        )
  ) {

    names(
      df
    )[
      names(
        df
      ) ==
        old
    ] <- new
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
  names(
    meta
  )
)


if (
  length(
    missing_meta
  ) >
    0
) {

  stop(
    "Missing metadata columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}


meta <- meta %>%

  mutate(

    SampleID =
      as.character(
        SampleID
      ),

    Clade =
      factor(
        Clade,
        levels = c(
          "A",
          "B",
          "C"
        )
      ),

    Stage =
      factor(
        Stage
      ),

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


if (
  anyDuplicated(
    meta$SampleID
  ) > 0
) {

  stop(
    "Duplicated SampleID in metadata."
  )
}


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


############################################################
# 5. READ GENUS
############################################################

genus_raw <- read_table_flexible(
  genus_file
)


if (
  !"SampleID" %in%
    names(
      genus_raw
    )
) {

  stop(
    "Genus table must contain SampleID."
  )
}


genus_raw$SampleID <-
  as.character(
    genus_raw$SampleID
  )


genus_mat <- genus_raw %>%

  select(
    -SampleID
  ) %>%

  to_numeric_matrix()


rownames(
  genus_mat
) <- genus_raw$SampleID


genus_mat[
  is.na(
    genus_mat
  )
] <- 0


genus_mat <- remove_zero_features(
  genus_mat
)


# Use a within-genus composition for distance analyses.
genus_rel <- normalize_rows(
  genus_mat
)


############################################################
# 6. READ KO
############################################################

ko_raw <- read_table_flexible(
  ko_file
)


if (
  !"SampleID" %in%
    names(
      ko_raw
    )
) {

  stop(
    "KO table must contain SampleID."
  )
}


ko_raw$SampleID <-
  as.character(
    ko_raw$SampleID
  )


ko_cols <- grep(
  "^K[0-9]{5}$",
  names(
    ko_raw
  ),
  value = TRUE
)


if (
  length(
    ko_cols
  ) ==
    0
) {

  stop(
    "No Kxxxxx columns detected in KO table."
  )
}


ko_original <- ko_raw %>%

  select(
    all_of(
      ko_cols
    )
  ) %>%

  to_numeric_matrix()


rownames(
  ko_original
) <- ko_raw$SampleID


ko_original[
  is.na(
    ko_original
  )
] <- 0


ko_original <- remove_zero_features(
  ko_original
)


# For functional composition / distances.
ko_rel <- normalize_rows(
  ko_original
)


############################################################
# 7. READ CAZy LEVEL 1
############################################################

cazy1_raw <- read_table_flexible(
  cazy_l1_file
)


# Uploaded format:
# CAZy_Modules | GH | GT | CBM | CE | PL | AA

sample_col_cazy1 <- names(
  cazy1_raw
)[1]


names(
  cazy1_raw
)[1] <- "SampleID"


cazy1_raw$SampleID <-
  as.character(
    cazy1_raw$SampleID
  )


cazy1_original <- cazy1_raw %>%

  select(
    -SampleID
  ) %>%

  to_numeric_matrix()


rownames(
  cazy1_original
) <- cazy1_raw$SampleID


cazy1_original[
  is.na(
    cazy1_original
  )
] <- 0


cazy1_original <- remove_zero_features(
  cazy1_original
)


# Preserve the original summed CAZyme relative signal for descriptive QC.
cazy1_total_signal <- rowSums(
  cazy1_original
)


# For composition-level analyses, normalize within the six CAZy classes.
cazy1_rel <- normalize_rows(
  cazy1_original
)


############################################################
# 8. READ CAZy LEVEL 2
############################################################

cazy2_raw <- read_table_flexible(
  cazy_l2_file
)


feature_col_cazy2 <- names(cazy2_raw)[1]

cazy2_features <- as.character(
  cazy2_raw[[feature_col_cazy2]]
)


cazy2_sample_cols <- intersect(
  names(cazy2_raw),
  meta$SampleID
)

cat("\nCAZy level 2:\n", "Detected sample columns:",length(cazy2_sample_cols), "\n")


missing_cazy2_samples <- setdiff(
  meta$SampleID,
  cazy2_sample_cols
)

if (length(missing_cazy2_samples) > 0) {
  
  cat(
    "Samples in metadata but absent from CAZy L2:\n"
  )
  
  print(missing_cazy2_samples)
}

extra_cazy2_cols <- setdiff(
  names(cazy2_raw),
  c(
    feature_col_cazy2,
    cazy2_sample_cols
  )
)

if (length(extra_cazy2_cols) > 0) {
  cat(
    "Non-sample columns ignored in CAZy L2:\n"
  )
  print(extra_cazy2_cols)
}

if (length(cazy2_sample_cols) == 0) {
  stop(
    "No CAZy level 2 sample columns match metadata$SampleID."
  )
}


cazy2_values <- cazy2_raw %>%
  
  select(
    all_of(cazy2_sample_cols)
  ) %>%
  
  to_numeric_matrix()

rownames(cazy2_values) <- make.unique(
  cazy2_features
)

# ----------------------------------------------------------
# 8.4 CAZy family × sample
#          ↓
#     sample × CAZy family
# ----------------------------------------------------------

cazy2_abs <- t(as.matrix(cazy2_values))

storage.mode(cazy2_abs) <- "numeric"

cazy2_abs[
  is.na(cazy2_abs)
] <- 0

cazy2_abs <- remove_zero_features(
  cazy2_abs
)

cazy2_sample_total <- rowSums(
  cazy2_abs,
  na.rm = TRUE
)

zero_cazy2_samples <- names(
  cazy2_sample_total
)[
  cazy2_sample_total <= 0
]

if (length(zero_cazy2_samples) > 0) {
  
  cat(
    "\nWARNING: CAZy level 2 contains samples with zero total abundance:\n"
  )
  
  print(zero_cazy2_samples)
  

  cazy2_abs <- cazy2_abs[
    cazy2_sample_total > 0,
    ,
    drop = FALSE
  ]
}

if (any(rowSums(cazy2_abs) <= 0)) {
  
  stop(
    "CAZy L2 still contains zero-total samples after filtering."
  )
}

# ----------------------------------------------------------
# 8.7 absolute abundance -> relative abundance
# ----------------------------------------------------------

cazy2_rel <- sweep(cazy2_abs, 1, rowSums(cazy2_abs), "/")


cat("\nCAZy L2 relative abundance row sums:\n")

print(summary(rowSums(cazy2_rel)))

cat("\nCAZy L2 dimensions:\n")

print(dim(cazy2_rel))


stopifnot(
  all(
    abs(
      rowSums(cazy2_rel) - 1
    ) < 1e-8
  )
)


write.csv(
  data.frame(
    SampleID = rownames(cazy2_rel),
    cazy2_rel,
    check.names = FALSE
  ),
  file.path(
    outdir,
    "00_CAZy_level2_converted_relative_abundance.csv"
  ),
  row.names = FALSE
)


############################################################
# 9. READ COG / EGGNOG
############################################################

cog_raw <- readr::read_tsv(
  cog_file,
  show_col_types = FALSE,
  progress = FALSE
)


cog_raw <- as.data.frame(
  cog_raw,
  check.names = FALSE
)


# Uploaded format:
# first column = ortholog / COG ID
# remaining columns = samples

feature_col_cog <- names(
  cog_raw
)[1]


cog_features <- as.character(cog_raw[[feature_col_cog]])


cog_values <- cog_raw %>%
  select(
    -all_of(
      feature_col_cog
    )
  ) %>%
  to_numeric_matrix()


rownames(
  cog_values
) <- cog_features


if (
  drop_cog_others
) {

  keep_cog_rows <-
    tolower(
      rownames(
        cog_values
      )
    ) !=
    "others"

  cog_values <-
    cog_values[
      keep_cog_rows,
      ,
      drop = FALSE
    ]
}


if (
  cog_feature_mode ==
    "COG_only"
) {

  keep_cog <-
    grepl(
      "^COG[0-9]+$",
      rownames(
        cog_values
      ),
      ignore.case = FALSE
    )

  cog_values <-
    cog_values[
      keep_cog,
      ,
      drop = FALSE
    ]
}


# transpose to sample x feature
cog_abs <- t(
  as.matrix(
    cog_values
  )
)


storage.mode(
  cog_abs
) <-
  "numeric"


cog_abs[
  is.na(
    cog_abs
  )
] <- 0


cog_abs <- remove_zero_features(
  cog_abs
)


# Convert to within-COG relative composition.
cog_rel <- normalize_rows(
  cog_abs
)


write.csv(
  data.frame(
    SampleID =
      rownames(
        cog_rel
      ),
    cog_rel,
    check.names = FALSE
  ),
  file.path(
    outdir,
    "00_COG_converted_relative_abundance.csv"
  ),
  row.names = FALSE
)


############################################################
# 10. ALIGN ALL DATA LAYERS TO COMMON SAMPLES
############################################################

sample_sets <- list(
  Metadata =
    meta$SampleID,

  Genus =
    rownames(
      genus_rel
    ),

  KO =
    rownames(
      ko_rel
    ),

  COG =
    rownames(
      cog_rel
    ),

  CAZy_L1 =
    rownames(
      cazy1_rel
    ),

  CAZy_L2 =
    rownames(
      cazy2_rel
    )
)


common_samples <- Reduce(
  intersect,
  sample_sets
)


if (
  length(
    common_samples
  ) <
    10
) {

  stop(
    "Fewer than 10 samples are shared across all layers."
  )
}


# Preserve metadata order.
common_samples <- meta$SampleID[
  meta$SampleID %in%
    common_samples
]


meta <- meta %>%

  filter(
    SampleID %in%
      common_samples
  )


genus_rel <-
  genus_rel[
    common_samples,
    ,
    drop = FALSE
  ]


ko_rel <-
  ko_rel[
    common_samples,
    ,
    drop = FALSE
  ]


ko_original <-
  ko_original[
    common_samples,
    ,
    drop = FALSE
  ]


cog_rel <-
  cog_rel[
    common_samples,
    ,
    drop = FALSE
  ]


cog_abs <-
  cog_abs[
    common_samples,
    ,
    drop = FALSE
  ]


cazy1_rel <-
  cazy1_rel[
    common_samples,
    ,
    drop = FALSE
  ]


cazy1_original <-
  cazy1_original[
    common_samples,
    ,
    drop = FALSE
  ]


cazy2_rel <-
  cazy2_rel[
    common_samples,
    ,
    drop = FALSE
  ]


cazy2_abs <-
  cazy2_abs[
    common_samples,
    ,
    drop = FALSE
  ]


stopifnot(
  all(
    meta$SampleID ==
      rownames(
        genus_rel
      )
  ),

  all(
    meta$SampleID ==
      rownames(
        ko_rel
      )
  ),

  all(
    meta$SampleID ==
      rownames(
        cog_rel
      )
  ),

  all(
    meta$SampleID ==
      rownames(
        cazy1_rel
      )
  ),

  all(
    meta$SampleID ==
      rownames(
        cazy2_rel
      )
  )
)


############################################################
# 11. QC SUMMARY
############################################################

qc_layers <- data.frame(

  Layer =
    c(
      "Genus",
      "KO",
      "COG",
      "CAZy_L1",
      "CAZy_L2"
    ),

  Samples =
    c(
      nrow(
        genus_rel
      ),
      nrow(
        ko_rel
      ),
      nrow(
        cog_rel
      ),
      nrow(
        cazy1_rel
      ),
      nrow(
        cazy2_rel
      )
    ),

  Features =
    c(
      ncol(
        genus_rel
      ),
      ncol(
        ko_rel
      ),
      ncol(
        cog_rel
      ),
      ncol(
        cazy1_rel
      ),
      ncol(
        cazy2_rel
      )
    )
)


write.csv(
  qc_layers,
  file.path(
    outdir,
    "00_integrated_layer_QC.csv"
  ),
  row.names = FALSE
)


print(
  qc_layers
)


############################################################
# 12. PREPARE ANALYSIS LAYERS
############################################################

layer_mats <- list(
  Genus = genus_rel,
  KO = ko_rel,
  COG = cog_rel,
  CAZy_L1 = cazy1_rel,
  CAZy_L2 = cazy2_rel
)


filter_beta_matrix <- function(
    mat,
    prevalence_cut
) {

  prev <- colMeans(
    mat > 0
  )

  keep <-
    prev >=
    prevalence_cut

  if (
    sum(
      keep
    ) <
      2
  ) {

    warning(
      "Filtering retained fewer than 2 features; using all non-zero features."
    )

    return(
      mat
    )
  }

  mat[
    ,
    keep,
    drop = FALSE
  ]
}


layer_beta <- list()


for (
  nm in
    names(
      layer_mats
    )
) {

  layer_beta[[nm]] <- filter_beta_matrix(layer_mats[[nm]],prevalence_beta[[nm]]  )
}


############################################################
# PART I. GENERIC MULTIVARIATE ANALYSIS
############################################################

############################################################
# 13. FUNCTION FOR PCOA / PERMANOVA / PERMDISP / dbRDA / VPA
############################################################

run_multivariate_layer <- function(
    mat,
    layer_name,
    meta
) {

  cat(
    "\n=================================================\n"
  )

  cat(
    "MULTIVARIATE LAYER:",
    layer_name,
    "\n"
  )

  cat(
    "Features:",
    ncol(
      mat
    ),
    "\n"
  )


  layer_dir <- file.path(
    outdir,
    paste0(
      "01_multivariate_",
      layer_name
    )
  )


  dir.create(
    layer_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # --------------------------------------------------------
  # Bray-Curtis
  # --------------------------------------------------------

  bray <- vegan::vegdist(
    mat,
    method = "bray"
  )


  # --------------------------------------------------------
  # PCoA
  # --------------------------------------------------------

  pco <- cmdscale(
    bray,
    eig = TRUE,
    k = 2,
    add = TRUE
  )


  pco_df <- data.frame(

    SampleID =
      rownames(
        pco$points
      ),

    Axis1 =
      pco$points[
        ,
        1
      ],

    Axis2 =
      pco$points[
        ,
        2
      ]
  ) %>%

    left_join(
      meta,
      by = "SampleID"
    )


  positive_eig <-
    pco$eig[
      pco$eig > 0
    ]


  pc1 <-
    100 *
    positive_eig[1] /
    sum(
      positive_eig
    )


  pc2 <-
    100 *
    positive_eig[2] /
    sum(
      positive_eig
    )


  p_pcoa <- ggplot(
    pco_df,
    aes(
      Axis1,
      Axis2,
      colour = Clade
    )
  ) +

    stat_ellipse(
      aes(
        fill = Clade
      ),
      geom = "polygon",
      alpha = 0.10,
      colour = NA,
      level = 0.95
    ) +

    geom_point(
      size = 2.8,
      alpha = 0.9
    ) +

    scale_colour_manual(
      values =
        nature_clade
    ) +

    scale_fill_manual(
      values =
        nature_clade
    ) +

    labs(

      x =
        paste0(
          "PCoA1 (",
          round(
            pc1,
            1
          ),
          "%)"
        ),

      y =
        paste0(
          "PCoA2 (",
          round(
            pc2,
            1
          ),
          "%)"
        ),

      colour = "Clade",
      fill = "Clade",

      title =
        paste0(
          layer_name,
          " Bray-Curtis PCoA"
        )
    ) +

    theme_nature()


  ggsave(
    file.path(
      layer_dir,
      "PCoA.pdf"
    ),
    p_pcoa,
    width = 6.5,
    height = 5.2,
    device = cairo_pdf
  )


  ggsave(
    file.path(
      layer_dir,
      "PCoA.png"
    ),
    p_pcoa,
    width = 6.5,
    height = 5.2,
    dpi = 600
  )


  # --------------------------------------------------------
  # Sequential PERMANOVA
  # --------------------------------------------------------

  set.seed(
    20260813
  )


  perm_seq <- vegan::adonis2(

    bray ~
      Clade +
      Stage +
      Temperature_z +
      Humidity_z +
      Altitude_z +
      Latitude_z +
      Longitude_z,

    data =
      meta,

    permutations =
      nperm,

    by =
      "terms"
  )


  # --------------------------------------------------------
  # Marginal PERMANOVA
  # --------------------------------------------------------

  set.seed(
    20260813
  )


  perm_margin <- vegan::adonis2(

    bray ~
      Clade +
      Stage +
      Temperature_z +
      Humidity_z +
      Altitude_z +
      Latitude_z +
      Longitude_z,

    data =
      meta,

    permutations =
      nperm,

    by =
      "margin"
  )


  write.csv(
    as.data.frame(
      perm_seq
    ),
    file.path(
      layer_dir,
      "PERMANOVA_sequential.csv"
    )
  )


  write.csv(
    as.data.frame(
      perm_margin
    ),
    file.path(
      layer_dir,
      "PERMANOVA_marginal.csv"
    )
  )


  # --------------------------------------------------------
  # PERMDISP
  # --------------------------------------------------------

  disp <- vegan::betadisper(
    bray,
    meta$Clade,
    type = "median"
  )


  set.seed(
    20260813
  )


  disp_test <- vegan::permutest(
    disp,
    permutations =
      nperm
  )


  capture.output(
    disp_test,
    file = file.path(
      layer_dir,
      "PERMDISP_Clade.txt"
    )
  )


  # --------------------------------------------------------
  # dbRDA
  # --------------------------------------------------------

  db <- vegan::dbrda(

    bray ~
      Clade +
      Stage +
      Temperature_z +
      Humidity_z +
      Altitude_z +
      Latitude_z +
      Longitude_z,

    data =
      meta,

    sqrt.dist =
      TRUE
  )


  set.seed(
    20260813
  )


  db_overall <- anova(
    db,
    permutations =
      nperm
  )


  set.seed(
    20260813
  )


  db_seq <- anova(
    db,
    by = "terms",
    permutations =
      nperm
  )


  set.seed(
    20260813
  )


  db_margin <- anova(
    db,
    by = "margin",
    permutations =
      nperm
  )


  write.csv(
    as.data.frame(
      db_overall
    ),
    file.path(
      layer_dir,
      "dbRDA_overall.csv"
    )
  )


  write.csv(
    as.data.frame(
      db_seq
    ),
    file.path(
      layer_dir,
      "dbRDA_sequential.csv"
    )
  )


  write.csv(
    as.data.frame(
      db_margin
    ),
    file.path(
      layer_dir,
      "dbRDA_marginal.csv"
    )
  )


  # --------------------------------------------------------
  # VIF
  # --------------------------------------------------------

  vif_values <- vegan::vif.cca(
    db
  )


  write.csv(

    data.frame(
      Variable =
        names(
          vif_values
        ),

      VIF =
        as.numeric(
          vif_values
        )
    ),

    file.path(
      layer_dir,
      "dbRDA_VIF.csv"
    ),

    row.names = FALSE
  )


  # --------------------------------------------------------
  # dbRDA plot
  # --------------------------------------------------------

  db_sites <- as.data.frame(
    vegan::scores(
      db,
      display = "sites",
      choices = 1:2,
      scaling = 2
    )
  )


  names(
    db_sites
  )[1:2] <- c(
    "dbRDA1",
    "dbRDA2"
  )


  db_sites$SampleID <-
    rownames(
      db_sites
    )


  db_sites <- db_sites %>%

    left_join(
      meta,
      by = "SampleID"
    )


  db_bp <- tryCatch(

    as.data.frame(
      vegan::scores(
        db,
        display = "bp",
        choices = 1:2,
        scaling = 2
      )
    ),

    error =
      function(e) {
        data.frame()
      }
  )


  if (
    nrow(
      db_bp
    ) > 0
  ) {

    names(
      db_bp
    )[1:2] <-
      c(
        "dbRDA1",
        "dbRDA2"
      )

    db_bp$Variable <-
      rownames(
        db_bp
      )


    db_bp <- db_bp %>%

      filter(
        Variable %in%
          c(
            "Temperature_z",
            "Humidity_z",
            "Altitude_z",
            "Latitude_z",
            "Longitude_z"
          )
      ) %>%

      mutate(
        Variable =
          recode(
            Variable,
            Temperature_z =
              "Temperature",
            Humidity_z =
              "Humidity",
            Altitude_z =
              "Altitude",
            Latitude_z =
              "Latitude",
            Longitude_z =
              "Longitude"
          )
      )
  }


  p_db <- ggplot(
    db_sites,
    aes(
      dbRDA1,
      dbRDA2,
      colour = Clade
    )
  ) +

    stat_ellipse(
      aes(
        fill = Clade
      ),
      geom = "polygon",
      alpha = 0.10,
      colour = NA,
      level = 0.95
    ) +

    geom_point(
      size = 2.8,
      alpha = 0.9
    )


  if (
    nrow(
      db_bp
    ) > 0
  ) {

    site_range <- max(
      abs(
        c(
          db_sites$dbRDA1,
          db_sites$dbRDA2
        )
      ),
      na.rm = TRUE
    )


    bp_range <- max(
      abs(
        c(
          db_bp$dbRDA1,
          db_bp$dbRDA2
        )
      ),
      na.rm = TRUE
    )


    mult <-
      ifelse(
        bp_range > 0,
        0.65 *
          site_range /
          bp_range,
        1
      )


    p_db <- p_db +

      geom_segment(
        data =
          db_bp,
        aes(
          x = 0,
          y = 0,
          xend =
            dbRDA1 *
            mult,
          yend =
            dbRDA2 *
            mult
        ),
        inherit.aes = FALSE,
        arrow =
          arrow(
            length =
              unit(
                0.18,
                "cm"
              )
          ),
        colour = "black",
        linewidth = 0.6
      ) +

      ggrepel::geom_text_repel(
        data =
          db_bp,
        aes(
          x =
            dbRDA1 *
            mult *
            1.1,
          y =
            dbRDA2 *
            mult *
            1.1,
          label =
            Variable
        ),
        inherit.aes = FALSE,
        colour = "black",
        size = 3.3
      )
  }


  p_db <- p_db +

    scale_colour_manual(
      values =
        nature_clade
    ) +

    scale_fill_manual(
      values =
        nature_clade
    ) +

    labs(
      x = "dbRDA1",
      y = "dbRDA2",
      colour = "Clade",
      fill = "Clade",
      title =
        paste0(
          layer_name,
          " dbRDA"
        )
    ) +

    theme_nature()


  ggsave(
    file.path(
      layer_dir,
      "dbRDA.pdf"
    ),
    p_db,
    width = 7,
    height = 5.5,
    device = cairo_pdf
  )


  ggsave(
    file.path(
      layer_dir,
      "dbRDA.png"
    ),
    p_db,
    width = 7,
    height = 5.5,
    dpi = 600
  )


  # --------------------------------------------------------
  # Variation partitioning
  # --------------------------------------------------------

  host_mm <- model.matrix(
    ~ Clade + Stage,
    data = meta
  )[
    ,
    -1,
    drop = FALSE
  ]


  env_mm <- as.matrix(
    meta[
      ,
      c(
        "Temperature_z",
        "Humidity_z",
        "Altitude_z"
      )
    ]
  )


  space_mm <- as.matrix(
    meta[
      ,
      c(
        "Longitude_z",
        "Latitude_z"
      )
    ]
  )


  vpa <- vegan::varpart(
    bray,
    host_mm,
    env_mm,
    space_mm,
    sqrt.dist = TRUE
  )


  capture.output(
    vpa,
    file = file.path(
      layer_dir,
      "VPA_full_output.txt"
    )
  )


  vpa_ind <- as.data.frame(
    vpa$part$indfract
  )


  write.csv(
    cbind(
      Fraction =
        rownames(
          vpa_ind
        ),
      vpa_ind
    ),
    file.path(
      layer_dir,
      "VPA_individual_fractions.csv"
    ),
    row.names = FALSE
  )


  # Pure fraction significance tests

  pure_host <- vegan::dbrda(

    bray ~
      Clade +
      Stage +
      Condition(
        Temperature_z,
        Humidity_z,
        Altitude_z,
        Longitude_z,
        Latitude_z
      ),

    data =
      meta,

    sqrt.dist =
      TRUE
  )


  pure_env <- vegan::dbrda(

    bray ~
      Temperature_z +
      Humidity_z +
      Altitude_z +
      Condition(
        Clade,
        Stage,
        Longitude_z,
        Latitude_z
      ),

    data =
      meta,

    sqrt.dist =
      TRUE
  )


  pure_space <- vegan::dbrda(

    bray ~
      Longitude_z +
      Latitude_z +
      Condition(
        Clade,
        Stage,
        Temperature_z,
        Humidity_z,
        Altitude_z
      ),

    data =
      meta,

    sqrt.dist =
      TRUE
  )


  set.seed(
    20260813
  )

  pure_host_test <- anova(
    pure_host,
    permutations =
      nperm
  )


  set.seed(
    20260813
  )

  pure_env_test <- anova(
    pure_env,
    permutations =
      nperm
  )


  set.seed(
    20260813
  )

  pure_space_test <- anova(
    pure_space,
    permutations =
      nperm
  )


  capture.output(
    pure_host_test,
    file = file.path(
      layer_dir,
      "VPA_pure_host_test.txt"
    )
  )


  capture.output(
    pure_env_test,
    file = file.path(
      layer_dir,
      "VPA_pure_environment_test.txt"
    )
  )


  capture.output(
    pure_space_test,
    file = file.path(
      layer_dir,
      "VPA_pure_space_test.txt"
    )
  )


  # --------------------------------------------------------
  # Extract Clade marginal PERMANOVA
  # --------------------------------------------------------

  perm_margin_df <-
    as.data.frame(
      perm_margin
    )


  perm_margin_df$Term <-
    rownames(
      perm_margin_df
    )


  clade_row <- perm_margin_df %>%

    filter(
      Term ==
        "Clade"
    )


  clade_summary <-
    if (
      nrow(
        clade_row
      ) ==
        1
    ) {

      data.frame(Layer = layer_name, Clade_R2 = clade_row$R2, Clade_P =clade_row[["Pr(>F)"]])

    } else {

      data.frame(
        Layer =
          layer_name,

        Clade_R2 =
          NA_real_,

        Clade_P =
          NA_real_
      )
    }


  # --------------------------------------------------------
  # Extract VPA pure fractions
  # --------------------------------------------------------

  adj_col <- grep(
    "Adj",
    names(
      vpa_ind
    ),
    value = TRUE
  )[1]


  get_vpa_fraction <- function(
      code
  ) {

    idx <- grep(
      code,
      rownames(
        vpa_ind
      ),
      fixed = TRUE
    )

    if (
      length(
        idx
      ) ==
        0
    ) {

      return(
        NA_real_
      )
    }

    as.numeric(
      vpa_ind[
        idx[1],
        adj_col
      ]
    )
  }


  vpa_summary <- data.frame(

    Layer =
      layer_name,

    Host =
      get_vpa_fraction(
        "[a]"
      ),

    Environment =
      get_vpa_fraction(
        "[b]"
      ),

    Space =
      get_vpa_fraction(
        "[c]"
      )
  )


  return(
    list(
      distance =
        bray,

      pcoa =
        pco_df,

      permanova_seq =
        perm_seq,

      permanova_margin =
        perm_margin,

      dbrda =
        db,

      vpa =
        vpa,

      clade_summary =
        clade_summary,

      vpa_summary =
        vpa_summary
    )
  )
}


############################################################
# 14. RUN MULTIVARIATE ANALYSIS ON ALL LAYERS
############################################################

multi_results <- list()


for (
  nm in
    names(
      layer_beta
    )
) {

  multi_results[[nm]] <- run_multivariate_layer(mat =
      layer_beta[[nm]],

    layer_name =
      nm,

    meta =
      meta
  )
}


############################################################
# PART II. CAZy LEVEL 1 DESCRIPTIVE COMPOSITION
############################################################

############################################################
# 15. CAZy LEVEL 1 STACKED COMPOSITION
############################################################

cazy1_long <- data.frame(
  SampleID =
    rownames(
      cazy1_rel
    ),
  cazy1_rel,
  check.names = FALSE
) %>%

  pivot_longer(
    cols =
      -SampleID,
    names_to =
      "CAZy_Class",
    values_to =
      "Relative_abundance"
  ) %>%

  left_join(
    meta %>%
      select(
        SampleID,
        Clade
      ),
    by =
      "SampleID"
  )


cazy1_group <- cazy1_long %>%

  group_by(
    Clade,
    CAZy_Class
  ) %>%

  summarise(
    Mean =
      mean(
        Relative_abundance
      ),
    .groups =
      "drop"
  )


p_cazy1 <- ggplot(
  cazy1_group,
  aes(
    x = Clade,
    y = Mean,
    fill = CAZy_Class
  )
) +

  geom_col(
    width = 0.7,
    colour = "black",
    linewidth = 0.15
  ) +

  scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +

  labs(
    x = "Clade",
    y =
      "Mean relative abundance\nwithin CAZy repertoire",
    fill =
      "CAZy class",
    title =
      "CAZy level-1 functional composition"
  ) +

  theme_nature()


save_fig(
  p_cazy1,
  "02_CAZy_level1_composition",
  width = 7,
  height = 5
)


############################################################
# 16. CAZy TOTAL SIGNAL
############################################################

cazy_signal_df <- data.frame(
  SampleID =
    common_samples,

  CAZy_total_signal =
    rowSums(
      cazy1_original
    )
) %>%

  left_join(
    meta %>%
      select(
        SampleID,
        Clade
      ),
    by =
      "SampleID"
  )


p_cazy_signal <- ggplot(
  cazy_signal_df,
  aes(
    Clade,
    CAZy_total_signal,
    fill = Clade
  )
) +

  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +

  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.75
  ) +

  scale_fill_manual(
    values =
      nature_clade
  ) +

  labs(
    x = "Clade",
    y =
      "Total supplied CAZy level-1 signal",
    title =
      "Overall CAZyme signal"
  ) +

  theme_nature() +

  theme(
    legend.position =
      "none"
  )


save_fig(
  p_cazy_signal,
  "02_CAZy_total_signal",
  width = 5.5,
  height = 4.5
)


############################################################
# PART III. DIFFERENTIAL FUNCTIONAL ANALYSIS
############################################################

############################################################
# 17. GENERIC MaAsLin2 FUNCTION
############################################################

run_maaslin_layer <- function(
    mat_original,
    layer_name,
    prevalence_cut,
    mean_cut,
    meta
) {

  cat(
    "\n-----------------------------------------------\n"
  )

  cat(
    "MaAsLin2:",
    layer_name,
    "\n"
  )


  prev <- colMeans(
    mat_original >
      0
  )


  mean_ab <- colMeans(
    mat_original
  )


  keep <-
    prev >=
    prevalence_cut &
    mean_ab >=
    mean_cut


  mat_da <- mat_original[
    ,
    keep,
    drop = FALSE
  ]


  cat(
    "Retained features:",
    ncol(
      mat_da
    ),
    "\n"
  )


  if (
    ncol(
      mat_da
    ) ==
      0
  ) {

    warning(
      "No features retained for ",
      layer_name
    )

    return(
      NULL
    )
  }


  layer_dir <- file.path(
    outdir,
    paste0(
      "03_MaAsLin2_",
      layer_name,
      "_",
      format(
        Sys.time(),
        "%Y%m%d_%H%M%S"
      )
    )
  )


  dir.create(
    layer_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # Test write access
  writeLines(
    "test",
    file.path(
      layer_dir,
      "write_test.txt"
    )
  )


  unlink(
    file.path(
      layer_dir,
      "write_test.txt"
    )
  )


  meta_maas <- meta %>%

    select(
      SampleID,
      Clade,
      Stage,
      Temperature_z,
      Humidity_z,
      Altitude_z,
      Latitude_z,
      Longitude_z
    ) %>%

    as.data.frame()


  rownames(
    meta_maas
  ) <-
    meta_maas$SampleID


  meta_maas$SampleID <-
    NULL


  data_maas <- as.data.frame(
    mat_da
  )


  fit <- Maaslin2::Maaslin2(

    input_data =
      data_maas,

    input_metadata =
      meta_maas,

    output =
      layer_dir,

    fixed_effects =
      c(
        "Clade",
        "Stage",
        "Temperature_z",
        "Humidity_z",
        "Altitude_z",
        "Latitude_z",
        "Longitude_z"
      ),

    random_effects =
      NULL,

    normalization =
      "NONE",

    transform =
      "LOG",

    analysis_method =
      "LM",

    correction =
      "BH",

    standardize =
      FALSE,

    min_abundance =
      0,

    min_prevalence =
      0,

    max_significance =
      fdr_cut,

    reference =
      c(
        "Clade,A"
      ),

    plot_heatmap =
      FALSE,

    plot_scatter =
      FALSE
  )


  res_file <- file.path(
    layer_dir,
    "all_results.tsv"
  )


  if (
    !file.exists(
      res_file
    )
  ) {

    warning(
      "No all_results.tsv for ",
      layer_name
    )

    return(
      NULL
    )
  }


  res <- read.delim(
    res_file,
    check.names = FALSE
  )


  clade_res <- res %>%

    filter(
      metadata ==
        "Clade",
      qval <
        fdr_cut
    ) %>%

    mutate(

      Contrast =
        case_when(
          value == "B" ~
            "B vs A",
          value == "C" ~
            "C vs A",
          TRUE ~
            paste0(
              value,
              " vs A"
            )
        ),

      Direction =
        case_when(
          value == "B" &
            coef > 0 ~
            "B enriched",

          value == "B" &
            coef < 0 ~
            "A enriched",

          value == "C" &
            coef > 0 ~
            "C enriched",

          value == "C" &
            coef < 0 ~
            "A enriched",

          TRUE ~
            NA_character_
        )
    )


  write.csv(
    clade_res,
    file.path(
      outdir,
      paste0(
        "03_",
        layer_name,
        "_MaAsLin2_significant_Clade.csv"
      )
    ),
    row.names = FALSE
  )


  # Top effects
  if (
    nrow(
      clade_res
    ) >
      0
  ) {

    top <- clade_res %>%

      group_by(
        Contrast
      ) %>%

      arrange(
        qval,
        desc(
          abs(
            coef
          )
        )
      ) %>%

      slice_head(
        n =
          top_n_features
      ) %>%

      ungroup()


    top$Display <-
      factor(
        top$feature,
        levels =
          rev(
            unique(
              top$feature
            )
          )
      )


    p <- ggplot(
      top,
      aes(
        coef,
        Display,
        fill = Direction
      )
    ) +

      geom_vline(
        xintercept = 0,
        linetype = 2,
        colour = "grey60"
      ) +

      geom_col(
        width = 0.7,
        colour = "black",
        linewidth = 0.15
      ) +

      facet_wrap(
        ~ Contrast,
        scales = "free_y"
      ) +

      scale_fill_manual(
        values = c(
          "A enriched" =
            "#E64B35",
          "B enriched" =
            "#4DBBD5",
          "C enriched" =
            "#00A087"
        )
      ) +

      labs(
        x =
          "MaAsLin2 coefficient",
        y = NULL,
        fill =
          "Higher in",
        title =
          paste0(
            layer_name,
            " differential functions"
          )
      ) +

      theme_nature(
        base_size = 9.5
      )


    save_fig(
      p,
      paste0(
        "03_",
        layer_name,
        "_MaAsLin2_top_effects"
      ),
      width = 9,
      height = 8
    )
  }


  list(
    filtered_matrix =
      mat_da,
    full_results =
      res,
    clade_results =
      clade_res,
    output_dir =
      layer_dir
  )
}


############################################################
# 18. PREPARE ORIGINAL-SCALE FUNCTIONAL MATRICES
############################################################

# KO: retain supplied KO values for modelling
# COG: converted sample-wise relative abundances
# CAZy L1: retain supplied level-1 relative values
# CAZy L2: REQUIRED relative conversion

functional_da_mats <- list(

  KO =
    ko_original,

  COG =
    cog_rel,

  CAZy_L1 =
    cazy1_original,

  CAZy_L2 =
    cazy2_rel
)


maaslin_results <- list()


for (
  nm in
    names(
      functional_da_mats
    )
) {

  maaslin_results[[nm]] <- run_maaslin_layer(mat_original =
      functional_da_mats[[nm]],

    layer_name =nm,
    prevalence_cut =
      prevalence_da[[nm]],

    mean_cut =
      mean_abundance_da[[nm]],

    meta =
      meta
  )
}


############################################################
# PART IV. KEGG PATHWAY ANALYSIS
############################################################

############################################################
# 19. GET / LOAD KEGG KO -> PATHWAY MAPPING
############################################################

kegg_mapping_file <- file.path(
  outdir,
  "04_KEGG_KO_pathway_mapping.csv"
)


if (
  run_kegg_mapping &&
  !file.exists(
    kegg_mapping_file
  )
) {

  message(
    "\nDownloading KEGG KO-pathway mapping..."
  )


  links <- tryCatch(

    KEGGREST::keggLink(
      "pathway",
      "ko"
    ),

    error =
      function(e) {

        warning(
          "KEGG mapping failed: ",
          e$message
        )

        NULL
      }
  )


  if (
    !is.null(
      links
    )
  ) {

    mapping <- data.frame(

      KO =
        sub(
          "^ko:",
          "",
          names(
            links
          )
        ),

      Pathway =
        sub(
          "^path:",
          "",
          unname(
            links
          )
        ),

      stringsAsFactors =
        FALSE
    ) %>%

      filter(
        grepl(
          "^ko[0-9]{5}$",
          Pathway
        )
      ) %>%

      distinct()


    pnames_raw <- tryCatch(

      KEGGREST::keggList(
        "pathway",
        "ko"
      ),

      error =
        function(e) {
          NULL
        }
    )


    if (
      !is.null(
        pnames_raw
      )
    ) {

      pnames <- data.frame(

        Pathway =
          sub(
            "^path:",
            "",
            names(
              pnames_raw
            )
          ),

        Pathway_name =
          unname(
            pnames_raw
          ),

        stringsAsFactors =
          FALSE
      )


      mapping <- mapping %>%

        left_join(
          pnames,
          by =
            "Pathway"
        )
    }


    write.csv(
      mapping,
      kegg_mapping_file,
      row.names = FALSE
    )
  }
}


############################################################
# 20. PATHWAY ABUNDANCE + PATHWAY ENRICHMENT
############################################################

if (
  file.exists(
    kegg_mapping_file
  )
) {
  
  ##########################################################
  # 20.1 Read KO -> KEGG pathway mapping
  ##########################################################
  
  mapping <- read.csv(
    kegg_mapping_file,
    stringsAsFactors = FALSE
  )
  
  
  ##########################################################
  # 20.2 Keep only KOs present in current dataset
  ##########################################################
  
  mapping <- mapping %>%
    
    filter(
      KO %in%
        colnames(
          ko_original
        )
    ) %>%
    
    # VERY IMPORTANT:
    # remove exact duplicate KO-pathway mappings
    distinct(
      KO,
      Pathway,
      .keep_all = TRUE
    )
  
  
  ##########################################################
  # 20.3 QC mapping
  ##########################################################
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "KEGG pathway mapping QC\n"
  )
  
  cat(
    "============================================\n"
  )
  
  cat(
    "Mapped KOs:",
    n_distinct(
      mapping$KO
    ),
    "\n"
  )
  
  cat(
    "Mapped pathways:",
    n_distinct(
      mapping$Pathway
    ),
    "\n"
  )
  
  cat(
    "KO-pathway links:",
    nrow(
      mapping
    ),
    "\n"
  )
  
  
  ##########################################################
  # Check exact duplicate mappings
  ##########################################################
  
  duplicate_mapping <- mapping %>%
    
    count(
      KO,
      Pathway
    ) %>%
    
    filter(
      n > 1
    )
  
  
  if (
    nrow(
      duplicate_mapping
    ) > 0
  ) {
    
    warning(
      "Duplicated KO-pathway mappings remain."
    )
    
  } else {
    
    cat(
      "No duplicated KO-pathway pairs detected.\n"
    )
  }
  
  
  ##########################################################
  # 20.4 How many pathways can each KO belong to?
  ##########################################################
  
  ko_pathway_degree <- mapping %>%
    
    count(
      KO,
      name =
        "N_pathways"
    ) %>%
    
    arrange(
      desc(
        N_pathways
      )
    )
  
  
  write.csv(
    ko_pathway_degree,
    file.path(
      outdir,
      "04_KEGG_KO_number_of_pathways.csv"
    ),
    row.names = FALSE
  )
  
  
  cat(
    "\nNumber of pathways per KO:\n"
  )
  
  print(
    summary(
      ko_pathway_degree$N_pathways
    )
  )
  
  
  ##########################################################
  # 20.5 Convert KO table to long format
  ##########################################################
  
  ko_long <- data.frame(
    
    SampleID =
      rownames(
        ko_original
      ),
    
    ko_original,
    
    check.names =
      FALSE
  ) %>%
    
    pivot_longer(
      
      cols =
        -SampleID,
      
      names_to =
        "KO",
      
      values_to =
        "Abundance"
    )
  
  
  ##########################################################
  # 20.6 Join KO abundance with KEGG pathway annotations
  #
  # many-to-many is EXPECTED:
  #
  # one KO occurs in many samples
  # AND
  # one KO can map to multiple pathways
  ##########################################################
  
  path_long <- ko_long %>%
    
    inner_join(
      
      mapping,
      
      by =
        "KO",
      
      relationship =
        "many-to-many"
    )
  
  
  ##########################################################
  # 20.7 QC joined table
  ##########################################################
  
  cat(
    "\nKO abundance rows:",
    nrow(
      ko_long
    ),
    "\n"
  )
  
  cat(
    "KO-pathway abundance rows after join:",
    nrow(
      path_long
    ),
    "\n"
  )
  
  
  ##########################################################
  # 20.8 Pathway abundance scores
  ##########################################################
  
  path_scores <- path_long %>%
    
    group_by(
      
      SampleID,
      
      Pathway,
      
      Pathway_name
    ) %>%
    
    summarise(
      
      Pathway_sum =
        sum(
          Abundance,
          na.rm = TRUE
        ),
      
      N_KOs =
        n_distinct(
          KO
        ),
      
      .groups =
        "drop"
    )
  
  
  ##########################################################
  # 20.9 Export pathway abundance
  ##########################################################
  
  write.csv(
    
    path_scores,
    
    file.path(
      outdir,
      "04_KEGG_pathway_scores.csv"
    ),
    
    row.names =
      FALSE
  )
  
  
  ##########################################################
  # 20.10 Optional:
  # Create pathway × sample matrix
  ##########################################################
  
  pathway_matrix <- path_scores %>%
    
    select(
      SampleID,
      Pathway,
      Pathway_sum
    ) %>%
    
    pivot_wider(
      
      names_from =
        Pathway,
      
      values_from =
        Pathway_sum,
      
      values_fill =
        0
    )
  
  
  write.csv(
    
    pathway_matrix,
    
    file.path(
      outdir,
      "04_KEGG_pathway_abundance_matrix.csv"
    ),
    
    row.names =
      FALSE
  )
  
  
  ##########################################################
  # 20.11 Enrichment of differential KOs
  ##########################################################
  
  ko_da_res <-
    maaslin_results$KO
  
  
  if (
    !is.null(
      ko_da_res
    ) &&
    nrow(
      ko_da_res$clade_results
    ) >
    0
  ) {
    
    ########################################################
    # KOs actually tested by MaAsLin2
    ########################################################
    
    tested_kos <-
      colnames(
        ko_da_res$filtered_matrix
      )
    
    
    ########################################################
    # Mapping background
    ########################################################
    
    mapping_tested <- mapping %>%
      
      filter(
        KO %in%
          tested_kos
      ) %>%
      
      distinct(
        KO,
        Pathway,
        .keep_all = TRUE
      )
    
    
    ########################################################
    # Remove broad overview pathways if requested
    ########################################################
    
    if (
      remove_broad_kegg_pathways
    ) {
      
      mapping_tested <- mapping_tested %>%
        
        filter(
          !Pathway %in%
            broad_kegg_pathways
        )
    }
    
    
    ########################################################
    # Statistical background
    #
    # IMPORTANT:
    # only mapped and tested KOs
    ########################################################
    
    background_kos <- unique(
      mapping_tested$KO
    )
    
    
    cat(
      "\nKEGG enrichment background KOs:",
      length(
        background_kos
      ),
      "\n"
    )
    
    
    ########################################################
    # Significant KO sets
    ########################################################
    
    sig_sets <- ko_da_res$clade_results %>%
      
      transmute(
        
        KO =
          feature,
        
        Contrast,
        
        Direction
      ) %>%
      
      filter(
        KO %in%
          background_kos
      )
    
    
    ########################################################
    # Enrichment loop
    ########################################################
    
    enrichment_list <- list()
    
    counter <- 1
    
    
    for (
      contrast in
      unique(
        sig_sets$Contrast
      )
    ) {
      
      dirs <- unique(
        sig_sets$Direction[
          sig_sets$Contrast ==
            contrast
        ]
      )
      
      
      for (
        direction in
        dirs
      ) {
        
        sig_kos <- sig_sets %>%
          
          filter(
            
            Contrast ==
              contrast,
            
            Direction ==
              direction
          ) %>%
          
          pull(
            KO
          ) %>%
          
          unique()
        
        
        cat(
          "\n",
          contrast,
          " / ",
          direction,
          ": ",
          length(
            sig_kos
          ),
          " significant mapped KOs\n",
          sep = ""
        )
        
        
        if (
          length(
            sig_kos
          ) <
          2
        ) {
          
          next
        }
        
        
        for (
          pw in
          unique(
            mapping_tested$Pathway
          )
        ) {
          
          ##################################################
          # KOs belonging to this pathway
          ##################################################
          
          pw_kos <- mapping_tested %>%
            
            filter(
              Pathway ==
                pw
            ) %>%
            
            pull(
              KO
            ) %>%
            
            unique()
          
          
          ##################################################
          # Fisher table
          ##################################################
          
          a <- length(
            intersect(
              sig_kos,
              pw_kos
            )
          )
          
          b <-
            length(
              sig_kos
            ) -
            a
          
          c <- length(
            setdiff(
              pw_kos,
              sig_kos
            )
          )
          
          d <-
            length(
              background_kos
            ) -
            a -
            b -
            c
          
          
          if (
            min(
              a,
              b,
              c,
              d
            ) <
            0
          ) {
            
            next
          }
          
          
          ##################################################
          # Skip pathways with no significant KO
          #
          # This greatly reduces unnecessary Fisher tests.
          ##################################################
          
          if (
            a ==
            0
          ) {
            
            next
          }
          
          
          ##################################################
          # Fisher over-representation
          ##################################################
          
          ft <- fisher.test(
            
            matrix(
              
              c(
                a,
                b,
                c,
                d
              ),
              
              nrow =
                2,
              
              byrow =
                TRUE
            ),
            
            alternative =
              "greater"
          )
          
          
          ##################################################
          # Save result
          ##################################################
          
          enrichment_list[[counter]] <- data.frame(
            
            Contrast =
              contrast,
            
            Direction =
              direction,
            
            Pathway =
              pw,
            
            Significant_KOs_in_pathway =
              a,
            
            Significant_KOs_total =
              length(
                sig_kos
              ),
            
            Background_KOs_in_pathway =
              length(
                pw_kos
              ),
            
            Background_KOs_total =
              length(
                background_kos
              ),
            
            OddsRatio =
              unname(
                ft$estimate
              ),
            
            P_value =
              ft$p.value,
            
            stringsAsFactors =
              FALSE
          )
          
          
          counter <-
            counter +
            1
        }
      }
    }
    
    
    ########################################################
    # 20.12 Combine enrichment results
    ########################################################
    
    if (
      length(
        enrichment_list
      ) >
      0
    ) {
      
      enrich <- bind_rows(
        enrichment_list
      ) %>%
        
        group_by(
          Contrast,
          Direction
        ) %>%
        
        mutate(
          
          Q_value =
            p.adjust(
              P_value,
              method =
                "BH"
            )
        ) %>%
        
        ungroup() %>%
        
        left_join(
          
          mapping %>%
            
            select(
              Pathway,
              Pathway_name
            ) %>%
            
            distinct(
              Pathway,
              .keep_all = TRUE
            ),
          
          by =
            "Pathway"
        ) %>%
        
        mutate(
          
          FoldEnrichment =
            
            (
              Significant_KOs_in_pathway /
                Significant_KOs_total
            ) /
            (
              Background_KOs_in_pathway /
                Background_KOs_total
            )
        )
      
      
      ######################################################
      # Export all enrichment results
      ######################################################
      
      write.csv(
        
        enrich,
        
        file.path(
          outdir,
          "04_KEGG_differential_KO_pathway_enrichment.csv"
        ),
        
        row.names =
          FALSE
      )
      
      
      ######################################################
      # 20.13 Significant enrichment
      ######################################################
      
      enrich_sig <- enrich %>%
        
        filter(
          
          Q_value <
            fdr_cut,
          
          Significant_KOs_in_pathway >=
            2
        )
      
      
      write.csv(
        
        enrich_sig,
        
        file.path(
          outdir,
          "04_KEGG_significant_pathway_enrichment.csv"
        ),
        
        row.names =
          FALSE
      )
      
      
      ######################################################
      # 20.14 Plot
      ######################################################
      
      if (
        nrow(
          enrich_sig
        ) >
        0
      ) {
        
        enrich_plot <- enrich_sig %>%
          
          group_by(
            Contrast,
            Direction
          ) %>%
          
          arrange(
            Q_value,
            desc(
              FoldEnrichment
            )
          ) %>%
          
          slice_head(
            n =
              15
          ) %>%
          
          ungroup() %>%
          
          mutate(
            
            Display =
              ifelse(
                is.na(
                  Pathway_name
                ),
                Pathway,
                Pathway_name
              ),
            
            minus_log_q =
              -log10(
                pmax(
                  Q_value,
                  1e-300
                )
              )
          )
        
        
        p_enrich <- ggplot(
          
          enrich_plot,
          
          aes(
            x =
              FoldEnrichment,
            
            y =
              reorder(
                Display,
                FoldEnrichment
              )
          )
        ) +
          
          geom_point(
            
            aes(
              
              size =
                Significant_KOs_in_pathway,
              
              fill =
                minus_log_q
            ),
            
            shape =
              21,
            
            colour =
              "black",
            
            stroke =
              0.25
          ) +
          
          facet_grid(
            
            Direction ~
              Contrast,
            
            scales =
              "free_y",
            
            space =
              "free_y"
          ) +
          
          scale_fill_gradient(
            
            low =
              "#4DBBD5",
            
            high =
              "#E64B35"
          ) +
          
          scale_size_continuous(
            range =
              c(
                2.5,
                8
              )
          ) +
          
          labs(
            
            x =
              "Fold enrichment",
            
            y =
              NULL,
            
            size =
              "Differential KOs",
            
            fill =
              "-log10(q)",
            
            title =
              "KEGG pathway enrichment of differential KOs"
          ) +
          
          theme_nature(
            base_size =
              9
          )
        
        
        save_fig(
          
          p_enrich,
          
          "04_KEGG_pathway_enrichment",
          
          width =
            11,
          
          height =
            9
        )
      }
      
    } else {
      
      message(
        "No KEGG pathway enrichment results were generated."
      )
    }
  }
}


############################################################
# PART V. TAXONOMY-FUNCTION CONCORDANCE
############################################################

############################################################
# 21. PAIRWISE MANTEL AMONG ALL LAYERS
############################################################

distance_list <- lapply(
  multi_results,
  function(x) {
    x$distance
  }
)


mantel_pairs <- combn(
  names(
    distance_list
  ),
  2,
  simplify = FALSE
)


mantel_results <- purrr::map_dfr(
  mantel_pairs,
  function(pair) {

    set.seed(
      20260813
    )


    m <- vegan::mantel(

      distance_list[[pair[1]]],

      distance_list[[pair[2]]],

      method =
        "spearman",

      permutations =
        nperm
    )


    data.frame(

      Layer1 =
        pair[1],

      Layer2 =
        pair[2],

      Mantel_r =
        as.numeric(
          m$statistic
        ),

      P_value =
        m$signif
    )
  }
)


mantel_results <- mantel_results %>%

  mutate(
    Q_value =
      p.adjust(
        P_value,
        method = "BH"
      )
  )


write.csv(
  mantel_results,
  file.path(
    outdir,
    "05_cross_layer_Mantel_results.csv"
  ),
  row.names = FALSE
)


############################################################
# 22. MANTEL MATRIX
############################################################

all_layers <- names(
  distance_list
)


mantel_mat <- matrix(
  1,
  nrow =
    length(
      all_layers
    ),
  ncol =
    length(
      all_layers
    ),
  dimnames =
    list(
      all_layers,
      all_layers
    )
)


for (
  i in
    seq_len(
      nrow(
        mantel_results
      )
    )
) {

  a <-
    mantel_results$Layer1[i]

  b <-
    mantel_results$Layer2[i]

  mantel_mat[
    a,
    b
  ] <-
    mantel_results$Mantel_r[i]

  mantel_mat[
    b,
    a
  ] <-
    mantel_results$Mantel_r[i]
}


write.csv(
  mantel_mat,
  file.path(
    outdir,
    "05_cross_layer_Mantel_matrix.csv"
  )
)


pheatmap::pheatmap(

  mantel_mat,

  cluster_rows =
    FALSE,

  cluster_cols =
    FALSE,

  display_numbers =
    TRUE,

  number_format =
    "%.2f",

  border_color =
    "white",

  main =
    "Taxonomic-functional concordance (Mantel r)",

  filename =
    file.path(
      outdir,
      "05_cross_layer_Mantel_heatmap.pdf"
    ),

  width = 7,
  height = 6
)


############################################################
# 23. PROCRUSTES / PROTEST BETWEEN GENUS AND EACH FUNCTION
############################################################

procrustes_results <- list()

proc_counter <- 1


for (
  functional_layer in
    c(
      "KO",
      "COG",
      "CAZy_L1",
      "CAZy_L2"
    )
) {

  genus_dist <-
    distance_list$Genus


  func_dist <-
    distance_list[[functional_layer]]


  # Use enough PCoA axes to capture multivariate structure,
  # but remain below sample count.
  k_axes <- min(
    10,
    nrow(
      meta
    ) - 2
  )


  genus_pco <- cmdscale(
    genus_dist,
    k =
      k_axes,
    eig = TRUE,
    add = TRUE
  )$points


  func_pco <- cmdscale(
    func_dist,
    k =
      k_axes,
    eig = TRUE,
    add = TRUE
  )$points


  proc <- vegan::procrustes(
    genus_pco,
    func_pco,
    symmetric = TRUE
  )


  set.seed(
    20260813
  )


  protest <- vegan::protest(
    genus_pco,
    func_pco,
    permutations =
      nperm
  )


  procrustes_results[[proc_counter]] <- data.frame(

    Taxonomy =
      "Genus",

    Function =
      functional_layer,

    Procrustes_correlation =
      sqrt(
        1 -
          proc$ss
      ),

    M2 =
      proc$ss,

    P_value =
      protest$signif
  )


  proc_counter <-
    proc_counter + 1
}


procrustes_df <- bind_rows(
  procrustes_results
) %>%

  mutate(
    Q_value =
      p.adjust(
        P_value,
        method = "BH"
      )
  )


write.csv(
  procrustes_df,
  file.path(
    outdir,
    "05_taxonomy_function_PROTEST_results.csv"
  ),
  row.names = FALSE
)


p_proc <- ggplot(
  procrustes_df,
  aes(
    Function,
    Procrustes_correlation,
    fill = Function
  )
) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label =
        paste0(
          "q=",
          signif(
            Q_value,
            2
          )
        )
    ),
    vjust = -0.35,
    size = 3.5
  ) +

  scale_fill_manual(
    values =
      layer_colours[
        names(
          layer_colours
        ) !=
          "Genus"
      ]
  ) +

  labs(
    x = NULL,
    y =
      "Procrustes correlation",
    title =
      "Taxonomy-function concordance"
  ) +

  theme_nature() +

  theme(
    legend.position =
      "none"
  )


save_fig(
  p_proc,
  "05_taxonomy_function_Procrustes_summary",
  width = 6.5,
  height = 4.8
)


############################################################
# PART VI. CROSS-LAYER SYNTHESIS
############################################################

############################################################
# 24. MARGINAL CLADE R2 COMPARISON
############################################################

clade_r2 <- bind_rows(
  lapply(
    multi_results,
    function(x) {
      x$clade_summary
    }
  )
)


write.csv(
  clade_r2,
  file.path(
    outdir,
    "06_cross_layer_Clade_marginal_R2.csv"
  ),
  row.names = FALSE
)


p_r2 <- ggplot(
  clade_r2,
  aes(
    Layer,
    Clade_R2,
    fill = Layer
  )
) +

  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label =
        paste0(
          "P=",
          signif(
            Clade_P,
            2
          )
        )
    ),
    vjust = -0.35,
    size = 3.4
  ) +

  scale_fill_manual(
    values =
      layer_colours
  ) +

  scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand =
      expansion(
        mult =
          c(
            0,
            0.15
          )
      )
  ) +

  labs(
    x = NULL,
    y =
      "Marginal PERMANOVA R² for Clade",
    title =
      "Independent Clade effect across taxonomic and functional layers"
  ) +

  theme_nature() +

  theme(
    legend.position =
      "none"
  )


save_fig(
  p_r2,
  "06_cross_layer_Clade_R2",
  width = 7,
  height = 5
)


############################################################
# 25. VPA PURE-FRACTION COMPARISON
############################################################

vpa_cross <- bind_rows(
  lapply(
    multi_results,
    function(x) {
      x$vpa_summary
    }
  )
) %>%

  pivot_longer(
    cols =
      c(
        Host,
        Environment,
        Space
      ),
    names_to =
      "Fraction",
    values_to =
      "Adjusted_R2"
  )


write.csv(
  vpa_cross,
  file.path(
    outdir,
    "06_cross_layer_VPA_pure_fractions.csv"
  ),
  row.names = FALSE
)


p_vpa_cross <- ggplot(
  vpa_cross,
  aes(
    Layer,
    Adjusted_R2,
    fill = Fraction
  )
) +

  geom_col(
    position =
      position_dodge(
        width = 0.75
      ),
    width = 0.68,
    colour = "black",
    linewidth = 0.18
  ) +

  scale_fill_manual(
    values =
      vpa_colours
  ) +

  labs(
    x = NULL,
    y =
      "Pure adjusted R²",
    fill =
      "Explanatory set",
    title =
      "Host, environmental and spatial contributions across layers"
  ) +

  theme_nature()


save_fig(
  p_vpa_cross,
  "06_cross_layer_VPA_pure_fractions",
  width = 8,
  height = 5.2
)


############################################################
# 26. DIFFERENTIAL FEATURE COUNTS ACROSS FUNCTIONAL LAYERS
############################################################

da_counts <- purrr::imap_dfr(
  maaslin_results,
  function(x, nm) {

    if (
      is.null(
        x
      )
    ) {

      return(
        data.frame(
          Layer =
            nm,
          Contrast =
            NA_character_,
          Significant_features =
            0
        )
      )
    }


    if (
      nrow(
        x$clade_results
      ) ==
        0
    ) {

      return(
        data.frame(
          Layer =
            nm,
          Contrast =
            NA_character_,
          Significant_features =
            0
        )
      )
    }


    x$clade_results %>%

      count(
        Contrast,
        name =
          "Significant_features"
      ) %>%

      mutate(
        Layer =
          nm
      ) %>%

      select(
        Layer,
        Contrast,
        Significant_features
      )
  }
)


write.csv(
  da_counts,
  file.path(
    outdir,
    "06_functional_DA_feature_counts.csv"
  ),
  row.names = FALSE
)


############################################################
# PART VII. REPRODUCIBILITY
############################################################

############################################################
# 27. SESSION INFO
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)


############################################################
# 28. ANALYSIS NOTES
############################################################

notes <- c(

  paste0(
    "Common samples analysed: ",
    length(
      common_samples
    )
  ),

  paste0(
    "COG feature mode: ",
    cog_feature_mode
  ),

  paste0(
    "COG Others removed: ",
    drop_cog_others
  ),

  "CAZy level 2 absolute abundance was transposed to sample x family and converted to sample-wise relative abundance.",

  "CAZy level 1 supplied values were retained for differential modelling; a within-CAZy TSS composition was used for Bray-Curtis/dbRDA.",

  "COG uploaded data were treated as counts/abundances rather than pre-normalized relative abundances and converted to within-COG sample-wise relative abundance.",

  "Genus, KO, COG and CAZy compositional Bray-Curtis distances use sample-wise unit-sum profiles.",

  "Genus is included primarily as the taxonomic reference layer for cross-layer concordance.",

  "Functional MaAsLin2 models adjust for Clade, Stage, Temperature, Humidity, Altitude, Latitude and Longitude.",

  "All multivariate layers use the same sequential/marginal PERMANOVA, PERMDISP, dbRDA and VPA framework.",

  "Cross-layer Mantel and Procrustes analyses test concordance, not causal relationships.",

  "KEGG pathway scores overlap because a KO can map to multiple pathways.",

  "KEGG over-representation uses KOs retained for MaAsLin2 as the background."
)


writeLines(
  notes,
  file.path(
    outdir,
    "00_analysis_notes.txt"
  )
)


############################################################
# 29. FINAL CONSOLE SUMMARY
############################################################

cat("\n============================================================\n")

cat("INTEGRATED FUNCTIONAL MICROBIOME ANALYSIS COMPLETED\n")

cat("============================================================\n")

cat(  "Common samples:",
  length(
    common_samples
  ),
  "\n"
)

print(
  qc_layers
)

cat(
  "\nMain output directory:\n",
  normalizePath(
    outdir
  ),
  "\n"
)

cat("\nPriority result files:\n")

cat(
  "  00_integrated_layer_QC.csv\n",
  "  00_CAZy_level2_converted_relative_abundance.csv\n",
  "  00_COG_converted_relative_abundance.csv\n",
  "  01_multivariate_*/PERMANOVA_marginal.csv\n",
  "  01_multivariate_*/PERMDISP_Clade.txt\n",
  "  01_multivariate_*/dbRDA_marginal.csv\n",
  "  01_multivariate_*/dbRDA_VIF.csv\n",
  "  01_multivariate_*/VPA_individual_fractions.csv\n",
  "  03_*_MaAsLin2_significant_Clade.csv\n",
  "  04_KEGG_differential_KO_pathway_enrichment.csv\n",
  "  05_cross_layer_Mantel_results.csv\n",
  "  05_taxonomy_function_PROTEST_results.csv\n",
  "  06_cross_layer_Clade_marginal_R2.csv\n",
  "  06_cross_layer_VPA_pure_fractions.csv\n"
)

cat("============================================================\n")
