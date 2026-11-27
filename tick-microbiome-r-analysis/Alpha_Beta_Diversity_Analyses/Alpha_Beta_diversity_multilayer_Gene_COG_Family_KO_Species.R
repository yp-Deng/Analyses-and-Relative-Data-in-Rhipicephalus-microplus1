############################################################
# MULTI-LAYER ALPHA + BETA DIVERSITY ANALYSIS
#
# Layers:
#   Gene
#   COG
#   Family
#   KO
#   Species
#
# Main improvements:
#   1. ALL five data layers are analysed automatically.
#   2. Clade and Stage are fitted in the SAME statistical model.
#   3. If every Clade × Stage cell has sufficient replication,
#      the model automatically includes Clade × Stage interaction.
#   4. Alpha:
#        robust LM + HC3 covariance
#        Type II / III tests
#        EMM pairwise comparisons + BH correction
#   5. Beta:
#        Bray-Curtis
#        Lingoes-corrected PCoA
#        sequential + marginal PERMANOVA
#        PERMDISP for Clade and Stage
#        optional pairwise PERMANOVA + BH
#   6. Sample matching is always by SampleID.
#   7. Final output:
#        Alpha_diversity_multilayer_panel.pdf/png
#        Beta_diversity_multilayer_panel.pdf/png
#
############################################################

rm(list = ls())
set.seed(20260813)

############################################################
# 0. USER SETTINGS
############################################################

outdir <- "Alpha_Beta_multilayer_results"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 0A. METADATA
############################################################

metadata_file <- "metadata.xlsx"

metadata_candidates <- c(
  "metadata.xlsx",
  "meta_clade_stage.txt",
  "meta.txt"
)
############################################################
# 0B. ALPHA FILES
############################################################
# These names follow your original analysis files.

alpha_files <- c(
  Gene    = "gene_alpha.txt",
  COG     = "COG_alpha.txt",
  Family  = "family_alpha.txt",
  KO      = "KO_alpha.txt",
  Species = "species_alpha.txt"
)


# Main alpha metric.
#
# Your original code used:
# m = "shannon_e"
#
# Add other existing metrics if present:
# c("shannon_e", "richness", "chao1")
alpha_metrics <- c(
  "shannon_e",
  "richness",
  "chao1"
)


############################################################
# 0C. BETA FILES
############################################################
# Change these names ONLY if your actual files have different names.
#
# The reader automatically detects:
#   samples × features
# OR
#   features × samples
#
# based on SampleID matching with metadata.

beta_files <- c(
  Gene    = "Gene.xlsx",
  COG     = "COG.xlsx",
  Family  = "Family.xlsx",
  KO      = "KO.xlsx",
  Species = "Species.xlsx"
)


############################################################
# 0D. MODEL SETTINGS
############################################################

allow_interaction <- TRUE

min_cell_n_for_interaction <- 3

nperm <- 9999

p_adjust_method <- "BH"


############################################################
# 0E. BETA FILTERING
############################################################
#
# Set to 0 if you do NOT want prevalence filtering.
#
# 0.10 = a feature must be present in >=10% of samples.

beta_prevalence_cut <- c(
  Gene    = 0.25,
  COG     = 0.10,
  Family  = 0.10,
  KO      = 0.10,
  Species = 0.10
)


############################################################
# 1. PACKAGES
############################################################

packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "vegan",
  "ape",
  "car",
  "emmeans",
  "sandwich",
  "scales"
)


to_install <- packages[
  !vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (
  length(
    to_install
  ) > 0
) {

  install.packages(
    to_install
  )
}


suppressPackageStartupMessages({

  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(ape)
  library(car)
  library(emmeans)
  library(sandwich)
  library(scales)

})


############################################################
# 2. COLOUR / SHAPE SETTINGS
############################################################

nature_clade <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)


layer_colours <- c(
  "Gene"    = "#3C5488",
  "COG"     = "#00A087",
  "Family"  = "#E64B35",
  "KO"      = "#4DBBD5",
  "Species" = "#F39B7F"
)


stage_shape_pool <- c(
  21,
  22,
  24,
  23,
  25,
  8
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

      strip.background =
        element_blank(),

      strip.text =
        element_text(
          face = "bold",
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
        )
    )
}


save_plot <- function(
    plot,
    name,
    width,
    height
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
# 3. HELPER FUNCTIONS
############################################################

resolve_metadata <- function() {

  hit <- metadata_candidates[
    file.exists(
      metadata_candidates
    )
  ]


  if (
    length(
      hit
    ) ==
      0
  ) {

    stop(
      "No metadata file found.\n",
      "Expected one of:\n",
      paste(
        metadata_candidates,
        collapse = "\n"
      )
    )
  }


  hit[1]
}


read_metadata_flexible <- function(
    file
) {

  ext <- tolower(
    tools::file_ext(
      file
    )
  )


  if (
    ext %in%
      c(
        "xlsx",
        "xls"
      )
  ) {

    x <- readxl::read_excel(
      file
    )

  } else {

    x <- read.delim(
      file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }


  x <- as.data.frame(
    x,
    check.names = FALSE
  )


  # If SampleID is not explicitly named,
  # use the first column.
  if (
    !"SampleID" %in%
      names(
        x
      )
  ) {

    names(
      x
    )[1] <- "SampleID"
  }


  # Compatibility with older metadata that used "group".
  if (
    !"Clade" %in%
      names(
        x
      ) &&
    "group" %in%
      names(
        x
      )
  ) {

    names(
      x
    )[
      names(
        x
      ) ==
        "group"
    ] <- "Clade"
  }


  x
}


read_alpha_file <- function(
    file,
    layer
) {

  if (
    !file.exists(
      file
    )
  ) {

    warning(
      "Alpha file missing and skipped: ",
      file
    )

    return(
      NULL
    )
  }


  x <- read.table(

    file,

    sep = "\t",

    header = TRUE,

    stringsAsFactors = FALSE,

    row.names = 1,

    check.names = FALSE,

    comment.char = ""
  )


  x$SampleID <-
    rownames(
      x
    )


  x$Layer <-
    layer


  x
}


read_beta_file <- function(
    file,
    sample_ids,
    layer
) {

  if (
    !file.exists(
      file
    )
  ) {

    warning(
      "Beta file missing and skipped: ",
      file
    )

    return(
      NULL
    )
  }


  ext <- tolower(
    tools::file_ext(
      file
    )
  )


  # First try tab-delimited format, including old *.xls text files.
  x <- tryCatch(

    read.table(

      file,

      header = TRUE,

      row.names = 1,

      sep = "\t",

      comment.char = "",

      check.names = FALSE,

      stringsAsFactors = FALSE
    ),

    error =
      function(e) {

        if (
          ext %in%
            c(
              "xlsx",
              "xls"
            )
        ) {

          y <- readxl::read_excel(
            file
          )


          y <- as.data.frame(
            y,
            check.names = FALSE
          )


          rownames(
            y
          ) <- as.character(
            y[[1]]
          )


          y[
            ,
            -1,
            drop = FALSE
          ]

        } else {

          stop(e)
        }
      }
  )


  x <- as.data.frame(
    x,
    check.names = FALSE
  )


  ##########################################################
  # DETECT ORIENTATION USING METADATA SAMPLE IDS
  ##########################################################

  row_match <- sum(
    rownames(
      x
    ) %in%
      sample_ids
  )


  col_match <- sum(
    colnames(
      x
    ) %in%
      sample_ids
  )


  cat(
    "\nBeta orientation check:",
    layer,
    "\n"
  )

  cat(
    "Sample IDs matched in rows:",
    row_match,
    "\n"
  )

  cat(
    "Sample IDs matched in columns:",
    col_match,
    "\n"
  )


  if (
    row_match >
      col_match
  ) {

    # samples are rows
    mat <- x

  } else if (
    col_match >
      row_match
  ) {

    # samples are columns
    mat <- t(
      as.matrix(
        x
      )
    )

    mat <- as.data.frame(
      mat,
      check.names = FALSE
    )

  } else {

    stop(
      "Cannot determine orientation of beta file: ",
      file,
      "\nSampleID matches are equal in rows and columns."
    )
  }


  ##########################################################
  # FORCE NUMERIC
  ##########################################################

  mat <- mat %>%

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


  if (
    anyNA(
      mat
    )
  ) {

    stop(
      "NA/non-numeric values detected after reading ",
      file
    )
  }


  if (
    any(
      mat <
        0
    )
  ) {

    stop(
      "Negative abundance detected in ",
      file
    )
  }


  mat
}

read_gene_beta_fast <- function(
    file,
    sample_ids
) {
  
  cat("\nReading large Gene beta table...\n")
  
  # 直接用 fread，速度通常明显快于 read.table
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }
  
  gene_raw <- data.table::fread(
    file,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = TRUE
  )
  
  # 第一列作为 feature ID
  feature_id <- gene_raw[[1]]
  
  gene_raw[[1]] <- NULL
  
  # 只保留 metadata 中真正存在的样本列
  keep_samples <- intersect(
    colnames(gene_raw),
    sample_ids
  )
  
  if (length(keep_samples) < 3) {
    stop(
      "Too few Gene sample columns matched metadata."
    )
  }
  
  cat(
    "Matched Gene samples:",
    length(keep_samples),
    "\n"
  )
  
  # 只保留这69个样本列
  gene_raw <- gene_raw[
    ,
    keep_samples,
    drop = FALSE
  ]
  
  # fread通常已经读成numeric
  gene_mat <- as.matrix(gene_raw)
  
  storage.mode(gene_mat) <- "numeric"
  
  rownames(gene_mat) <- feature_id
  
  # 当前是:
  # features × samples
  #
  # 转成:
  # samples × features
  gene_mat <- t(gene_mat)
  
  # 按metadata顺序排列
  gene_mat <- gene_mat[
    sample_ids[
      sample_ids %in%
        rownames(gene_mat)
    ],
    ,
    drop = FALSE
  ]
  
  # 删除全0 feature
  keep_nonzero <- colSums(
    gene_mat
  ) > 0
  
  gene_mat <- gene_mat[
    ,
    keep_nonzero,
    drop = FALSE
  ]
  
  cat(
    "Gene matrix after removing zero features:",
    nrow(gene_mat),
    "samples ×",
    ncol(gene_mat),
    "features\n"
  )
  
  return(gene_mat)
}

format_p <- function(
    p
) {

  if (
    is.na(
      p
    )
  ) {

    return(
      "NA"
    )
  }


  if (
    p <
      0.0001
  ) {

    return(
      "<0.0001"
    )
  }


  sprintf(
    "%.4f",
    p
  )
}


############################################################
# 4. READ METADATA
############################################################

metadata_file <- resolve_metadata()


meta <- read_metadata_flexible(
  metadata_file
)


required_meta <- c(
  "SampleID",
  "Clade",
  "Stage"
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
    "Metadata must contain: ",
    paste(
      required_meta,
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
      )
  )


if (
  anyDuplicated(
    meta$SampleID
  ) >
    0
) {

  stop(
    "Duplicated SampleID in metadata."
  )
}


############################################################
# 5. DETERMINE JOINT MODEL
############################################################

cell_table <- with(
  meta,
  table(
    Clade,
    Stage
  )
)


write.csv(
  as.data.frame(
    cell_table
  ),
  file.path(
    outdir,
    "00_Clade_Stage_sample_counts.csv"
  ),
  row.names = FALSE
)


cat(
  "\nClade × Stage sample distribution:\n"
)

print(
  cell_table
)


positive_cells <- cell_table[
  cell_table >
    0
]


use_interaction <-
  allow_interaction &&
  length(
    positive_cells
  ) ==
    length(
      cell_table
    ) &&
  min(
    positive_cells
  ) >=
    min_cell_n_for_interaction


model_text <- ifelse(
  use_interaction,
  "Clade * Stage",
  "Clade + Stage"
)


cat(
  "\nJoint model selected:",
  model_text,
  "\n"
)


############################################################
# PART I. MULTI-LAYER ALPHA DIVERSITY
############################################################

############################################################
# 6. READ ALL ALPHA FILES
############################################################

alpha_list <- list()


for (
  layer in
    names(
      alpha_files
    )
) {

  alpha_list[[layer]] <- read_alpha_file(
    alpha_files[[layer]],
    layer
  )
}


alpha_list <- alpha_list[
  !vapply(
    alpha_list,
    is.null,
    logical(
      1
    )
  )
]


if (
  length(
    alpha_list
  ) ==
    0
) {

  stop(
    "No alpha-diversity files were found."
  )
}


############################################################
# 7. COMBINE AVAILABLE ALPHA DATA
############################################################

alpha_all <- bind_rows(
  alpha_list
)


available_metrics <- intersect(
  alpha_metrics,
  names(
    alpha_all
  )
)


if (
  length(
    available_metrics
  ) ==
    0
) {

  stop(
    "None of the requested alpha metrics were found."
  )
}


alpha_long <- alpha_all %>%

  select(
    SampleID,
    Layer,
    all_of(
      available_metrics
    )
  ) %>%

  pivot_longer(
    cols =
      all_of(
        available_metrics
      ),
    names_to =
      "Metric",
    values_to =
      "Value"
  ) %>%

  inner_join(
    meta %>%
      select(
        SampleID,
        Clade,
        Stage
      ),
    by =
      "SampleID"
  )


############################################################
# 8. ALPHA JOINT MODELS FOR EVERY LAYER × METRIC
############################################################

alpha_tests <- list()
alpha_pairwise <- list()
alpha_counter <- 1
pair_counter <- 1


old_contrasts <- options(
  "contrasts"
)


options(
  contrasts = c(
    "contr.sum",
    "contr.poly"
  )
)


for (
  layer in
    unique(
      alpha_long$Layer
    )
) {

  for (
    metric in
      unique(
        alpha_long$Metric
      )
  ) {

    dat <- alpha_long %>%

      filter(
        Layer ==
          layer,
        Metric ==
          metric,
        is.finite(
          Value
        )
      )


    if (
      nrow(
        dat
      ) <
        6
    ) {

      next
    }


    if (
      use_interaction
    ) {

      form <-
        Value ~
        Clade *
        Stage

    } else {

      form <-
        Value ~
        Clade +
        Stage
    }


    fit <- lm(
      form,
      data =
        dat
    )


    vcov_hc3 <- sandwich::vcovHC(
      fit,
      type =
        "HC3"
    )


    if (
      use_interaction
    ) {

      test <- car::Anova(
        fit,
        type = 3,
        vcov. =
          vcov_hc3
      )

    } else {

      test <- car::Anova(
        fit,
        type = 2,
        vcov. =
          vcov_hc3
      )
    }


    test_df <- as.data.frame(
      test
    )


    test_df$Term <-
      rownames(
        test_df
      )


    rownames(
      test_df
    ) <- NULL


    test_df$Layer <-
      layer


    test_df$Metric <-
      metric


    alpha_tests[[alpha_counter]] <- test_df


    alpha_counter <-
      alpha_counter +
      1


    ########################################################
    # EMM: CLADE
    ########################################################

    emm_clade <- emmeans::emmeans(
      fit,
      ~ Clade,
      vcov. =
        vcov_hc3
    )


    pair_clade <- as.data.frame(
      pairs(
        emm_clade,
        adjust =
          p_adjust_method
      )
    ) %>%

      mutate(
        Layer =
          layer,
        Metric =
          metric,
        Factor =
          "Clade"
      )


    alpha_pairwise[[pair_counter]] <- pair_clade


    pair_counter <-
      pair_counter +
      1


    ########################################################
    # EMM: STAGE
    ########################################################

    emm_stage <- emmeans::emmeans(
      fit,
      ~ Stage,
      vcov. =
        vcov_hc3
    )


    pair_stage <- as.data.frame(
      pairs(
        emm_stage,
        adjust =
          p_adjust_method
      )
    ) %>%

      mutate(
        Layer =
          layer,
        Metric =
          metric,
        Factor =
          "Stage"
      )


    alpha_pairwise[[pair_counter]] <- pair_stage


    pair_counter <-
      pair_counter +
      1


    ########################################################
    # SIMPLE EFFECTS IF INTERACTION
    ########################################################

    if (
      use_interaction
    ) {

      simple_clade <- as.data.frame(

        pairs(

          emmeans::emmeans(
            fit,
            ~ Clade |
              Stage,
            vcov. =
              vcov_hc3
          ),

          adjust =
            p_adjust_method
        )
      ) %>%

        mutate(
          Layer =
            layer,
          Metric =
            metric,
          Factor =
            "Clade_within_Stage"
        )


      alpha_pairwise[[pair_counter]] <- simple_clade


      pair_counter <-
        pair_counter +
        1


      simple_stage <- as.data.frame(

        pairs(

          emmeans::emmeans(
            fit,
            ~ Stage |
              Clade,
            vcov. =
              vcov_hc3
          ),

          adjust =
            p_adjust_method
        )
      ) %>%

        mutate(
          Layer =
            layer,
          Metric =
            metric,
          Factor =
            "Stage_within_Clade"
        )


      alpha_pairwise[[pair_counter]] <- simple_stage


      pair_counter <-
        pair_counter +
        1
    }
  }
}


options(
  contrasts =
    old_contrasts$contrasts
)


alpha_test_df <- bind_rows(
  alpha_tests
)


alpha_pair_df <- bind_rows(
  alpha_pairwise
)


write.csv(
  alpha_test_df,
  file.path(
    outdir,
    "01_ALPHA_joint_model_tests_all_layers.csv"
  ),
  row.names = FALSE
)


write.csv(
  alpha_pair_df,
  file.path(
    outdir,
    "01_ALPHA_EMMeans_pairwise_BH_all_layers.csv"
  ),
  row.names = FALSE
)


############################################################
# 9. EXTRACT ALPHA P VALUES FOR PLOT LABELS
############################################################

pcol_alpha <- grep(
  "Pr\\(",
  names(
    alpha_test_df
  ),
  value = TRUE
)[1]


alpha_effect_labels <- alpha_test_df %>%

  filter(
    Term %in%
      c(
        "Clade",
        "Stage",
        "Clade:Stage"
      )
  ) %>%

  transmute(
    Layer,
    Metric,
    Term,
    P =
      .data[[pcol_alpha]]
  ) %>%

  pivot_wider(
    names_from =
      Term,
    values_from =
      P
  )


if (
  !"Clade:Stage" %in%
    names(
      alpha_effect_labels
    )
) {

  alpha_effect_labels[["Clade:Stage"]] <- NA_real_
}


alpha_effect_labels <- alpha_effect_labels %>%

  mutate(

    Label =
      if (
        use_interaction
      ) {

        paste0(
          "Clade P ",
          vapply(
            Clade,
            format_p,
            character(
              1
            )
          ),
          "\nStage P ",
          vapply(
            Stage,
            format_p,
            character(
              1
            )
          ),
          "\nC×S P ",
          vapply(
            `Clade:Stage`,
            format_p,
            character(
              1
            )
          )
        )

      } else {

        paste0(
          "Clade P ",
          vapply(
            Clade,
            format_p,
            character(
              1
            )
          ),
          "\nStage P ",
          vapply(
            Stage,
            format_p,
            character(
              1
            )
          )
        )
      }
  )


############################################################
# 10. ALPHA PANEL — CLADE VIEW
############################################################

stage_levels <- levels(
  droplevels(
    meta$Stage
  )
)


stage_shapes <- setNames(
  stage_shape_pool[
    seq_along(
      stage_levels
    )
  ],
  stage_levels
)


alpha_clade_label_pos <- alpha_long %>%

  group_by(
    Layer,
    Metric
  ) %>%

  summarise(
    X =
      1,

    Y =
      max(
        Value,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  left_join(
    alpha_effect_labels %>%
      select(
        Layer,
        Metric,
        Label
      ),
    by =
      c(
        "Layer",
        "Metric"
      )
  )


p_alpha_clade <- ggplot(
  alpha_long,
  aes(
    x =
      Clade,
    y =
      Value
  )
) +

  geom_boxplot(
    aes(
      colour =
        Clade
    ),
    width =
      0.55,
    outlier.shape =
      NA,
    fill =
      "white",
    linewidth =
      0.6
  ) +

  geom_jitter(
    aes(
      fill =
        Clade,
      shape =
        Stage
    ),
    width =
      0.10,
    size =
      2,
    alpha =
      0.80,
    colour =
      "black",
    stroke =
      0.25
  ) +

  geom_text(
    data =
      alpha_clade_label_pos,
    aes(
      x =
        X,
      y =
        Y,
      label =
        Label
    ),
    inherit.aes =
      FALSE,
    hjust =
      0,
    vjust =
      1,
    size =
      2.7
  ) +

  facet_grid(
    Metric ~ Layer,
    scales =
      "free_y"
  ) +

  scale_colour_manual(
    values =
      nature_clade
  ) +

  scale_fill_manual(
    values =
      nature_clade
  ) +

  scale_shape_manual(
    values =
      stage_shapes
  ) +

  labs(
    x =
      "Clade",
    y =
      "Alpha diversity",
    colour =
      "Clade",
    fill =
      "Clade",
    shape =
      "Stage",
    title =
      "Alpha diversity: Clade view"
  ) +

  theme_nature(
    base_size = 9.5
  )


############################################################
# 11. ALPHA PANEL — STAGE VIEW
############################################################

p_alpha_stage <- ggplot(
  alpha_long,
  aes(
    x =
      Stage,
    y =
      Value
  )
) +

  geom_boxplot(
    width =
      0.55,
    outlier.shape =
      NA,
    fill =
      "white",
    colour =
      "black",
    linewidth =
      0.6
  ) +

  geom_jitter(
    aes(
      fill =
        Clade,
      shape =
        Stage
    ),
    width =
      0.10,
    size =
      2,
    alpha =
      0.80,
    colour =
      "black",
    stroke =
      0.25
  ) +

  facet_grid(
    Metric ~ Layer,
    scales =
      "free_y"
  ) +

  scale_fill_manual(
    values =
      nature_clade
  ) +

  scale_shape_manual(
    values =
      stage_shapes
  ) +

  labs(
    x =
      "Stage",
    y =
      "Alpha diversity",
    fill =
      "Clade",
    shape =
      "Stage",
    title =
      "Alpha diversity: Stage view"
  ) +

  theme_nature(
    base_size =
      9.5
  ) +

  theme(
    legend.position =
      "none"
  )


############################################################
# 12. FINAL ALPHA PANEL
############################################################

alpha_panel <-

  p_alpha_clade /

  p_alpha_stage +

  plot_layout(
    heights =
      c(
        1,
        1
      ),
    guides =
      "collect"
  ) +

  plot_annotation(
    tag_levels =
      "A"
  ) &

  theme(
    legend.position =
      "right"
  )


save_plot(
  alpha_panel,
  "Alpha_diversity_multilayer_panel",
  width = 14,
  height = 9
)


############################################################
# PART II. MULTI-LAYER BETA DIVERSITY
############################################################

############################################################
# 13. GENERIC BETA ANALYSIS FUNCTION
############################################################

run_beta_layer <- function(
    file,
    layer,
    meta,
    prevalence_cut
) {

  if (layer == "Gene") {
    
    mat <- read_gene_beta_fast(
      file = file,
      sample_ids = meta$SampleID
    )
    
  } else {
    
    mat <- read_beta_file(
      file = file,
      sample_ids = meta$SampleID,
      layer = layer
    )
  }


  if (
    is.null(
      mat
    )
  ) {

    return(
      NULL
    )
  }


  ##########################################################
  # ALIGN SAMPLE IDS
  ##########################################################

  common <- intersect(
    meta$SampleID,
    rownames(
      mat
    )
  )


  sub_meta <- meta %>%

    filter(
      SampleID %in%
        common
    )


  mat <- mat[
    sub_meta$SampleID,
    ,
    drop = FALSE
  ]


  stopifnot(
    all(
      rownames(
        mat
      ) ==
        sub_meta$SampleID
    )
  )


  ##########################################################
  # REMOVE ZERO FEATURES
  ##########################################################

  mat <- mat[
    ,
    colSums(
      mat
    ) >
      0,
    drop = FALSE
  ]


  ##########################################################
  # PREVALENCE FILTER
  ##########################################################

  prevalence <- colMeans(
    mat >
      0
  )


  keep <-
    prevalence >=
    prevalence_cut


  if (
    sum(
      keep
    ) >=
      2
  ) {

    mat <- mat[
      ,
      keep,
      drop = FALSE
    ]
  }


  ##########################################################
  # BRAY-CURTIS
  ##########################################################

  bray <- vegan::vegdist(
    mat,
    method = "bray"
  )
  if (layer == "Gene") {
    
    mat_for_pairwise <- NULL
    
    gc()
  }


  ##########################################################
  # PCoA
  ##########################################################

  pco <- ape::pcoa(
    bray,
    correction =
      "lingoes"
  )


  points <- data.frame(

    SampleID =
      rownames(
        pco$vectors
      ),

    PCoA1 =
      pco$vectors[
        ,
        1
      ],

    PCoA2 =
      pco$vectors[
        ,
        2
      ]
  ) %>%

    left_join(
      sub_meta,
      by =
        "SampleID"
    )


  pc1 <-
    100 *
    pco$values$Relative_eig[
      1
    ]


  pc2 <-
    100 *
    pco$values$Relative_eig[
      2
    ]


  ##########################################################
  # PERMANOVA MODEL
  ##########################################################

  if (
    use_interaction
  ) {

    form <-
      bray ~
      Clade *
      Stage

  } else {

    form <-
      bray ~
      Clade +
      Stage
  }


  set.seed(
    20260813
  )


  perm_seq <- vegan::adonis2(
    form,
    data =
      sub_meta,
    permutations =
      nperm,
    by =
      "terms"
  )


  set.seed(
    20260813
  )


  perm_margin <- vegan::adonis2(
    form,
    data =
      sub_meta,
    permutations =
      nperm,
    by =
      "margin"
  )


  ##########################################################
  # PERMDISP CLADE
  ##########################################################

  bd_clade <- vegan::betadisper(
    bray,
    group =
      sub_meta$Clade,
    type =
      "median"
  )


  set.seed(
    20260813
  )


  pt_clade <- vegan::permutest(
    bd_clade,
    permutations =
      nperm
  )


  ##########################################################
  # PERMDISP STAGE
  ##########################################################

  bd_stage <- vegan::betadisper(
    bray,
    group =
      sub_meta$Stage,
    type =
      "median"
  )


  set.seed(
    20260813
  )


  pt_stage <- vegan::permutest(
    bd_stage,
    permutations =
      nperm
  )


  ##########################################################
  # SUMMARY TABLE
  ##########################################################

  margin_df <- as.data.frame(
    perm_margin
  )


  margin_df$Term <-
    rownames(
      margin_df
    )


  effect_df <- margin_df %>%

    filter(
      Term %in%
        c(
          "Clade",
          "Stage",
          "Clade:Stage"
        )
    ) %>%

    transmute(

      Layer =
        layer,

      Term,

      R2 =
        R2,

      P =
        `Pr(>F)`
    )


  dispersion_df <- data.frame(

    Layer =
      layer,

    Factor =
      c(
        "Clade",
        "Stage"
      ),

    F =
      c(
        pt_clade$tab[
          1,
          "F"
        ],
        pt_stage$tab[
          1,
          "F"
        ]
      ),

    P =
      c(
        pt_clade$tab[
          1,
          "Pr(>F)"
        ],
        pt_stage$tab[
          1,
          "Pr(>F)"
        ]
      )
  )


  ##########################################################
  # WRITE LAYER-SPECIFIC RESULTS
  ##########################################################

  layer_dir <- file.path(
    outdir,
    paste0(
      "02_BETA_",
      layer
    )
  )


  dir.create(
    layer_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  write.csv(
    as.matrix(
      bray
    ),
    file.path(
      layer_dir,
      "BrayCurtis_distance.csv"
    )
  )


  write.csv(
    points,
    file.path(
      layer_dir,
      "PCoA_coordinates.csv"
    ),
    row.names = FALSE
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


  capture.output(
    pt_clade,
    file = file.path(
      layer_dir,
      "PERMDISP_Clade.txt"
    )
  )


  capture.output(
    pt_stage,
    file = file.path(
      layer_dir,
      "PERMDISP_Stage.txt"
    )
  )


  ##########################################################
  # RETURN
  ##########################################################

  list(

    Layer =
      layer,

    Matrix =
      mat,

    Metadata =
      sub_meta,

    Distance =
      bray,

    Points =
      points,

    PC1 =
      pc1,

    PC2 =
      pc2,

    PERMANOVA_seq =
      perm_seq,

    PERMANOVA_margin =
      perm_margin,

    Effect_summary =
      effect_df,

    Dispersion_summary =
      dispersion_df
  )
}


############################################################
# 14. RUN ALL FIVE BETA LAYERS
############################################################

beta_results <- list()


for (
  layer in
    names(
      beta_files
    )
) {

  if (
    !file.exists(
      beta_files[[layer]]
    )
  ) {

    warning(
      "Beta file not found; skipping ",
      layer,
      ": ",
      beta_files[[layer]]
    )

    next
  }


  beta_results[[layer]] <- run_beta_layer(

    file =
      beta_files[[layer]],

    layer =
      layer,

    meta =
      meta,

    prevalence_cut =
      beta_prevalence_cut[[layer]]
  )
}


if (
  length(
    beta_results
  ) ==
    0
) {

  stop(
    "No beta-diversity datasets were successfully analysed."
  )
}


############################################################
# 15. COMBINE BETA SUMMARY TABLES
############################################################

beta_effects <- bind_rows(
  lapply(
    beta_results,
    function(x) {
      x$Effect_summary
    }
  )
)


beta_dispersion <- bind_rows(
  lapply(
    beta_results,
    function(x) {
      x$Dispersion_summary
    }
  )
)


write.csv(
  beta_effects,
  file.path(
    outdir,
    "02_BETA_marginal_PERMANOVA_effects_all_layers.csv"
  ),
  row.names = FALSE
)


write.csv(
  beta_dispersion,
  file.path(
    outdir,
    "02_BETA_PERMDISP_all_layers.csv"
  ),
  row.names = FALSE
)


############################################################
# 16. BUILD PCoA PANELS
############################################################

pcoa_plots <- list()


for (
  layer in
    names(
      beta_results
    )
) {

  x <- beta_results[[layer]]


  sub_stage_levels <- levels(
    droplevels(
      x$Metadata$Stage
    )
  )


  sub_shapes <- setNames(
    stage_shape_pool[
      seq_along(
        sub_stage_levels
      )
    ],
    sub_stage_levels
  )


  p <- ggplot(
    x$Points,
    aes(
      PCoA1,
      PCoA2,
      colour =
        Clade,
      shape =
        Stage
    )
  ) +

    stat_ellipse(
      aes(
        group =
          Clade,
        fill =
          Clade
      ),
      geom =
        "polygon",
      alpha =
        0.08,
      colour =
        NA,
      level =
        0.95
    ) +

    geom_point(
      aes(
        fill =
          Clade
      ),
      size =
        2.7,
      stroke =
        0.3,
      alpha =
        0.9
    ) +

    geom_vline(
      xintercept =
        0,
      linetype =
        3,
      colour =
        "grey75"
    ) +

    geom_hline(
      yintercept =
        0,
      linetype =
        3,
      colour =
        "grey75"
    ) +

    scale_colour_manual(
      values =
        nature_clade
    ) +

    scale_fill_manual(
      values =
        nature_clade
    ) +

    scale_shape_manual(
      values =
        sub_shapes
    ) +

    labs(

      x =
        paste0(
          "PCoA1 (",
          round(
            x$PC1,
            1
          ),
          "%)"
        ),

      y =
        paste0(
          "PCoA2 (",
          round(
            x$PC2,
            1
          ),
          "%)"
        ),

      title =
        layer,

      colour =
        "Clade",

      fill =
        "Clade",

      shape =
        "Stage"
    ) +

    theme_nature(
      base_size =
        9
    )


  pcoa_plots[[layer]] <- p
}


############################################################
# 17. MARGINAL PERMANOVA R2 SUMMARY
############################################################

beta_effects_plot <- beta_effects %>%

  mutate(

    Term =
      factor(
        Term,
        levels =
          c(
            "Clade",
            "Stage",
            "Clade:Stage"
          )
      ),

    Layer =
      factor(
        Layer,
        levels =
          names(
            layer_colours
          )
      ),

    P_label =
      ifelse(
        P <
          0.0001,
        "<0.0001",
        format(
          signif(
            P,
            2
          ),
          scientific =
            TRUE
        )
      )
  )


p_beta_r2 <- ggplot(
  beta_effects_plot,
  aes(
    x =
      Layer,
    y =
      R2,
    fill =
      Term
  )
) +

  geom_col(
    position =
      position_dodge(
        width =
          0.78
      ),
    width =
      0.70,
    colour =
      "black",
    linewidth =
      0.2
  ) +

  scale_fill_manual(
    values =
      c(
        "Clade" =
          "#3C5488",
        "Stage" =
          "#E64B35",
        "Clade:Stage" =
          "#00A087"
      ),
    drop =
      FALSE
  ) +

  scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy =
          1
      )
  ) +

  labs(
    x =
      NULL,
    y =
      "Marginal PERMANOVA R²",
    fill =
      "Effect",
    title =
      "Independent Clade and Stage effects"
  ) +

  theme_nature(
    base_size =
      9.5
  )


############################################################
# 18. PERMDISP SUMMARY
############################################################

beta_disp_plot <- beta_dispersion %>%

  mutate(
    Significant =
      P <
      0.05
  )


p_beta_disp <- ggplot(
  beta_disp_plot,
  aes(
    x =
      Layer,
    y =
      F,
    fill =
      Factor
  )
) +

  geom_col(
    position =
      position_dodge(
        width =
          0.78
      ),
    width =
      0.70,
    colour =
      "black",
    linewidth =
      0.2
  ) +

  geom_text(
    aes(
      label =
        ifelse(
          Significant,
          "*",
          ""
        ),
      group =
        Factor
    ),
    position =
      position_dodge(
        width =
          0.78
      ),
    vjust =
      -0.35,
    size =
      4
  ) +

  scale_fill_manual(
    values =
      c(
        "Clade" =
          "#3C5488",
        "Stage" =
          "#E64B35"
      )
  ) +

  labs(
    x =
      NULL,
    y =
      "PERMDISP F",
    fill =
      "Grouping factor",
    title =
      "Homogeneity of multivariate dispersion"
  ) +

  theme_nature(
    base_size =
      9.5
  )


############################################################
# 19. FINAL BETA PANEL
############################################################

# Arrange all available PCoAs automatically.
pcoa_combined <- wrap_plots(
  pcoa_plots,
  ncol =
    min(
      3,
      length(
        pcoa_plots
      )
    ),
  guides =
    "collect"
)


beta_summary_row <-
  p_beta_r2 |
  p_beta_disp


beta_panel <-

  pcoa_combined /

  beta_summary_row +

  plot_layout(
    heights =
      c(
        2,
        1
      ),
    guides =
      "collect"
  ) +

  plot_annotation(
    tag_levels =
      "A"
  ) &

  theme(
    legend.position =
      "right"
  )


save_plot(
  beta_panel,
  "Beta_diversity_multilayer_panel",
  width = 14,
  height = 11
)


############################################################
# PART III. OPTIONAL PAIRWISE BETA TESTS
############################################################

pairwise_permanova_from_distance <- function(
    full_distance,
    meta,
    factor_name,
    permutations = 9999
) {
  
  levs <- levels(
    droplevels(
      meta[[factor_name]]
    )
  )
  
  if (length(levs) < 2) {
    return(data.frame())
  }
  
  pairs <- combn(
    levs,
    2,
    simplify = FALSE
  )
  
  full_mat <- as.matrix(
    full_distance
  )
  
  out <- lapply(
    pairs,
    function(pair) {
      
      idx <- meta[[factor_name]] %in% pair
      
      sub_meta <- droplevels(
        meta[
          idx,
          ,
          drop = FALSE
        ]
      )
      
      sample_ids <- sub_meta$SampleID
      
      sub_dist <- as.dist(
        full_mat[
          sample_ids,
          sample_ids,
          drop = FALSE
        ]
      )
      
      form <- as.formula(
        paste0(
          "sub_dist ~ ",
          factor_name
        )
      )
      
      set.seed(20260813)
      
      ad <- vegan::adonis2(
        form,
        data = sub_meta,
        permutations = permutations
      )
      
      data.frame(
        Factor = factor_name,
        Comparison = paste(
          pair,
          collapse = " vs "
        ),
        R2 = ad$R2[1],
        P = ad[["Pr(>F)"]][1]
      )
    }
  )
  
  out <- bind_rows(
    out
  )
  
  out$Q <- p.adjust(
    out$P,
    method = "BH"
  )
  
  out
}


pairwise_beta_all <- list()
pw_counter <- 1


for (
  layer in
    names(
      beta_results
    )
) {

  x <- beta_results[[layer]]


  pw_clade <- pairwise_permanova_from_distance(
    x$Distance,
    x$Metadata,
    "Clade",
    nperm
  )
  
  pw_stage <- pairwise_permanova_from_distance(
    x$Distance,
    x$Metadata,
    "Stage",
    nperm
  )


  pairwise_beta_all[[pw_counter]] <- bind_rows(
    pw_clade,
    pw_stage
  ) %>%

    mutate(
      Layer =
        layer
    )


  pw_counter <-
    pw_counter +
    1
}


write.csv(
  bind_rows(
    pairwise_beta_all
  ),
  file.path(
    outdir,
    "02_BETA_pairwise_PERMANOVA_BH_all_layers.csv"
  ),
  row.names = FALSE
)


############################################################
# PART IV. REPRODUCIBILITY
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)


notes <- c(

  paste0(
    "Joint model = ",
    model_text
  ),

  "All five requested layers are supported: Gene, COG, Family, KO and Species.",

  "Missing layer files are skipped with a warning rather than stopping all analyses.",

  "Alpha analysis uses HC3 robust covariance and BH-adjusted estimated marginal means.",

  "The original unadjusted LSD.test workflow is not used.",

  "Beta analysis uses Bray-Curtis dissimilarity and Lingoes-corrected PCoA.",

  paste0(
    "PERMANOVA/PERMDISP permutations = ",
    nperm
  ),

  "Marginal PERMANOVA is used for the independent Clade and Stage effects.",

  "PERMDISP is evaluated separately for Clade and Stage.",

  "Beta input orientation is automatically determined by matching metadata SampleIDs.",

  "Every abundance matrix is explicitly reordered to metadata SampleID before statistics.",

  "Final Alpha and Beta diversity figures are exported as separate multilayer panels."
)


writeLines(
  notes,
  file.path(
    outdir,
    "00_analysis_notes.txt"
  )
)


cat("\n===========================================================\n")

cat("MULTI-LAYER ALPHA/BETA ANALYSIS COMPLETED\n")

cat("===========================================================\n")

cat(
  "Joint model:",
  model_text,
  "\n"
)

cat(
  "Alpha layers analysed:",
  paste(
    unique(
      alpha_long$Layer
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Beta layers analysed:",
  paste(
    names(
      beta_results
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nAlpha panel:\n",
  file.path(
    outdir,
    "Alpha_diversity_multilayer_panel.pdf"
  ),
  "\n"
)

cat(
  "\nBeta panel:\n",
  file.path(
    outdir,
    "Beta_diversity_multilayer_panel.pdf"
  ),
  "\n"
)

cat("===========================================================\n")
