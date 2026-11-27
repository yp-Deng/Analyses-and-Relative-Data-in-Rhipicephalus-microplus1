############################################################
# Fig_S22_NatureCom_Final.R
#
# Nature Communications-style Supplementary Figure
#
# Panel a:
#   VF gene presence/absence heatmap
#
# Panel b:
#   VF repertoire size by genome
#
# Panel c:
#   Genus-specific VF enrichment:
#   A1_final vs Genus1 references
#   A3_final vs Genus2 references
#   G1_final vs Genus3 references
#
# Group membership is determined automatically from
# the three IQ-TREE treefiles.
############################################################


############################################################
# 0. Packages
############################################################

library(tidyverse)
library(ape)
library(pheatmap)
library(patchwork)
library(grid)
library(ggrepel)
library(scales)


############################################################
# 1. Working directory and paths
############################################################

WD <- "G:/comparetive_genomics/comparative_genomics"

setwd(WD)


TREE_DIR <- "three_genus_results"

VF_FILE <-
  "VF_BGC_results_V4/VF_summary/VF_gene_matrix_strain.csv"


OUTDIR <- "Fig_S22_NatureCom_Final"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


TABLEDIR <- file.path(
  OUTDIR,
  "tables"
)

dir.create(
  TABLEDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 2. Nature-style palette
############################################################

COL_GENUS <- c(
  "Genus1" = "#3C5488",
  "Genus2" = "#00A087",
  "Genus3" = "#E64B35"
)

COL_STATUS <- c(
  "Reference" = "#BDBDBD",
  "Novel"     = "#F39B7F"
)

COL_VF <- c(
  "0" = "#FFFFFF",
  "1" = "#C44E52"
)


############################################################
# 3. Novel species mapping
############################################################

novel_map <- c(
  "Genus1" = "A1_final",
  "Genus2" = "A3_final",
  "Genus3" = "G1_final"
)


############################################################
# 4. Automatically read tree tips
############################################################

get_tree_tips <- function(genus) {
  
  treefiles <- list.files(
    file.path(TREE_DIR, genus),
    pattern = "\\.treefile$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(treefiles) == 0) {
    stop(
      paste0(
        "No treefile found for ",
        genus
      )
    )
  }
  
  if (length(treefiles) > 1) {
    
    message(
      genus,
      ": multiple treefiles found; using ",
      treefiles[1]
    )
    
  }
  
  tr <- ape::read.tree(
    treefiles[1]
  )
  
  if (is.null(tr$tip.label)) {
    stop(
      paste0(
        "No tip labels detected in ",
        genus
      )
    )
  }
  
  data.frame(
    Genome = tr$tip.label,
    Genus = genus,
    Tree_order = seq_along(tr$tip.label),
    stringsAsFactors = FALSE
  )
}


tree_groups <- bind_rows(
  get_tree_tips("Genus1"),
  get_tree_tips("Genus2"),
  get_tree_tips("Genus3")
)


############################################################
# Check duplicate genomes between trees
############################################################

duplicated_tips <- tree_groups %>%
  count(Genome) %>%
  filter(n > 1)


if (nrow(duplicated_tips) > 0) {
  
  stop(
    paste(
      "Some genomes occur in more than one tree:",
      paste(
        duplicated_tips$Genome,
        collapse = ", "
      )
    )
  )
  
}


############################################################
# Check novel species
############################################################

for (g in names(novel_map)) {
  
  nsp <- novel_map[[g]]
  
  expected <- tree_groups %>%
    filter(
      Genus == g,
      Genome == nsp
    )
  
  if (nrow(expected) == 0) {
    
    stop(
      paste0(
        nsp,
        " was not found in ",
        g,
        " tree."
      )
    )
    
  }
  
}


tree_groups <- tree_groups %>%
  mutate(
    Status = ifelse(
      Genome %in% novel_map,
      "Novel",
      "Reference"
    )
  )


write.csv(
  tree_groups,
  file.path(
    TABLEDIR,
    "Tree_based_genome_groups.csv"
  ),
  row.names = FALSE
)


cat("\nTree-based groups:\n")

print(
  table(
    tree_groups$Genus
  )
)


############################################################
# 5. Load strain-level VF matrix
############################################################

vf <- read.csv(
  VF_FILE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

vf <- as.matrix(vf)

storage.mode(vf) <- "numeric"


############################################################
# Basic matrix checks
############################################################

if (!all(vf %in% c(0, 1))) {
  stop(
    "VF matrix contains values other than 0/1."
  )
}


missing_tree_genomes <- setdiff(
  tree_groups$Genome,
  colnames(vf)
)


if (length(missing_tree_genomes) > 0) {
  
  stop(
    paste(
      "These tree genomes are absent from the VF matrix:",
      paste(
        missing_tree_genomes,
        collapse = ", "
      )
    )
  )
  
}


############################################################
# Keep only genomes represented in the three trees
############################################################

vf <- vf[
  ,
  tree_groups$Genome,
  drop = FALSE
]


############################################################
# Remove genes absent from every selected genome
############################################################

vf <- vf[
  rowSums(vf) > 0,
  ,
  drop = FALSE
]


cat(
  "\nVF matrix used:",
  nrow(vf),
  "genes x",
  ncol(vf),
  "genomes\n"
)


############################################################
# 6. Genus/tree-based genome order
############################################################

genus_order <- c(
  "Genus1",
  "Genus2",
  "Genus3"
)


tree_groups$Genus <- factor(
  tree_groups$Genus,
  levels = genus_order
)


tree_groups <- tree_groups %>%
  arrange(
    Genus,
    Tree_order
  )


genome_order <- tree_groups$Genome


vf <- vf[
  ,
  genome_order,
  drop = FALSE
]


############################################################
# 7. Panel a
# VF presence/absence heatmap
############################################################

annotation_col <- tree_groups %>%
  select(
    Genome,
    Genus,
    Status
  ) %>%
  column_to_rownames(
    "Genome"
  )


############################################################
# Optional:
# Put highly variable VF genes first
############################################################

vf_freq <- rowMeans(vf)

vf_variability <- vf_freq *
  (1 - vf_freq)


row_order <- order(
  vf_variability,
  decreasing = TRUE
)

vf_heat <- vf[
  row_order,
  ,
  drop = FALSE
]


############################################################
# Novel genome labels
############################################################

display_names <- colnames(vf_heat)

display_names[
  display_names %in% novel_map
] <- paste0(
  display_names[
    display_names %in% novel_map
  ],
  " *"
)


############################################################
# Draw pheatmap silently as grob
############################################################

ph <- pheatmap(
  vf_heat,
  
  color = c(
    "#FFFFFF",
    "#C44E52"
  ),
  
  breaks = c(
    -0.5,
    0.5,
    1.5
  ),
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_rownames = FALSE,
  show_colnames = TRUE,
  
  labels_col = display_names,
  
  annotation_col = annotation_col,
  
  annotation_colors = list(
    Genus = COL_GENUS,
    Status = COL_STATUS
  ),
  
  annotation_names_col = FALSE,
  
  border_color = NA,
  
  fontsize = 7,
  fontsize_col = 6,
  
  angle_col = 90,
  
  legend_breaks = c(
    0,
    1
  ),
  
  legend_labels = c(
    "Absent",
    "Present"
  ),
  
  main = "",
  
  silent = TRUE
)


pA <- wrap_elements(
  full = ph$gtable
)


############################################################
# Save panel a separately
############################################################

pdf(
  file.path(
    OUTDIR,
    "Fig_S22a_VF_heatmap.pdf"
  ),
  width = 10,
  height = 8
)

grid.newpage()

grid.draw(
  ph$gtable
)

dev.off()


############################################################
# 8. Panel b
# VF repertoire size
############################################################

vf_number <- data.frame(
  Genome = colnames(vf),
  VF_number = colSums(vf),
  stringsAsFactors = FALSE
)


vf_number <- vf_number %>%
  left_join(
    tree_groups,
    by = "Genome"
  )


############################################################
# Check classification
############################################################

if (any(is.na(vf_number$Genus))) {
  
  warning(
    "Some VF genomes were not matched to tree tips: ",
    paste(
      vf_number$Genome[is.na(vf_number$Genus)],
      collapse = ", "
    )
  )
}


############################################################
# Keep tree-classified genomes
############################################################

vf_number <- vf_number %>%
  filter(
    !is.na(Genus)
  )


############################################################
# Order genomes within each genus by VF number
############################################################

vf_number <- vf_number %>%
  mutate(
    Genus = factor(
      Genus,
      levels = c(
        "Genus1",
        "Genus2",
        "Genus3"
      )
    )
  ) %>%
  arrange(
    Genus,
    VF_number
  )


vf_number$Genome_plot <- factor(
  vf_number$Genome,
  levels = vf_number$Genome
)


############################################################
# Novel species table
############################################################

novel_vf <- vf_number %>%
  filter(
    Status == "Novel"
  )


print(novel_vf)


############################################################
# Plot
############################################################

pB <- ggplot(
  vf_number,
  aes(
    x = VF_number,
    y = Genome_plot,
    fill = Status
  )
) +
  
  geom_col(
    width = 0.72
  ) +
  
  ##########################################################
# Highlight novel species
##########################################################

geom_point(
  data = novel_vf,
  aes(
    x = VF_number,
    y = Genome_plot
  ),
  inherit.aes = FALSE,
  shape = 21,
  size = 2.8,
  stroke = 0.7,
  fill = "white",
  colour = "#B2182B"
) +
  
  ##########################################################
# Label the three novel species
#
# IMPORTANT:
# y = Genome_plot is explicitly supplied
##########################################################

geom_text(
  data = novel_vf,
  aes(
    x = VF_number,
    y = Genome_plot,
    label = Genome
  ),
  inherit.aes = FALSE,
  hjust = -0.15,
  size = 2.6,
  fontface = "bold",
  colour = "#B2182B"
) +
  
  facet_grid(
    rows = vars(Genus),
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Reference" = "#BDBDBD",
      "Novel" = "#00A087"
    )
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.23
      )
    )
  ) +
  
  labs(
    x = "Predicted virulence-associated homologues",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 8
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text.y = element_text(
      face = "bold",
      size = 8,
      angle = 0
    ),
    
    axis.text.y = element_text(
      size = 5.8
    ),
    
    axis.text.x = element_text(
      size = 7
    ),
    
    axis.title.x = element_text(
      size = 8
    ),
    
    legend.position = "top",
    
    legend.justification = "left",
    
    panel.spacing.y = grid::unit(
      0.45,
      "lines"
    )
  )


############################################################
# Save Panel b
############################################################

ggsave(
  file.path(
    OUTDIR,
    "Fig_S22b_VF_repertoire.pdf"
  ),
  pB,
  width = 5.2,
  height = 8
)


write.csv(
  vf_number,
  file.path(
    TABLEDIR,
    "VF_number_by_genome.csv"
  ),
  row.names = FALSE
)


############################################################
# 9. Genus-specific Fisher enrichment
############################################################

run_vf_enrichment <- function(
    genus,
    novel_genome,
    vf_matrix,
    groups
) {
  
  genus_genomes <- groups %>%
    filter(
      Genus == genus
    ) %>%
    pull(
      Genome
    )
  
  references <- setdiff(
    genus_genomes,
    novel_genome
  )
  
  
  if (!novel_genome %in% genus_genomes) {
    
    stop(
      paste0(
        novel_genome,
        " is not present in ",
        genus
      )
    )
    
  }
  
  
  if (length(references) < 1) {
    
    stop(
      paste0(
        genus,
        " has no reference genomes."
      )
    )
    
  }
  
  
  results <- vector(
    "list",
    nrow(vf_matrix)
  )
  
  
  for (i in seq_len(nrow(vf_matrix))) {
    
    gene <- rownames(
      vf_matrix
    )[i]
    
    
    novel_present <-
      as.integer(
        vf_matrix[
          i,
          novel_genome
        ] > 0
      )
    
    
    ref_present <-
      sum(
        vf_matrix[
          i,
          references,
          drop = TRUE
        ] > 0
      )
    
    
    ref_absent <-
      length(references) -
      ref_present
    
    
    novel_absent <-
      1 -
      novel_present
    
    
    ########################################################
    # 2 x 2:
    #
    #                Present   Absent
    # Novel             a        b
    # Reference         c        d
    ########################################################
    
    tab <- matrix(
      c(
        novel_present,
        novel_absent,
        ref_present,
        ref_absent
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        Group = c(
          "Novel",
          "Reference"
        ),
        VF = c(
          "Present",
          "Absent"
        )
      )
    )
    
    
    ft <- fisher.test(
      tab,
      alternative = "greater"
    )
    
    
    ########################################################
    # Haldane-Anscombe corrected OR
    #
    # avoids zero / infinite values for plotting
    ########################################################
    
    a <- novel_present + 0.5
    b <- novel_absent + 0.5
    c <- ref_present + 0.5
    d <- ref_absent + 0.5
    
    
    corrected_or <-
      (a * d) /
      (b * c)
    
    
    results[[i]] <- data.frame(
      
      Genus = genus,
      
      Novel_genome = novel_genome,
      
      VF_gene = gene,
      
      Novel_present = novel_present,
      
      Reference_present = ref_present,
      
      Reference_total = length(
        references
      ),
      
      Reference_frequency =
        ref_present /
        length(
          references
        ),
      
      Fisher_OR =
        unname(
          ft$estimate
        ),
      
      Corrected_OR =
        corrected_or,
      
      log2_OR =
        log2(
          corrected_or
        ),
      
      P_value =
        ft$p.value,
      
      stringsAsFactors = FALSE
      
    )
    
  }
  
  
  result <- bind_rows(
    results
  )
  
  
  result$FDR <- p.adjust(
    result$P_value,
    method = "BH"
  )
  
  
  result
  
}


############################################################
# Run all three independent comparisons
############################################################

enrich_all <- bind_rows(
  
  run_vf_enrichment(
    genus = "Genus1",
    novel_genome = "A1_final",
    vf_matrix = vf,
    groups = tree_groups
  ),
  
  run_vf_enrichment(
    genus = "Genus2",
    novel_genome = "A3_final",
    vf_matrix = vf,
    groups = tree_groups
  ),
  
  run_vf_enrichment(
    genus = "Genus3",
    novel_genome = "G1_final",
    vf_matrix = vf,
    groups = tree_groups
  )
  
)


############################################################
# Significance / descriptive categories
############################################################

enrich_all <- enrich_all %>%
  mutate(
    
    neglog10FDR =
      -log10(
        pmax(
          FDR,
          .Machine$double.xmin
        )
      ),
    
    Significance = case_when(
      
      FDR < 0.05 &
        Corrected_OR > 2
      ~ "FDR < 0.05",
      
      P_value < 0.05 &
        Corrected_OR > 2
      ~ "Nominal P < 0.05",
      
      TRUE
      ~ "Not significant"
      
    )
    
  )


write.csv(
  enrich_all,
  file.path(
    TABLEDIR,
    "VF_genus_specific_enrichment_all.csv"
  ),
  row.names = FALSE
)


############################################################
# Separate result tables
############################################################

for (g in genus_order) {
  
  write.csv(
    
    enrich_all %>%
      filter(
        Genus == g
      ),
    
    file.path(
      TABLEDIR,
      paste0(
        g,
        "_VF_enrichment.csv"
      )
    ),
    
    row.names = FALSE
    
  )
  
}


############################################################
# 10. Select labels for Panel c
#
# Prefer FDR-significant genes.
# If none, label 5 genes with smallest P in each genus.
############################################################

label_data <- enrich_all %>%
  group_by(Genus) %>%
  group_modify(
    ~ {
      
      sig <- .x %>%
        filter(
          FDR < 0.05,
          Corrected_OR > 2
        ) %>%
        arrange(
          FDR,
          desc(
            Corrected_OR
          )
        )
      
      if (nrow(sig) > 0) {
        
        head(
          sig,
          8
        )
        
      } else {
        
        .x %>%
          filter(
            Novel_present == 1
          ) %>%
          arrange(
            P_value,
            desc(
              Corrected_OR
            )
          ) %>%
          head(
            5
          )
        
      }
      
    }
  ) %>%
  ungroup()


############################################################
# Panel c volcano/effect-size plot
############################################################

ymax_data <- max(
  enrich_all$neglog10FDR,
  na.rm = TRUE
)

# 给标签预留一些空间
ymax_plot <- max(
  0.15,
  ymax_data * 1.55
)


############################################################
# Label only the most informative genes
#
# 每个属：
# 优先 nominal P < 0.05；
# 如果不足，则按 P 值 + effect size 选前 5 个
############################################################

label_data <- enrich_all %>%
  group_by(Genus) %>%
  group_modify(
    ~ {
      
      sig_nominal <- .x %>%
        filter(
          P_value < 0.05,
          Corrected_OR > 1
        ) %>%
        arrange(
          P_value,
          desc(abs(log2_OR))
        )
      
      if(nrow(sig_nominal) >= 5){
        
        head(
          sig_nominal,
          5
        )
        
      } else {
        
        .x %>%
          filter(
            Novel_present == 1
          ) %>%
          arrange(
            P_value,
            desc(abs(log2_OR))
          ) %>%
          head(5)
        
      }
      
    }
  ) %>%
  ungroup()


############################################################
# Plot
############################################################

############################################################
# Panel c
# Genus-specific VF effect-size / enrichment landscape
# Zoomed version
############################################################


enrich_all <- enrich_all %>%
  mutate(
    
    ########################################################
    # Raw Fisher P for Y axis
    ########################################################
    
    neglog10P =
      -log10(
        pmax(
          P_value,
          .Machine$double.xmin
        )
      ),
    
    
    ########################################################
    # BH-FDR status
    #
    # IMPORTANT:
    # Colour represents multiple-testing corrected status,
    # not raw P-value significance.
    ########################################################
    
    FDR_status = case_when(
      
      FDR < 0.05 &
        Corrected_OR > 2
      ~ "FDR < 0.05, enriched",
      
      FDR < 0.05 &
        Corrected_OR < 0.5
      ~ "FDR < 0.05, depleted",
      
      TRUE
      ~ "FDR ≥ 0.05"
      
    )
    
  )


############################################################
# 2. Check results
############################################################

cat(
  "\nPanel c FDR status:\n"
)

print(
  table(
    enrich_all$Genus,
    enrich_all$FDR_status
  )
)


############################################################
# 3. Select genes for annotation
#
# Strategy:
#
# 1) Prefer FDR < 0.05 genes
# 2) If none are FDR-significant:
#    select nominal P < 0.05 genes
# 3) If still too few:
#    select the most informative genes according to
#    Fisher P + |log2 OR|
############################################################

label_data <- enrich_all %>%
  
  group_by(
    Genus
  ) %>%
  
  group_modify(
    
    ~ {
      
      ######################################################
      # FDR significant
      ######################################################
      
      fdr_sig <- .x %>%
        
        filter(
          FDR < 0.05
        ) %>%
        
        arrange(
          FDR,
          P_value,
          desc(
            abs(
              log2_OR
            )
          )
        )
      
      
      if(
        nrow(
          fdr_sig
        ) > 0
      ){
        
        return(
          head(
            fdr_sig,
            6
          )
        )
        
      }
      
      
      ######################################################
      # Nominal Fisher P significant
      ######################################################
      
      nominal <- .x %>%
        
        filter(
          P_value < 0.05
        ) %>%
        
        arrange(
          P_value,
          desc(
            abs(
              log2_OR
            )
          )
        )
      
      
      if(
        nrow(
          nominal
        ) > 0
      ){
        
        return(
          head(
            nominal,
            6
          )
        )
        
      }
      
      
      ######################################################
      # No significant genes:
      # show strongest lineage-associated candidates
      ######################################################
      
      .x %>%
        
        filter(
          Novel_present == 1
        ) %>%
        
        arrange(
          P_value,
          desc(
            abs(
              log2_OR
            )
          )
        ) %>%
        
        head(
          5
        )
      
    }
    
  ) %>%
  
  ungroup()


############################################################
# Save labelled genes
############################################################

write.csv(
  
  label_data,
  
  file.path(
    TABLEDIR,
    "Fig_S22c_labelled_VFs.csv"
  ),
  
  row.names = FALSE
  
)



############################################################
# 4. Determine Y-axis range automatically
############################################################

ymax_data <- max(
  
  enrich_all$neglog10P,
  
  na.rm = TRUE
  
)


############################################################
# Give annotation labels some additional room
############################################################

ymax_plot <- max(
  
  1.5,
  
  ymax_data * 1.20
  
)


cat(
  "\nMaximum -log10(P):",
  ymax_data,
  "\n"
)


############################################################
# 5. Panel c plot
############################################################

pC <- ggplot(
  
  enrich_all,
  
  aes(
    x = log2_OR,
    y = neglog10P
  )
  
) +
  
  
  ##########################################################
# Fisher P = 0.05 threshold
##########################################################

geom_hline(
  
  yintercept =
    -log10(
      0.05
    ),
  
  linewidth = 0.40,
  
  linetype = 2,
  
  colour = "#777777"
  
) +
  
  
  ##########################################################
# OR = 1
##########################################################

geom_vline(
  
  xintercept = 0,
  
  linewidth = 0.35,
  
  colour = "#777777"
  
) +
  
  
  ##########################################################
# OR = 2 and OR = 0.5
##########################################################

geom_vline(
  
  xintercept = c(
    -1,
    1
  ),
  
  linewidth = 0.30,
  
  linetype = 3,
  
  colour = "#AAAAAA"
  
) +
  
  
  ##########################################################
# VF genes
##########################################################

geom_point(
  
  aes(
    colour = FDR_status
  ),
  
  size = 1.75,
  
  alpha = 0.78
  
) +
  
  
  ##########################################################
# Gene labels
##########################################################

geom_text_repel(
  
  data = label_data,
  
  aes(
    x = log2_OR,
    y = neglog10P,
    label = VF_gene
  ),
  
  inherit.aes = FALSE,
  
  size = 2.15,
  
  max.overlaps = Inf,
  
  box.padding = 0.40,
  
  point.padding = 0.25,
  
  min.segment.length = 0,
  
  segment.size = 0.25,
  
  segment.alpha = 0.65,
  
  seed = 123,
  
  show.legend = FALSE
  
) +
  
  
  ##########################################################
# Three independent genus comparisons
##########################################################

facet_wrap(
  
  ~ Genus,
  
  nrow = 1,
  
  scales = "free_x",
  
  labeller = as_labeller(
    
    c(
      
      "Genus1" =
        "A1_final vs Genus1 references",
      
      "Genus2" =
        "A3_final vs Genus2 references",
      
      "Genus3" =
        "G1_final vs Genus3 references"
      
    )
    
  )
  
) +
  
  
  ##########################################################
# BH-FDR status colours
##########################################################

scale_colour_manual(
  
  values = c(
    
    "FDR < 0.05, enriched" =
      "#B2182B",
    
    "FDR < 0.05, depleted" =
      "#2166AC",
    
    "FDR ≥ 0.05" =
      "#BDBDBD"
    
  ),
  
  drop = FALSE
  
) +
  
  
  ##########################################################
# Zoom to actual data region
##########################################################

coord_cartesian(
  
  ylim = c(
    0,
    ymax_plot
  ),
  
  clip = "off"
  
) +
  
  
  ##########################################################
# Axis expansion
##########################################################

scale_y_continuous(
  
  expand = expansion(
    
    mult = c(
      0.02,
      0.08
    )
    
  )
  
) +
  
  
  scale_x_continuous(
    
    expand = expansion(
      
      mult = c(
        0.06,
        0.10
      )
      
    )
    
  ) +
  
  
  ##########################################################
# Labels
##########################################################

labs(
  
  x = expression(
    log[2] *
      " corrected odds ratio"
  ),
  
  y = expression(
    -log[10] *
      " Fisher " *
      italic(P)
  ),
  
  colour = NULL
  
) +
  
  
  ##########################################################
# Nature Communications style
##########################################################

theme_classic(
  
  base_size = 9
  
) +
  
  
  theme(
    
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        
        size = 8.5,
        
        face = "bold"
        
      ),
    
    
    legend.position =
      "top",
    
    legend.justification =
      "left",
    
    legend.text =
      element_text(
        size = 7.5
      ),
    
    
    axis.title =
      element_text(
        size = 9
      ),
    
    axis.text =
      element_text(
        size = 8
      ),
    
    
    panel.spacing =
      grid::unit(
        1.1,
        "lines"
      ),
    
    
    plot.margin =
      margin(
        8,
        14,
        6,
        6
      )
    
  )



############################################################
# 6. Save Panel c individually
############################################################

ggsave(
  
  file.path(
    OUTDIR,
    "Fig_S22c_VF_Fisher_effect_size.pdf"
  ),
  
  pC,
  
  width = 11,
  
  height = 4.2
  
)



############################################################
# 7. Panel labels
############################################################

pB_final <- pB +
  
  labs(
    tag = "b"
  ) +
  
  theme(
    
    plot.tag =
      element_text(
        
        face = "bold",
        
        size = 12
        
      )
    
  )



pC_final <- pC +
  
  labs(
    tag = "c"
  ) +
  
  theme(
    
    plot.tag =
      element_text(
        
        face = "bold",
        
        size = 12
        
      )
    
  )



############################################################
# Panel a
############################################################

pA_final <- pA



############################################################
# 8. Right side:
#
# Increase Panel c height
#
# Previous:
# heights = c(1.25, 0.75)
#
# FINAL:
# heights = c(1.00, 1.05)
############################################################

right_panel <-
  
  pB_final /
  
  pC_final +
  
  plot_layout(
    
    heights = c(
      
      1.00,
      
      1.05
      
    )
    
  )



############################################################
# 9. Final Figure S22 layout
############################################################

final_fig <-
  
  pA_final |
  
  right_panel



final_fig <-
  
  final_fig +
  
  plot_layout(
    
    widths = c(
      
      1.20,
      
      1.10
      
    )
    
  )



############################################################
# 10. Save final PDF
############################################################

ggsave(
  
  file.path(
    
    OUTDIR,
    
    "Fig_S22_NatureCom_Final.pdf"
    
  ),
  
  final_fig,
  
  width = 15,
  
  height = 10,
  
  units = "in"
  
)



############################################################
# 11. Save final TIFF
############################################################

ggsave(
  
  file.path(
    
    OUTDIR,
    
    "Fig_S22_NatureCom_Final.tiff"
    
  ),
  
  final_fig,
  
  width = 15,
  
  height = 10,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw"
  
)


