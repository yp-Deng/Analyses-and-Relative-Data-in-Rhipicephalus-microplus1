############################################################
# R. microplus genus-level advanced microbiome analyses
#
# INPUTS
#   Genus.xlsx        : relative abundance, samples x genera
#   Genus_Count.xlsx  : raw genus counts, samples x genera
#   metadata.xlsx     : SampleID + Clade + Stage + environment
#
# Analyses
#   1. Procrustes / PROTEST
#      - microbiome Bray-Curtis ordination vs environmental ordination
#      - microbiome vs geographic-distance ordination
#   2. Microbial association network
#      - SPIEC-EASI (MB) + StARS
#   3. Upgraded differential abundance
#      - ANCOM-BC2 + MaAsLin2 consensus
#   4. Random forest
#      - A/B/C three-class repeated stratified CV
#      - one-vs-rest ROC/AUC
#      - feature-selection stability
#   5. Spike-in detection limit
#      - detection-probability curve
#      - LOD50 / LOD95
#      - ROC/AUC and Youden threshold
#
# IMPORTANT SPIKE-IN NOTE
# Genus.xlsx, Genus_Count.xlsx and metadata.xlsx do NOT contain
# genome-specific spike-in trial results. Therefore the spike-in module
# creates a required-data template and runs only when real spike-in
# results are supplied. It never fabricates detection data.
#
# Nature/NPG-like colours
# A = #E64B35; B = #4DBBD5; C = #00A087
############################################################

rm(list = ls())
set.seed(20260812)

############################
# 0. File resolution
############################

resolve_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find any of: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

genus_file <- "Genus.xlsx"
count_file <- "Genus_Count.xlsx"
metadata_file <- "metadata.xlsx"

# Real spike-in trial result file, if available.
spikein_file <- "SpikeIn_results.csv"

outdir <- "Advanced_genus_analyses"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

############################
# 1. Packages
############################

cran_pkgs <- c(
  "readxl", "dplyr", "tidyr", "tibble", "ggplot2", "ggrepel",
  "vegan", "geosphere", "igraph", "ggraph", "scales",
  "ranger", "pROC", "patchwork", "pheatmap", "car"
)

need_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(need_cran) > 0) install.packages(need_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("SpiecEasi", "ANCOMBC", "phyloseq", "Maaslin2")
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
  library(vegan)
  library(geosphere)
  library(igraph)
  library(ggraph)
  library(scales)
  library(ranger)
  library(pROC)
  library(patchwork)
  library(pheatmap)
  library(car)
  library(SpiecEasi)
  library(ANCOMBC)
  library(phyloseq)
  library(Maaslin2)
})

############################
# 2. Plot style
############################

nature_clade <- c(
  "A" = "#E64B35",
  "B" = "#4DBBD5",
  "C" = "#00A087"
)

nature_misc <- c(
  "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
  "#DC0000", "#7E6148", "#B09C85"
)

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.45
      )
    )
}

save_fig <- function(p, filename, width = 6, height = 5) {
  ggsave(
    file.path(outdir, paste0(filename, ".pdf")),
    p, width = width, height = height, device = cairo_pdf
  )
  ggsave(
    file.path(outdir, paste0(filename, ".png")),
    p, width = width, height = height, dpi = 600
  )
}

############################
# 3. Read and align data
############################

genus0 <- read_excel(genus_file)
count0 <- read_excel(count_file)
meta0  <- read_excel(metadata_file)

stopifnot("SampleID" %in% names(genus0))
stopifnot("SampleID" %in% names(count0))
stopifnot("SampleID" %in% names(meta0))

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

genus0$SampleID <- as.character(genus0$SampleID)
count0$SampleID <- as.character(count0$SampleID)

common_samples <- Reduce(
  intersect,
  list(meta$SampleID, genus0$SampleID, count0$SampleID)
)

if (length(common_samples) != nrow(meta)) {
  warning(
    "Not all metadata samples are present in both abundance tables. ",
    "Analyses use common samples only."
  )
}

meta <- meta %>%
  filter(SampleID %in% common_samples)

genus0 <- genus0[
  match(meta$SampleID, genus0$SampleID),
]
count0 <- count0[
  match(meta$SampleID, count0$SampleID),
]

stopifnot(all(meta$SampleID == genus0$SampleID))
stopifnot(all(meta$SampleID == count0$SampleID))

genus <- genus0 %>%
  column_to_rownames("SampleID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

counts <- count0 %>%
  column_to_rownames("SampleID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# Harmonize taxa that occur in both tables.
common_taxa <- intersect(colnames(genus), colnames(counts))
genus <- genus[, common_taxa, drop = FALSE]
counts <- counts[, common_taxa, drop = FALSE]

# Remove all-zero features.
keep_nonzero <- colSums(counts, na.rm = TRUE) > 0
counts <- counts[, keep_nonzero, drop = FALSE]
genus  <- genus[, colnames(counts), drop = FALSE]

# QC checks.
if (any(counts < 0, na.rm = TRUE)) stop("Negative raw counts found.")
if (any(abs(as.matrix(counts) - round(as.matrix(counts))) > 1e-8,
        na.rm = TRUE)) {
  stop("Genus_Count input contains non-integer values.")
}

# If Genus table is already relative abundance, retain it.
# Otherwise normalize to relative abundance.
rs <- rowSums(genus, na.rm = TRUE)
if (!all(abs(rs - 1) < 0.02)) {
  genus <- genus / rowSums(genus)
}

# Scale environmental variables.
meta <- meta %>%
  mutate(
    Temperature_z = as.numeric(scale(Temperature)),
    Humidity_z    = as.numeric(scale(Humidity)),
    Altitude_z    = as.numeric(scale(Altitude)),
    Longitude_z   = as.numeric(scale(Longitude)),
    Latitude_z    = as.numeric(scale(Latitude))
  )

write.csv(
  data.frame(
    SampleID = meta$SampleID,
    GenusRelativeSum = rowSums(genus),
    GenusCountSum = rowSums(counts)
  ),
  file.path(outdir, "00_input_QC.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(table(meta$Clade, meta$Stage)),
  file.path(outdir, "00_Clade_by_Stage.csv")
)

############################################################
# 4. PROCRUSTES ANALYSIS
############################################################

# Goal:
# Compare the genus-level microbiome ordination with:
# (A) measured environment
# (B) geographic distances
#
# Microbiome: Bray-Curtis PCoA
# Environment: Euclidean PCA after z-scaling
# Geography: Haversine distance PCoA
#
# PROTEST provides permutation significance.

bray <- vegdist(genus, method = "bray")

micro_pcoa <- cmdscale(
  bray,
  eig = TRUE,
  k = 2,
  add = TRUE
)

micro_scores <- as.data.frame(micro_pcoa$points)
colnames(micro_scores) <- c("M1", "M2")
rownames(micro_scores) <- meta$SampleID

# Environmental PCA
env_mat <- meta %>%
  select(
    Temperature_z, Humidity_z, Altitude_z,
    Longitude_z, Latitude_z
  ) %>%
  as.data.frame()

rownames(env_mat) <- meta$SampleID

env_pca <- prcomp(
  env_mat,
  center = FALSE,
  scale. = FALSE
)

env_scores <- env_pca$x[, 1:2, drop = FALSE]

# Symmetric Procrustes + PROTEST
proc_env <- procrustes(
  micro_scores,
  env_scores,
  symmetric = TRUE
)

set.seed(20260812)
prot_env <- protest(
  micro_scores,
  env_scores,
  permutations = 9999,
  symmetric = TRUE
)

capture.output(
  prot_env,
  file = file.path(outdir, "01A_PROTEST_environment.txt")
)

# Extract fitted environmental configuration in microbiome coordinates.
env_rot <- fitted(proc_env, truemean = FALSE)

proc_env_df <- data.frame(
  SampleID = meta$SampleID,
  Micro1 = proc_env$X[, 1],
  Micro2 = proc_env$X[, 2],
  Env1 = env_rot[, 1],
  Env2 = env_rot[, 2],
  Clade = meta$Clade
)

write.csv(
  proc_env_df,
  file.path(outdir, "01A_Procrustes_environment_coordinates.csv"),
  row.names = FALSE
)

r_env <- prot_env$t0
p_env <- prot_env$signif

p_proc_env <- ggplot(proc_env_df) +
  geom_segment(
    aes(
      x = Env1, y = Env2,
      xend = Micro1, yend = Micro2,
      colour = Clade
    ),
    alpha = 0.45,
    linewidth = 0.45,
    arrow = arrow(length = unit(0.06, "inches"))
  ) +
  geom_point(
    aes(Env1, Env2),
    shape = 1,
    size = 2.4,
    colour = "grey45"
  ) +
  geom_point(
    aes(Micro1, Micro2, fill = Clade),
    shape = 21,
    size = 2.8,
    colour = "black",
    stroke = 0.25
  ) +
  scale_fill_manual(values = nature_clade) +
  scale_colour_manual(values = nature_clade) +
  labs(
    x = "Procrustes axis 1",
    y = "Procrustes axis 2",
    fill = "Clade",
    colour = "Clade",
    title = sprintf(
      "Microbiome vs environment: Procrustes r = %.3f, P = %.4g",
      r_env, p_env
    )
  ) +
  theme_nature()

save_fig(p_proc_env, "01A_Procrustes_environment", 6.5, 5.5)

# Residual plot: samples with largest discordance
proc_resid <- data.frame(
  SampleID = meta$SampleID,
  Residual = residuals(proc_env),
  Clade = meta$Clade
) %>%
  arrange(desc(Residual))

write.csv(
  proc_resid,
  file.path(outdir, "01B_Procrustes_environment_residuals.csv"),
  row.names = FALSE
)

p_proc_res <- proc_resid %>%
  slice_head(n = 20) %>%
  mutate(
    SampleID = factor(
      SampleID,
      levels = rev(SampleID)
    )
  ) %>%
  ggplot(
    aes(
      x = Residual,
      y = SampleID,
      fill = Clade
    )
  ) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = nature_clade) +
  labs(
    x = "Procrustes residual",
    y = NULL,
    fill = "Clade",
    title = "Samples with greatest microbiome-environment discordance"
  ) +
  theme_nature()

save_fig(p_proc_res, "01B_Procrustes_environment_residuals", 6.2, 5.5)

# Geographic Procrustes using Haversine distances.
coords <- as.matrix(meta[, c("Longitude", "Latitude")])
geo_dist_mat <- geosphere::distm(
  coords,
  fun = geosphere::distHaversine
)
geo_dist <- as.dist(geo_dist_mat)

geo_pcoa <- cmdscale(
  geo_dist,
  eig = TRUE,
  k = 2,
  add = TRUE
)

geo_scores <- geo_pcoa$points
rownames(geo_scores) <- meta$SampleID

proc_geo <- procrustes(
  micro_scores,
  geo_scores,
  symmetric = TRUE
)

set.seed(20260812)
prot_geo <- protest(
  micro_scores,
  geo_scores,
  permutations = 9999,
  symmetric = TRUE
)

capture.output(
  prot_geo,
  file = file.path(outdir, "01C_PROTEST_geography.txt")
)

geo_rot <- fitted(proc_geo, truemean = FALSE)

proc_geo_df <- data.frame(
  SampleID = meta$SampleID,
  Micro1 = proc_geo$X[, 1],
  Micro2 = proc_geo$X[, 2],
  Geo1 = geo_rot[, 1],
  Geo2 = geo_rot[, 2],
  Clade = meta$Clade
)

p_proc_geo <- ggplot(proc_geo_df) +
  geom_segment(
    aes(
      x = Geo1, y = Geo2,
      xend = Micro1, yend = Micro2,
      colour = Clade
    ),
    alpha = 0.45,
    linewidth = 0.45
  ) +
  geom_point(aes(Geo1, Geo2), shape = 1, size = 2.4, colour = "grey45") +
  geom_point(
    aes(Micro1, Micro2, fill = Clade),
    shape = 21, size = 2.8, colour = "black", stroke = 0.25
  ) +
  scale_fill_manual(values = nature_clade) +
  scale_colour_manual(values = nature_clade) +
  labs(
    x = "Procrustes axis 1",
    y = "Procrustes axis 2",
    fill = "Clade",
    colour = "Clade",
    title = sprintf(
      "Microbiome vs geography: Procrustes r = %.3f, P = %.4g",
      prot_geo$t0, prot_geo$signif
    )
  ) +
  theme_nature()

save_fig(p_proc_geo, "01C_Procrustes_geography", 6.5, 5.5)

KT_Residual <- kruskal.test(
  Residual ~ Clade,
  data = proc_resid
)

PWT_Residual <- pairwise.wilcox.test(
  proc_resid$Residual,
  proc_resid$Clade,
  p.adjust.method = "BH"
)


############################################################
# 5. MICROBIAL ASSOCIATION NETWORK: SPIEC-EASI
############################################################

# With n = 69 and ~1900 genera, an unfiltered network is not statistically
# defensible. We therefore filter by prevalence and abundance and cap the
# network to the most informative 100 genera.
#
# This is an association / conditional-independence network.
# Do NOT describe edges as proven biological interaction.

net_prev_cut <- 0.20
net_mean_rel_cut <- 1e-4   # 0.01%
net_max_taxa <- 100

net_prev <- colMeans(counts > 0)
net_mean <- colMeans(genus)

net_taxa <- names(which(
  net_prev >= net_prev_cut &
    net_mean >= net_mean_rel_cut
))

if (length(net_taxa) < 20) {
  # Relax abundance threshold, keep prevalence criterion.
  net_taxa <- names(which(net_prev >= net_prev_cut))
}

# Cap feature number by mean relative abundance.
if (length(net_taxa) > net_max_taxa) {
  net_taxa <- net_taxa[
    order(net_mean[net_taxa], decreasing = TRUE)[1:net_max_taxa]
  ]
}

net_counts <- counts[, net_taxa, drop = FALSE]
net_counts <- data.matrix(net_counts)

write.csv(
  data.frame(
    Genus = net_taxa,
    Prevalence = net_prev[net_taxa],
    MeanRelativeAbundance = net_mean[net_taxa]
  ),
  file.path(outdir, "02A_network_taxa_used.csv"),
  row.names = FALSE
)

# SPIEC-EASI MB model with StARS stability selection.
set.seed(20260812)

se_mb <- SpiecEasi::spiec.easi(
  net_counts,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 30,
  sel.criterion = "stars",
  pulsar.select = TRUE,
  pulsar.params = list(
    rep.num = 50,
    seed = 20260812
  )
)

saveRDS(
  se_mb,
  file.path(outdir, "02B_SPIECEASI_model.rds")
)

adj <- SpiecEasi::getRefit(se_mb)
beta <- SpiecEasi::symBeta(
  SpiecEasi::getOptBeta(se_mb),
  mode = "maxabs"
)

# Ensure dimnames.
rownames(adj) <- colnames(adj) <- colnames(net_counts)
rownames(beta) <- colnames(beta) <- colnames(net_counts)

# Extract edges.
idx <- which(
  upper.tri(adj) & as.matrix(adj) != 0,
  arr.ind = TRUE
)

if (nrow(idx) > 0) {
  edges <- data.frame(
    From = rownames(adj)[idx[, 1]],
    To = colnames(adj)[idx[, 2]],
    Weight = beta[idx],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Sign = ifelse(Weight >= 0, "Positive", "Negative"),
      AbsWeight = abs(Weight)
    )
} else {
  edges <- data.frame(
    From = character(),
    To = character(),
    Weight = numeric(),
    Sign = character(),
    AbsWeight = numeric()
  )
}

write.csv(
  edges,
  file.path(outdir, "02C_SPIECEASI_edges.csv"),
  row.names = FALSE
)

# Assign node to clade where mean relative abundance is highest.
clade_means <- genus[, net_taxa, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(meta %>% select(SampleID, Clade), by = "SampleID") %>%
  pivot_longer(
    -c(SampleID, Clade),
    names_to = "Genus",
    values_to = "Abundance"
  ) %>%
  group_by(Genus, Clade) %>%
  summarise(Mean = mean(Abundance), .groups = "drop") %>%
  group_by(Genus) %>%
  slice_max(Mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(EnrichedClade = Clade)

g <- igraph::graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = data.frame(name = net_taxa)
)

nodes <- data.frame(
  Genus = V(g)$name,
  Degree = igraph::degree(g),
  Betweenness = igraph::betweenness(g, normalized = TRUE),
  MeanRelativeAbundance = net_mean[V(g)$name]
) %>%
  left_join(clade_means, by = "Genus")

write.csv(
  nodes,
  file.path(outdir, "02D_SPIECEASI_node_metrics.csv"),
  row.names = FALSE
)

# Network plot.
if (ecount(g) > 0) {
  nodes_plot <- nodes[
    match(V(g)$name, nodes$Genus),
    ,
    drop = FALSE
  ]
  
  # Check matching
  stopifnot(
    all(
      V(g)$name == nodes_plot$Genus
    )
  )
  
  V(g)$Degree <-
    nodes_plot$Degree
  
  V(g)$Betweenness <-
    nodes_plot$Betweenness
  
  V(g)$MeanRelativeAbundance <-
    nodes_plot$MeanRelativeAbundance
  
  V(g)$EnrichedClade <-
    as.character(
      nodes_plot$EnrichedClade
    )
  E(g)$Sign <-
    edges$Sign
  
  E(g)$AbsWeight <-
    edges$AbsWeight
  
  hub_threshold <-
    quantile(
      V(g)$Degree,
      probs = 0.90,
      na.rm = TRUE
    )
  
  V(g)$HubLabel <-
    ifelse(
      V(g)$Degree >= hub_threshold,
      V(g)$name,
      ""
    )
  set.seed(20260812)
  vertex_attr_names(g)
  head(
    data.frame(
      name = V(g)$name,
      Degree = V(g)$Degree,
      EnrichedClade = V(g)$EnrichedClade
    )
  )

  # Degree plot.
  p_net <-
    ggraph(
      g,
      layout = "fr"
    ) +
    geom_edge_link(
      aes(
        edge_colour = Sign,
        edge_width = AbsWeight
      ),
      alpha = 0.45
    ) +
    geom_node_point(
      aes(
        size = Degree,
        fill = EnrichedClade
      ),
      shape = 21,
      colour = "black",
      stroke = 0.25
    ) +
    geom_node_text(
      aes(
        label = HubLabel
      ),
      repel = TRUE,
      size = 3
    ) +
    
    scale_fill_manual(
      values = nature_clade,
      na.value = "grey70"
    ) +
    
    scale_edge_colour_manual(
      values = c(
        "Positive" = "#3C5488",
        "Negative" = "#E64B35"
      )
    ) +
    
    scale_edge_width(
      range = c(
        0.25,
        1.4
      )
    ) +
    
    scale_size_continuous(
      range = c(
        2.5,
        8
      )
    ) +
    
    labs(
      fill = "Highest mean\nabundance",
      size = "Degree",
      edge_colour = "Association",
      edge_width = "|weight|",
      title = "SPIEC-EASI genus association network"
    ) +
    
    theme_void(
      base_size = 11
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold"
        ),
      legend.position =
        "right"
    )
  save_fig(
    p_net,
    "02E_SPIECEASI_network",
    8.5,
    7.2
  )
  p_degree <-
    nodes %>%
    
    arrange(
      desc(Degree)
    ) %>%
    
    slice_head(
      n = 20
    ) %>%
    
    mutate(
      Genus = factor(
        Genus,
        levels = rev(Genus)
      )
    ) %>%
    
    ggplot(
      aes(
        x = Degree,
        y = Genus,
        fill = EnrichedClade
      )
    ) +
    
    geom_col(
      width = 0.7
    ) +
    
    scale_fill_manual(
      values = nature_clade,
      na.value = "grey70"
    ) +
    
    labs(
      x = "Degree",
      y = NULL,
      fill = "Highest mean\nabundance in",
      title = "Highly connected microbial taxa"
    ) +
    
    theme_nature()

  save_fig(p_degree, "02F_SPIECEASI_top_hubs", 6.2, 5.5)
}

############################################################
# 6. UPGRADED DIFFERENTIAL ABUNDANCE
#    ANCOM-BC2 + MaAsLin2 consensus
############################################################

# We use counts for ANCOM-BC2 and relative abundance for MaAsLin2.
# Clade is the biological variable of interest.
# Stage + measured environmental variables are adjustment covariates.
#
# Because Clade A contains only adults, we DO NOT fit Clade*Stage.

# ---------- ANCOM-BC2 ----------
otu <- phyloseq::otu_table(
  t(as.matrix(counts)),
  taxa_are_rows = TRUE
)

sam <- phyloseq::sample_data(
  data.frame(meta, row.names = meta$SampleID)
)

ps <- phyloseq::phyloseq(otu, sam)

set.seed(20260812)

ancom <- ANCOMBC::ancombc2(
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
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE
)

saveRDS(ancom, file.path(outdir, "03A_ANCOMBC2_object.rds"))

ancom_res <- as.data.frame(ancom$res)
if (!"taxon" %in% names(ancom_res)) ancom_res$taxon <- rownames(ancom_res)

write.csv(
  ancom_res,
  file.path(outdir, "03A_ANCOMBC2_primary.csv"),
  row.names = FALSE
)

if (!is.null(ancom$res_global)) {
  global_res <- as.data.frame(ancom$res_global)
  if (!"taxon" %in% names(global_res)) global_res$taxon <- rownames(global_res)
  write.csv(
    global_res,
    file.path(outdir, "03B_ANCOMBC2_global.csv"),
    row.names = FALSE
  )
}

if (!is.null(ancom$res_pair)) {
  pair_res <- as.data.frame(ancom$res_pair)
  if (!"taxon" %in% names(pair_res)) pair_res$taxon <- rownames(pair_res)
  write.csv(
    pair_res,
    file.path(outdir, "03C_ANCOMBC2_pairwise.csv"),
    row.names = FALSE
  )
}

extract_ancom_contrast <- function(res, level, reference = "A") {
  tag <- paste0("Clade", level)
  lfc_col <- paste0("lfc_", tag)
  se_col <- paste0("se_", tag)
  p_col <- paste0("p_", tag)
  q_col <- paste0("q_", tag)
  diff_col <- paste0("diff_", tag)

  if (!lfc_col %in% names(res)) return(NULL)

  data.frame(
    Genus = res$taxon,
    Contrast = paste0(level, " vs ", reference),
    ANCOM_LFC = res[[lfc_col]],
    ANCOM_SE = if (se_col %in% names(res)) res[[se_col]] else NA,
    ANCOM_P = if (p_col %in% names(res)) res[[p_col]] else NA,
    ANCOM_Q = if (q_col %in% names(res)) res[[q_col]] else NA,
    ANCOM_sig = if (diff_col %in% names(res)) {
      as.logical(res[[diff_col]])
    } else {
      res[[q_col]] < 0.05
    }
  )
}

ancom_tidy <- bind_rows(
  extract_ancom_contrast(ancom_res, "B"),
  extract_ancom_contrast(ancom_res, "C")
)

write.csv(
  ancom_tidy,
  file.path(outdir, "03D_ANCOMBC2_tidy_Clade.csv"),
  row.names = FALSE
)

# ---------- MaAsLin2 ----------
prev <- colMeans(genus > 0)
maas_taxa <- names(prev[prev >= 0.10])
maas_data <- genus[, maas_taxa, drop = FALSE]

maas_meta <- meta %>%
  select(
    SampleID, Clade, Stage,
    Temperature_z, Humidity_z, Altitude_z,
    Latitude_z, Longitude_z
  ) %>%
  column_to_rownames("SampleID")

maas_dir <- file.path(outdir, "03E_MaAsLin2")
dir.create(maas_dir, recursive = TRUE, showWarnings = FALSE)

maas_fit <- Maaslin2::Maaslin2(
  input_data = maas_data,
  input_metadata = maas_meta,
  output = maas_dir,
  fixed_effects = c(
    "Clade", "Stage",
    "Temperature_z", "Humidity_z", "Altitude_z",
    "Latitude_z", "Longitude_z"
  ),
  random_effects = NULL,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  standardize = FALSE,
  min_abundance = 0,
  min_prevalence = 0,
  reference = c("Clade,A", "Stage,Adult"),
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

maas_res <- maas_fit$results

write.csv(
  maas_res,
  file.path(outdir, "03E_MaAsLin2_all_results.csv"),
  row.names = FALSE
)

maas_clade <- maas_res %>%
  filter(
    metadata == "Clade",
    value %in% c("B", "C")
  ) %>%
  transmute(
    Genus = feature,
    Contrast = paste0(value, " vs A"),
    Maaslin_coef = coef,
    Maaslin_SE = stderr,
    Maaslin_P = pval,
    Maaslin_Q = qval,
    Maaslin_sig = qval < 0.05
  )

# ---------- Consensus ----------
consensus <- full_join(
  ancom_tidy,
  maas_clade,
  by = c("Genus", "Contrast")
) %>%
  mutate(
    ANCOM_sig = replace_na(ANCOM_sig, FALSE),
    Maaslin_sig = replace_na(Maaslin_sig, FALSE),
    Methods = case_when(
      ANCOM_sig & Maaslin_sig ~ "ANCOM-BC2 + MaAsLin2",
      ANCOM_sig ~ "ANCOM-BC2 only",
      Maaslin_sig ~ "MaAsLin2 only",
      TRUE ~ "Not significant"
    ),
    Robust = ANCOM_sig & Maaslin_sig
  )

write.csv(
  consensus,
  file.path(outdir, "03F_DA_consensus.csv"),
  row.names = FALSE
)

# ANCOM-BC2 forest plot for significant/robust taxa.
forest <- consensus %>%
  filter(ANCOM_sig, is.finite(ANCOM_LFC)) %>%
  group_by(Contrast) %>%
  slice_max(abs(ANCOM_LFC), n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Label = paste0(Contrast, "__", Genus),
    Label = factor(Label, levels = rev(unique(Label)))
  )

if (nrow(forest) > 0) {
  p_da <- ggplot(
    forest,
    aes(ANCOM_LFC, Label, colour = Methods)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      colour = "grey45"
    ) +
    geom_errorbar(
      aes(
        xmin = ANCOM_LFC - 1.96 * ANCOM_SE,
        xmax = ANCOM_LFC + 1.96 * ANCOM_SE
      ),
      height = 0.15,
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    geom_point(size = 2.8) +
    facet_wrap(~Contrast, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
    scale_colour_manual(
      values = c(
        "ANCOM-BC2 + MaAsLin2" = "#E64B35",
        "ANCOM-BC2 only" = "#3C5488",
        "MaAsLin2 only" = "#00A087",
        "Not significant" = "grey75"
      )
    ) +
    labs(
      x = "ANCOM-BC2 log fold change",
      y = NULL,
      colour = "Evidence",
      title = "Differential genus abundance with cross-method validation"
    ) +
    theme_nature()

  save_fig(p_da, "03G_DA_consensus_forest", 9.2, 7.0)
}

# Methods-support summary.
support_summary <- consensus %>%
  filter(Methods != "Not significant") %>%
  count(Contrast, Methods)

p_support <- ggplot(
  support_summary,
  aes(Contrast, n, fill = Methods)
) +
  geom_col(position = "stack", width = 0.65) +
  scale_fill_manual(
    values = c(
      "ANCOM-BC2 + MaAsLin2" = "#E64B35",
      "ANCOM-BC2 only" = "#3C5488",
      "MaAsLin2 only" = "#00A087"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of significant genera",
    fill = "Evidence",
    title = "Cross-method support for differential genera"
  ) +
  theme_nature()

save_fig(p_support, "03H_DA_method_support", 6.2, 4.8)

############################################################
# 7. RANDOM FOREST: THREE-CLASS + FEATURE STABILITY
############################################################

# To reduce p >> n and unstable importance:
# prevalence >= 10%, then retain up to 300 most variable genera.
rf_prev <- colMeans(genus > 0)
rf_taxa <- names(rf_prev[rf_prev >= 0.10])

rf_x <- genus[, rf_taxa, drop = FALSE]

# log transform; RF does not require scaling.
rf_x <- log10(rf_x + 1e-8)

if (ncol(rf_x) > 300) {
  vars <- apply(rf_x, 2, var)
  rf_x <- rf_x[
    ,
    names(sort(vars, decreasing = TRUE))[1:300],
    drop = FALSE
  ]
}

rf_dat <- data.frame(
  Clade = meta$Clade,
  rf_x,
  check.names = FALSE
)

# Manual stratified fold generator.
make_stratified_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  folds <- vector("list", k)
  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    split_idx <- split(
      idx,
      rep(1:k, length.out = length(idx))
    )
    for (i in seq_len(k)) {
      folds[[i]] <- c(folds[[i]], split_idx[[i]])
    }
  }
  folds
}

rf_repeats <- 20
rf_k <- 5
num_trees <- 1000
top_n_features <- 20

pred_store <- list()
stability_store <- list()
model_counter <- 1

bad_names <- colnames(rf_dat)[
  make.names(
    colnames(rf_dat)
  ) != colnames(rf_dat)
]

length(bad_names)


for (r in seq_len(rf_repeats)) {

  folds <- make_stratified_folds(
    rf_dat$Clade,
    k = rf_k,
    seed = 20260812 + r
  )

  for (f in seq_len(rf_k)) {

    test_idx <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(rf_dat)), test_idx)

    train <- rf_dat[train_idx, , drop = FALSE]
    test  <- rf_dat[test_idx, , drop = FALSE]

    p <- ncol(train) - 1

    fit <- ranger::ranger(
      dependent.variable.name = "Clade",
      data = train,
      probability = TRUE,
      num.trees = num_trees,
      mtry = max(1, floor(sqrt(p))),
      min.node.size = 3,
      importance = "permutation",
      respect.unordered.factors = "order",
      seed = 20260812 + model_counter
    )

    pr <- predict(fit, data = test)$predictions

    tmp_pred <- data.frame(
      SampleID = meta$SampleID[test_idx],
      Truth = as.character(test$Clade),
      Repeat = r,
      Fold = f,
      pr,
      check.names = FALSE
    )

    pred_store[[model_counter]] <- tmp_pred

    imp <- sort(
      fit$variable.importance,
      decreasing = TRUE
    )

    topf <- head(names(imp), top_n_features)

    stability_store[[model_counter]] <- data.frame(
      Genus = topf,
      Importance = imp[topf],
      Repeat = r,
      Fold = f
    )

    model_counter <- model_counter + 1
  }
}

rf_pred_all <- bind_rows(pred_store)
rf_stability_raw <- bind_rows(stability_store)

write.csv(
  rf_pred_all,
  file.path(outdir, "04A_RF_repeatedCV_predictions.csv"),
  row.names = FALSE
)

# Aggregate repeated out-of-fold probabilities per sample.
prob_cols <- levels(meta$Clade)

rf_pred <- rf_pred_all %>%
  group_by(SampleID, Truth) %>%
  summarise(
    across(all_of(prob_cols), mean),
    .groups = "drop"
  )

rf_pred$Predicted <- prob_cols[
  max.col(
    as.matrix(rf_pred[, prob_cols]),
    ties.method = "first"
  )
]

conf_mat <- table(
  Truth = rf_pred$Truth,
  Predicted = rf_pred$Predicted
)

write.csv(
  as.data.frame.matrix(conf_mat),
  file.path(outdir, "04B_RF_confusion_matrix.csv")
)

accuracy <- mean(rf_pred$Truth == rf_pred$Predicted)

# One-vs-rest ROC.
roc_list <- list()
roc_curve <- list()
auc_tab <- list()

for (cl in prob_cols) {

  truth_bin <- rf_pred$Truth == cl

  roc_obj <- pROC::roc(
    response = truth_bin,
    predictor = rf_pred[[cl]],
    levels = c(FALSE, TRUE),
    direction = "<",
    quiet = TRUE
  )

  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- as.numeric(pROC::ci.auc(roc_obj))

  auc_tab[[cl]] <- data.frame(
    Clade = cl,
    AUC = auc_val,
    CI_low = ci_val[1],
    CI_high = ci_val[3]
  )

  roc_curve[[cl]] <- data.frame(
    Clade = cl,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )

  roc_list[[cl]] <- roc_obj
}

auc_df <- bind_rows(auc_tab)
macro_auc <- mean(auc_df$AUC)

write.csv(
  auc_df,
  file.path(outdir, "04C_RF_AUC_by_Clade.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste0("Repeated CV accuracy = ", round(accuracy, 4)),
    paste0("Macro one-vs-rest AUC = ", round(macro_auc, 4))
  ),
  file.path(outdir, "04C_RF_overall_performance.txt")
)

roc_df <- bind_rows(roc_curve)

p_rf_roc <- ggplot(
  roc_df,
  aes(FPR, TPR, colour = Clade)
) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = 2, colour = "grey60"
  ) +
  geom_line(linewidth = 1.05) +
  scale_colour_manual(values = nature_clade) +
  coord_equal() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    colour = "Clade",
    title = sprintf(
      "Three-class random forest: accuracy %.2f, macro-AUC %.2f",
      accuracy, macro_auc
    )
  ) +
  theme_nature()

save_fig(p_rf_roc, "04D_RF_multiclass_ROC", 5.8, 5.4)

# Feature stability:
# frequency with which each genus appears among top 20 across CV models.
n_models <- length(stability_store)

rf_stability <- rf_stability_raw %>%
  group_by(Genus) %>%
  summarise(
    SelectionCount = n(),
    SelectionFrequency = n() / n_models,
    MeanImportance = mean(Importance, na.rm = TRUE),
    MedianImportance = median(Importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(SelectionFrequency), desc(MeanImportance))

write.csv(
  rf_stability,
  file.path(outdir, "04E_RF_feature_stability.csv"),
  row.names = FALSE
)

p_stability <- rf_stability %>%
  slice_head(n = 30) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus))) %>%
  ggplot(
    aes(SelectionFrequency, Genus, size = MeanImportance)
  ) +
  geom_point(
    fill = "#3C5488",
    colour = "black",
    shape = 21,
    alpha = 0.85
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    x = "Top-20 selection frequency across CV models",
    y = NULL,
    size = "Mean\nimportance",
    title = "Random-forest feature stability"
  ) +
  theme_nature()

save_fig(p_stability, "04F_RF_feature_stability", 7.0, 7.0)

# Final RF model trained on all data for final importance ranking.
rf_final <- ranger::ranger(
  dependent.variable.name = "Clade",
  data = rf_dat,
  probability = TRUE,
  num.trees = 2000,
  mtry = max(1, floor(sqrt(ncol(rf_dat) - 1))),
  min.node.size = 3,
  importance = "permutation",
  respect.unordered.factors = "order",
  seed = 20260812
)

final_imp <- data.frame(
  Genus = names(rf_final$variable.importance),
  Importance = as.numeric(rf_final$variable.importance)
) %>%
  left_join(
    rf_stability %>%
      select(Genus, SelectionFrequency),
    by = "Genus"
  ) %>%
  arrange(desc(Importance))

write.csv(
  final_imp,
  file.path(outdir, "04G_RF_final_importance.csv"),
  row.names = FALSE
)

############################################################
# 8. SPIKE-IN DETECTION LIMIT + ROC
############################################################

# REQUIRED real-data columns:
#   Isolate
#   BackgroundSample
#   SpikeFraction
#       e.g. 0, 0.000001, 0.00001, 0.0001, 0.001
#       NOTE: use fractions, so 0.001% = 0.00001
#   DetectedReads
#       genome-specific mapped/read count or another continuous detection score
#
# Optional:
#   Detected  (0/1)
#
# To estimate a reliable ROC, include zero-spike negative controls.
# To estimate a detection-probability curve / LOD95, use replicate
# backgrounds at each abundance level whenever possible.

if (!file.exists(spikein_file)) {

  spike_template <- data.frame(
    Isolate = character(),
    BackgroundSample = character(),
    SpikeFraction = numeric(),
    DetectedReads = numeric(),
    Detected = integer()
  )

  write.csv(
    spike_template,
    file.path(outdir, "05A_SpikeIn_results_TEMPLATE.csv"),
    row.names = FALSE
  )

  writeLines(
    c(
      "Spike-in ROC was not run because real genome-specific spike-in trial results were not supplied.",
      "Copy 05A_SpikeIn_results_TEMPLATE.csv to SpikeIn_results.csv and fill it with real results.",
      "Include zero-spike negative controls for ROC specificity.",
      "Use replicate metagenomic backgrounds at each spike abundance for detection-probability / LOD estimation.",
      "Do not infer ROC from Genus.xlsx or Genus_Count.xlsx because they do not contain spike-in truth labels."
    ),
    file.path(outdir, "05A_SpikeIn_README.txt")
  )

} else {

  spike <- read.csv(
    spikein_file,
    stringsAsFactors = FALSE
  )

  required <- c(
    "Isolate",
    "BackgroundSample",
    "SpikeFraction",
    "DetectedReads"
  )

  if (!all(required %in% names(spike))) {
    stop(
      "SpikeIn_results.csv must contain: ",
      paste(required, collapse = ", ")
    )
  }

  if (!"Detected" %in% names(spike)) {
    # Detection definition can be changed if laboratory/bioinformatic
    # validation uses another minimum-read criterion.
    min_detect_reads <- 1
    spike$Detected <- as.integer(
      spike$DetectedReads >= min_detect_reads
    )
  }

  spike <- spike %>%
    mutate(
      Truth = SpikeFraction > 0,
      Log10Spike = log10(
        pmax(SpikeFraction, 1e-10)
      )
    )

  write.csv(
    spike,
    file.path(outdir, "05B_SpikeIn_cleaned.csv"),
    row.names = FALSE
  )

  # Detection probability curves and LOD50 / LOD95.
  lod_results <- list()
  detection_curve <- list()

  for (iso in unique(spike$Isolate)) {

    dat <- spike %>%
      filter(Isolate == iso)

    if (
      length(unique(dat$Detected)) < 2 ||
      length(unique(dat$SpikeFraction)) < 3
    ) {
      next
    }

    fit <- glm(
      Detected ~ Log10Spike,
      family = binomial(),
      data = dat
    )

    grid <- data.frame(
      Log10Spike = seq(
        min(dat$Log10Spike),
        max(dat$Log10Spike),
        length.out = 300
      )
    )

    pr <- predict(
      fit,
      newdata = grid,
      type = "response"
    )

    detection_curve[[iso]] <- data.frame(
      Isolate = iso,
      Log10Spike = grid$Log10Spike,
      SpikeFraction = 10^grid$Log10Spike,
      DetectionProbability = pr
    )

    b <- coef(fit)

    lod50_log <- (qlogis(0.50) - b[1]) / b[2]
    lod95_log <- (qlogis(0.95) - b[1]) / b[2]

    lod_results[[iso]] <- data.frame(
      Isolate = iso,
      LOD50_fraction = 10^lod50_log,
      LOD95_fraction = 10^lod95_log,
      LOD50_percent = 100 * 10^lod50_log,
      LOD95_percent = 100 * 10^lod95_log
    )
  }

  lod_df <- bind_rows(lod_results)
  curve_df <- bind_rows(detection_curve)

  write.csv(
    lod_df,
    file.path(outdir, "05C_SpikeIn_LOD50_LOD95.csv"),
    row.names = FALSE
  )

  if (nrow(curve_df) > 0) {

    obs_summary <- spike %>%
      group_by(Isolate, SpikeFraction) %>%
      summarise(
        DetectionRate = mean(Detected),
        N = n(),
        .groups = "drop"
      )

    p_det <- ggplot() +
      geom_line(
        data = curve_df,
        aes(
          x = 100 * SpikeFraction,
          y = DetectionProbability,
          colour = Isolate
        ),
        linewidth = 1.0
      ) +
      geom_point(
        data = obs_summary,
        aes(
          x = 100 * SpikeFraction,
          y = DetectionRate,
          colour = Isolate,
          size = N
        ),
        alpha = 0.85
      ) +
      scale_x_log10() +
      scale_colour_manual(
        values = setNames(
          nature_misc[
            seq_len(min(length(unique(spike$Isolate)), length(nature_misc)))
          ],
          unique(spike$Isolate)[
            seq_len(min(length(unique(spike$Isolate)), length(nature_misc)))
          ]
        )
      ) +
      labs(
        x = "Spike-in relative abundance (%)",
        y = "Detection probability",
        colour = "Isolate",
        size = "Replicates",
        title = "Spike-in detection probability curves"
      ) +
      theme_nature()

    save_fig(p_det, "05D_SpikeIn_detection_probability", 7.0, 5.2)
  }

  # ROC based on continuous DetectedReads score:
  # Truth = real spike present (>0) vs zero-spike control.
  roc_results <- list()
  roc_curves <- list()

  for (iso in unique(spike$Isolate)) {

    dat <- spike %>%
      filter(Isolate == iso)

    if (length(unique(dat$Truth)) < 2) next

    ro <- pROC::roc(
      response = dat$Truth,
      predictor = dat$DetectedReads,
      levels = c(FALSE, TRUE),
      direction = "<",
      quiet = TRUE
    )

    auc_val <- as.numeric(pROC::auc(ro))
    ci_val <- as.numeric(pROC::ci.auc(ro))

    best <- pROC::coords(
      ro,
      x = "best",
      best.method = "youden",
      ret = c(
        "threshold",
        "sensitivity",
        "specificity"
      ),
      transpose = FALSE
    )

    roc_results[[iso]] <- data.frame(
      Isolate = iso,
      AUC = auc_val,
      CI_low = ci_val[1],
      CI_high = ci_val[3],
      YoudenThreshold = best$threshold,
      Sensitivity = best$sensitivity,
      Specificity = best$specificity
    )

    roc_curves[[iso]] <- data.frame(
      Isolate = iso,
      FPR = 1 - ro$specificities,
      TPR = ro$sensitivities
    )
  }

  spike_auc <- bind_rows(roc_results)
  spike_roc_df <- bind_rows(roc_curves)

  write.csv(
    spike_auc,
    file.path(outdir, "05E_SpikeIn_ROC_AUC.csv"),
    row.names = FALSE
  )

  if (nrow(spike_roc_df) > 0) {

    p_spike_roc <- ggplot(
      spike_roc_df,
      aes(FPR, TPR, colour = Isolate)
    ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2,
        colour = "grey60"
      ) +
      geom_line(linewidth = 1.05) +
      coord_equal() +
      labs(
        x = "1 - Specificity",
        y = "Sensitivity",
        colour = "Isolate",
        title = "Spike-in detection ROC curves"
      ) +
      theme_nature()

    save_fig(p_spike_roc, "05F_SpikeIn_ROC", 5.8, 5.4)
  }
}

############################################################
# 9. Final reproducibility files
############################################################

capture.output(
  sessionInfo(),
  file = file.path(outdir, "99_sessionInfo.txt")
)

writeLines(
  c(
    "PRIMARY INTERPRETATION GUIDE",
    "",
    "Procrustes:",
    "  Report PROTEST r and permutation P. A significant fit indicates concordance between ordination configurations, not causation.",
    "",
    "Network:",
    "  SPIEC-EASI edges are conditional associations, not experimentally proven interactions.",
    "  The global network is preferred because Clade A has only 9 samples, which is too small for a stable clade-specific high-dimensional network.",
    "",
    "Differential abundance:",
    "  Prioritize genera significant in ANCOM-BC2 and supported by MaAsLin2.",
    "  ANCOM-BC2 uses raw counts; MaAsLin2 uses relative abundance with multivariable adjustment.",
    "",
    "Random forest:",
    "  Report repeated-CV accuracy, class-specific AUCs, macro-AUC and top-20 selection frequency.",
    "  Stable biomarkers should recur across many CV models rather than appear in one full-data forest only.",
    "",
    "Spike-in:",
    "  ROC requires real zero-spike controls plus continuous detection scores.",
    "  LOD95 is most defensible when each spike level is tested across replicate metagenomic backgrounds."
  ),
  file.path(outdir, "00_ANALYSIS_INTERPRETATION_GUIDE.txt")
)

message("All available analyses completed.")
message("Output directory: ", normalizePath(outdir))
