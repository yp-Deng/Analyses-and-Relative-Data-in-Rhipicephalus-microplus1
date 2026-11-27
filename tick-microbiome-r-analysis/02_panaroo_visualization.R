############################################################
# Panaroo downstream visualization pipeline
#
# Input:
# three_genus_results/
#
# Output:
# Pan-genome composition
# Gene presence/absence PCoA
# Core genome phylogeny
#
############################################################


############################################################
# Packages
############################################################

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "vegan",
  "ape",
  "ggtree",
  "treeio"
)


for(p in packages){
  
  if(!require(p, character.only = TRUE)){
    install.packages(p)
  }
  
}


############################################################
# Directory
############################################################


result_dir <- "three_genus_results"

outdir <- "panaroo_visualization"

dir.create(
  outdir,
  showWarnings = FALSE
)



############################################################
# 1. Pan-genome statistics
############################################################


get_statistics <- function(file, genus){
  
  
  x <- readLines(file)
  
  
  core <- as.numeric(
    gsub(
      ".*\\s",
      "",
      grep("^Core genes",x,value=T)
    )
  )
  
  
  soft <- as.numeric(
    gsub(
      ".*\\s",
      "",
      grep("^Soft core genes",x,value=T)
    )
  )
  
  
  shell <- as.numeric(
    gsub(
      ".*\\s",
      "",
      grep("^Shell genes",x,value=T)
    )
  )
  
  
  cloud <- as.numeric(
    gsub(
      ".*\\s",
      "",
      grep("^Cloud genes",x,value=T)
    )
  )
  
  
  total <- as.numeric(
    gsub(
      ".*\\s",
      "",
      grep("^Total genes",x,value=T)
    )
  )
  
  
  data.frame(
    Genus=genus,
    Core=core,
    Soft_core=soft,
    Shell=shell,
    Cloud=cloud,
    Total=total
  )
  
}



stats <- bind_rows(
  
  get_statistics(
    paste0(result_dir,"/Genus1/panaroo/summary_statistics.txt"),
    "Genus1"
  ),
  
  
  get_statistics(
    paste0(result_dir,"/Genus2/panaroo/summary_statistics.txt"),
    "Genus2"
  ),
  
  
  get_statistics(
    paste0(result_dir,"/Genus3/panaroo/summary_statistics.txt"),
    "Genus3"
  )
  
)



write.csv(
  stats,
  paste0(outdir,"/pan_genome_statistics.csv"),
  row.names=FALSE
)



############################################################
# Plot pan-genome composition
############################################################


plot_data <- stats %>%
  select(
    Genus,
    Core,
    Soft_core,
    Shell,
    Cloud
  )%>%
  pivot_longer(
    cols=-Genus,
    names_to="Category",
    values_to="Genes"
  )



p1 <- ggplot(
  plot_data,
  aes(
    x=Genus,
    y=Genes,
    fill=Category
  )
)+
  geom_bar(
    stat="identity"
  )+
  theme_bw()+
  ylab(
    "Number of gene clusters"
  )+
  theme(
    text=element_text(size=14)
  )



ggsave(
  paste0(outdir,"/Pan_genome_composition.pdf"),
  p1,
  width=8,
  height=6
)



############################################################
# 2. Gene presence absence PCoA
############################################################


run_pcoa <- function(csvfile, genus){
  
  
  gene <- read.csv(
    csvfile,
    check.names=FALSE
  )
  
  
  # Panaroo annotation columns
  mat <- gene[,15:ncol(gene)]
  
  
  rownames(mat)<-gene[,1]
  
  
  mat <- as.matrix(mat)
  
  
  mat[is.na(mat)]<-0
  
  
  # presence absence
  mat[mat>0]<-1
  
  
  
  # strain x gene
  
  mat <- t(mat)
  
  
  
  dist <- vegdist(
    mat,
    method="jaccard"
  )
  
  
  
  pcoa <- cmdscale(
    dist,
    k=2
  )
  
  
  df <- data.frame(
    Sample=rownames(pcoa),
    PC1=pcoa[,1],
    PC2=pcoa[,2],
    Genus=genus
  )
  
  
  
  p <- ggplot(
    df,
    aes(
      PC1,
      PC2,
      label=Sample
    )
  )+
    geom_point(
      size=4
    )+
    geom_text(
      size=3,
      vjust=-0.5
    )+
    theme_bw()+
    ggtitle(genus)
  
  
  
  ggsave(
    paste0(
      outdir,
      "/",
      genus,
      "_gene_presence_absence_PCoA.pdf"
    ),
    p,
    width=8,
    height=6
  )
  
  
  
  write.csv(
    df,
    paste0(
      outdir,
      "/",
      genus,
      "_PCoA_coordinates.csv"
    ),
    row.names=FALSE
  )
  
  
}



run_pcoa(
  paste0(result_dir,
         "/Genus1/panaroo/gene_presence_absence.csv"),
  "Genus1"
)


run_pcoa(
  paste0(result_dir,
         "/Genus2/panaroo/gene_presence_absence.csv"),
  "Genus2"
)


run_pcoa(
  paste0(result_dir,
         "/Genus3/panaroo/gene_presence_absence.csv"),
  "Genus3"
)




############################################################
# 3. Core genome phylogenetic trees
############################################################


plot_tree <- function(treefile, genus){
  
  
  tree <- read.tree(treefile)
  
  
  
  p <- ggtree(
    tree
  )+
    geom_tiplab(
      size=3
    )+
    theme_tree2()
  
  
  
  ggsave(
    paste0(
      outdir,
      "/",
      genus,
      "_core_genome_tree.pdf"
    ),
    p,
    width=8,
    height=10
  )
  
  
  
}




plot_tree(
  paste0(result_dir,
         "/Genus1/phylogeny/Genus1.treefile"),
  "Genus1"
)



plot_tree(
  paste0(result_dir,
         "/Genus2/phylogeny/Genus2.treefile"),
  "Genus2"
)



plot_tree(
  paste0(result_dir,
         "/Genus3/phylogeny/Genus3.treefile"),
  "Genus3"
)



############################################################
# Finish
############################################################


cat(
  "
====================================
Panaroo visualization finished
====================================

Output directory:

panaroo_visualization/

Files:

1. Pan_genome_composition.pdf

2. Genus1_gene_presence_absence_PCoA.pdf

3. Genus2_gene_presence_absence_PCoA.pdf

4. Genus3_gene_presence_absence_PCoA.pdf

5. Genus1_core_genome_tree.pdf

6. Genus2_core_genome_tree.pdf

7. Genus3_core_genome_tree.pdf

====================================

"
)