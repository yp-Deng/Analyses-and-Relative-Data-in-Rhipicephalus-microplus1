#############################################################
# Panaroo pangenome analysis
#
# 1. Complete gene presence/absence heatmap
# 2. Pangenome accumulation curve
# 3. Core genome curve
#
# Three independent genera
#
#############################################################


library(tidyverse)
library(pheatmap)
library(vegan)
library(ggplot2)



###############################
# paths
###############################


PAN_DIR <- "three_genus_results"


OUTDIR <- "VF_BGC_Figures_V5/Panaroo_Pangenome_Figures"


dir.create(
  OUTDIR,
  showWarnings = FALSE
)



###############################
# genus list
###############################


genus_list <- c(
  "Genus1",
  "Genus2",
  "Genus3"
)



###############################
# Function:
# read Panaroo matrix
###############################


read_panaroo_matrix <- function(genus){
  
  
  file <- paste0(
    "three_genus_results/",
    genus,
    "/panaroo/gene_presence_absence.csv"
  )
  
  
  cat(
    "Reading:",
    file,
    "\n"
  )
  
  
  
  tab <- read.csv(
    file,
    check.names=FALSE,
    stringsAsFactors=FALSE
  )
  
  
  
  ############################
  # 自动识别菌株列
  ############################
  
  
  ignore_columns <- c(
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "Description",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Fragment ID",
    "QC"
  )
  
  
  
  genomes <- setdiff(
    colnames(tab),
    ignore_columns
  )
  
  
  
  cat(
    "Detected genomes:",
    length(genomes),
    "\n"
  )
  
  
  print(genomes)
  
  
  
  ############################
  # build matrix
  ############################
  
  
  mat <- matrix(
    0,
    nrow=nrow(tab),
    ncol=length(genomes)
  )
  
  
  
  rownames(mat)<-tab$Gene
  
  
  
  colnames(mat)<-genomes
  
  
  
  for(i in seq_along(genomes)){
    
    
    mat[,i] <-
      
      ifelse(
        is.na(tab[,genomes[i]]) |
          tab[,genomes[i]]=="",
        0,
        1
      )
    
    
  }
  
  
  
  mode(mat)<-"numeric"
  
  
  
  cat(
    "Matrix:",
    nrow(mat),
    "x",
    ncol(mat),
    "\n"
  )
  
  
  return(mat)
  
}

Genus_list <- c(
  "Genus1",
  "Genus2",
  "Genus3"
)


for(g in Genus_list){
  
  
  mat <- read_panaroo_matrix(g)
  
  
  write.csv(
    
    mat,
    
    paste0(
      g,
      "_gene_presence_absence_matrix.csv"
    )
    
  )
  
}



###############################
# 1. Gene presence absence heatmap
###############################
plot_heatmap <- function(mat, genus){
  
  
  cat(
    "Generating:",
    genus,
    "\n"
  )
  
  
  mat <- mat[
    rowSums(mat)>0,
    ,
    drop=FALSE
  ]
  
  
  cat(
    "Remaining genes:",
    nrow(mat),
    "\n"
  )
  
  
  
  if(nrow(mat)>2000){
    
    sdv <- apply(
      mat,
      1,
      sd
    )
    
    
    keep <- order(
      sdv,
      decreasing=TRUE
    )[1:2000]
    
    
    mat <- mat[
      keep,
      ,
      drop=FALSE
    ]
    
    
  }
  
  
  
  outfile <- paste0(
    genus,
    "_gene_presence_absence_heatmap.pdf"
  )
  
  
  pdf(
    outfile,
    width=10,
    height=12
  )
  
  
  
  heatmap(
    
    mat,
    
    Rowv=NA,
    
    Colv=NA,
    
    scale="none",
    
    col=c(
      "white",
      "#7B3294"
    ),
    
    labRow=FALSE,
    
    main=paste(
      genus,
      "gene presence/absence"
    )
    
  )
  
  
  
  dev.off()
  
  
  
  cat(
    "Finished:",
    outfile,
    "\n"
  )
  
  
}

for(g in c(
  "Genus1",
  "Genus2",
  "Genus3"
)){
  
  
  mat <- read_panaroo_matrix(g)
  
  
  plot_heatmap(
    mat,
    g
  )
  
  
}
###############################
# 2. Pangenome accumulation
###############################


pangenome_curve <- function(mat, genus){
  
  
  n <- ncol(mat)
  
  
  
  result <- data.frame()
  
  
  
  set.seed(123)
  
  
  
  # repeat random genome order
  
  
  for(i in 1:100){
    
    
    order <- sample(
      1:n
    )
    
    
    
    pan <- c()
    
    core <- c()
    
    
    
    for(j in 1:n){
      
      
      sub <- mat[,order[1:j],drop=FALSE]
      
      
      
      pan[j] <-
        
        sum(
          rowSums(sub)>0
        )
      
      
      
      core[j] <-
        
        sum(
          rowSums(sub)==j
        )
      
      
      
    }
    
    
    
    result <-
      
      rbind(
        
        result,
        
        data.frame(
          
          iteration=i,
          
          genome_number=1:n,
          
          pangenome=pan,
          
          core=core
          
        )
        
      )
    
    
  }
  
  
  
  summary <- result %>%
    
    group_by(
      genome_number
    )%>%
    
    summarise(
      
      pangenome_mean=
        mean(pangenome),
      
      pangenome_sd=
        sd(pangenome),
      
      core_mean=
        mean(core),
      
      core_sd=
        sd(core)
      
    )
  
  
  
  write.csv(
    
    summary,
    
    paste0(
      
      OUTDIR,
      
      "/",
      
      genus,
      
      "_pangenome_curve.csv"
      
    ),
    
    row.names=FALSE
    
  )
  
  
  
  #############################
  # plot
  #############################
  
  
  
  p1 <- ggplot(
    
    summary,
    
    aes(
      genome_number,
      pangenome_mean
    )
    
  )+
    
    geom_line(
      linewidth=1
    )+
    
    geom_ribbon(
      
      aes(
        
        ymin=pangenome_mean-pangenome_sd,
        
        ymax=pangenome_mean+pangenome_sd
        
      ),
      
      alpha=0.2
      
    )+
    
    
    theme_classic()+
    
    labs(
      
      x="Number of genomes",
      
      y="Pan-genome size",
      
      title=paste(
        genus,
        "pangenome accumulation"
      )
      
    )
  
  
  
  p2 <- ggplot(
    
    summary,
    
    aes(
      genome_number,
      core_mean
    )
    
  )+
    
    geom_line(
      
      linewidth=1
      
    )+
    
    
    geom_ribbon(
      
      aes(
        
        ymin=core_mean-core_sd,
        
        ymax=core_mean+core_sd
        
      ),
      
      alpha=0.2
      
    )+
    
    
    theme_classic()+
    
    
    labs(
      
      x="Number of genomes",
      
      y="Core genome size",
      
      title=paste(
        genus,
        "core genome decay"
      )
      
    )
  
  
  
  ggsave(
    
    paste0(
      
      OUTDIR,
      
      "/",
      
      genus,
      
      "_pangenome_accumulation.pdf"
      
    ),
    
    p1,
    
    width=7,
    
    height=5
    
  )
  
  
  
  ggsave(
    
    paste0(
      
      OUTDIR,
      
      "/",
      
      genus,
      
      "_core_genome_curve.pdf"
      
    ),
    
    p2,
    
    width=7,
    
    height=5
    
  )
  
  
  
}



###############################
# RUN
###############################


for(g in genus_list){
  
  
  cat(
    "====================\n"
  )
  
  cat(
    g,
    "\n"
  )
  
  
  mat <- read_panaroo_matrix(g)
  
  
  
  write.csv(
    
    mat,
    
    paste0(
      
      OUTDIR,
      
      "/",
      
      g,
      
      "_gene_presence_absence_matrix.csv"
      
    )
    
  )
  
  
  
  plot_heatmap(
    mat,
    g
  )
  
  
  
  pangenome_curve(
    mat,
    g
  )
  
  
  
}



cat(
  "
Finished!
"
)
