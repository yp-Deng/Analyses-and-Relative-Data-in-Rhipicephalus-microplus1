#############################################################
# Nature Communications Main Figure V4
#
# Tree-based genome classification
#
# Phylogeny
# VF repertoire
# BGC repertoire
# Panaroo core-shell-cloud
#
#############################################################


library(tidyverse)
library(ape)
library(patchwork)



############################
# PATH
############################


TREE_DIR <-
  "three_genus_results"


PAN_DIR <-
  "VF_BGC_Figures_V5/iTOL_VF_BGC_Panaroo"


VF_FILE <-
  "VF_BGC_results_V4/VF_summary/VF_gene_matrix_strain.csv"


BGC_FILE <-
  "VF_BGC_Figures_V5/tables/BGC_type_matrix.csv"



OUTDIR <-
  "MainFigure_V4"



dir.create(
  OUTDIR,
  recursive=TRUE,
  showWarnings=FALSE
)



############################
# 1. Extract genome groups
# from trees
############################


get_tree_tips <- function(genus){
  
  
  treefile <- list.files(
    
    paste0(
      TREE_DIR,
      "/",
      genus
    ),
    
    pattern="treefile$",
    
    recursive=TRUE,
    
    full.names=TRUE
    
  )
  
  
  if(length(treefile)==0){
    
    stop(
      paste(
        "No tree found:",
        genus
      )
    )
    
  }
  
  
  tree <- read.tree(
    treefile[1]
  )
  
  
  tips <- tree$tip.label
  
  
  data.frame(
    
    Genome=tips,
    
    Genus=genus
    
  )
  
}



genome_group <- bind_rows(
  
  get_tree_tips("Genus1"),
  
  get_tree_tips("Genus2"),
  
  get_tree_tips("Genus3")
  
)



write.csv(
  
  genome_group,
  
  paste0(
    OUTDIR,
    "/Genome_group_from_tree.csv"
  ),
  
  row.names=FALSE
  
)



cat(
  "Tree genome classification finished\n"
)



############################
# 2. Read Panaroo
############################


read_pan_itol <- function(genus){
  
  
  file <- paste0(
    
    PAN_DIR,
    
    "/",
    
    genus,
    
    "/Panaroo_core_shell_cloud_fixed.txt"
    
  )
  
  
  cat(
    "Reading:",
    file,
    "\n"
  )
  
  
  
  txt <- readLines(file)
  
  
  
  ############################
  # 找真正 DATA 行
  ############################
  
  
  start <- grep(
    "^DATA$",
    txt
  )
  
  
  
  if(length(start)==0){
    
    stop(
      paste(
        "No DATA section in",
        file
      )
    )
    
  }
  
  
  start <- start[1]+1
  
  
  
  ############################
  # read data
  ############################
  
  
  dat <- read.table(
    
    text=txt[start:length(txt)],
    
    sep="\t",
    
    header=FALSE,
    
    stringsAsFactors=FALSE,
    
    fill=TRUE
    
  )
  
  
  
  ############################
  # remove empty columns
  ############################
  
  
  dat <- dat[,1:4]
  
  
  
  colnames(dat)<-
    
    c(
      "Genome",
      "Core",
      "Shell",
      "Cloud"
    )
  
  
  
  dat$Core <- as.numeric(dat$Core)
  
  dat$Shell <- as.numeric(dat$Shell)
  
  dat$Cloud <- as.numeric(dat$Cloud)
  
  
  
  dat$Genus <- genus
  
  
  
  return(dat)
  
}



pan <- bind_rows(
  
  read_pan_itol("Genus1"),
  
  read_pan_itol("Genus2"),
  
  read_pan_itol("Genus3")
  
)



write.csv(
  
  pan,
  
  paste0(
    OUTDIR,
    "/Panaroo_core_shell_cloud_all.csv"
  ),
  
  row.names=FALSE
  
)



############################
# 3. Core Shell Cloud Figure
############################


pan_sum <- pan %>%
  
  group_by(Genus)%>%
  
  summarise(
    
    Core=sum(Core),
    
    Shell=sum(Shell),
    
    Cloud=sum(Cloud)
    
  )



pan_long <- pan_sum %>%
  
  pivot_longer(
    
    cols=c(
      Core,
      Shell,
      Cloud
    ),
    
    names_to="Category",
    
    values_to="Number"
    
  )



p1 <- ggplot(
  
  pan_long,
  
  aes(
    
    x=Genus,
    
    y=Number,
    
    fill=Category
    
  )
  
)+
  
  geom_bar(
    
    stat="identity",
    
    position="fill"
    
  )+
  
  
  scale_y_continuous(
    
    labels=scales::percent
    
  )+
  
  
  scale_fill_manual(
    
    values=c(
      
      Core="#2166AC",
      
      Shell="#67A9CF",
      
      Cloud="#B2182B"
      
    )
    
  )+
  
  
  theme_classic(
    
    base_size=14
    
  )+
  
  
  labs(
    
    x=NULL,
    
    y="Genome architecture"
    
  )



ggsave(
  
  paste0(
    OUTDIR,
    "/Core_shell_cloud.pdf"
  ),
  
  p1,
  
  width=5,
  
  height=4
  
)



############################
# 4. Accessory genome
############################


accessory <- pan_sum %>%
  
  mutate(
    
    Accessory=
      
      Shell+Cloud,
    
    Total=
      
      Core+Shell+Cloud,
    
    Accessory_percent=
      
      Accessory/Total*100
    
  )



p2 <- ggplot(
  
  accessory,
  
  aes(
    
    x=Genus,
    
    y=Accessory_percent
    
  )
  
)+
  
  geom_col(
    
    fill="#7B3294"
    
  )+
  
  geom_text(
    
    aes(
      label=paste0(
        round(Accessory_percent,1),
        "%"
      )
    ),
    
    vjust=-0.5
    
  )+
  
  
  theme_classic()+
  
  labs(
    
    x=NULL,
    
    y="Accessory genome (%)"
    
  )



ggsave(
  
  paste0(
    OUTDIR,
    "/Accessory_genome.pdf"
  ),
  
  p2,
  
  width=5,
  
  height=4
  
)



############################
# 5. VF abundance
############################


vf <- read.csv(
  
  VF_FILE,
  
  check.names=FALSE
  
)



rownames(vf)<-
  
  vf[,1]


vf<-vf[,-1]



vf_number <- data.frame(
  
  Genome=colnames(vf),
  
  VF_number=colSums(vf)
  
)



vf_number <- left_join(
  
  vf_number,
  
  genome_group,
  
  by="Genome"
  
)



write.csv(
  
  vf_number,
  
  paste0(
    OUTDIR,
    "/VF_number_tree_classified.csv"
  ),
  
  row.names=FALSE
  
)



p3 <- ggplot(
  
  vf_number,
  
  aes(
    
    x=Genus,
    
    y=VF_number
    
  )
  
)+
  
  
  geom_boxplot()+
  
  geom_jitter(
    
    width=0.15
    
  )+
  
  
  theme_classic(
    
    base_size=14
    
  )+
  
  
  labs(
    
    x=NULL,
    
    y="VF genes"
    
  )



ggsave(
  
  paste0(
    OUTDIR,
    "/VF_repertoire.pdf"
  ),
  
  p3,
  
  width=5,
  
  height=5
  
)



############################
# 6. BGC abundance
############################


bgc <- read.csv(
  
  BGC_FILE,
  
  check.names=FALSE
  
)



bgc_long <- bgc %>%
  
  pivot_longer(
    
    cols=-Genome,
    
    names_to="BGC_type",
    
    values_to="Number"
    
  )



bgc_long <- left_join(
  
  bgc_long,
  
  genome_group,
  
  by="Genome"
  
)



bgc_sum <- bgc_long %>%
  
  group_by(
    
    Genus,
    
    BGC_type
    
  )%>%
  
  summarise(
    
    Number=sum(Number),
    
    .groups="drop"
    
  )



p4 <- ggplot(
  
  bgc_sum,
  
  aes(
    
    x=Genus,
    
    y=Number,
    
    fill=BGC_type
    
  )
  
)+
  
  geom_col()+
  
  
  theme_classic(
    
    base_size=14
    
  )+
  
  
  labs(
    
    x=NULL,
    
    y="BGC clusters"
    
  )



ggsave(
  
  paste0(
    OUTDIR,
    "/BGC_landscape.pdf"
  ),
  
  p4,
  
  width=6,
  
  height=5
  
)



############################
# 7. Combine main figure
############################


main <-
  
  p1 /
  
  p2 |
  
  p3 /
  
  p4



ggsave(
  
  paste0(
    OUTDIR,
    "/NatureCom_MainFigure_V4.pdf"
  ),
  
  main,
  
  width=12,
  
  height=10
  
)



cat(
  "
================================
NatureCom Main Figure V4 finished
================================
"
)
