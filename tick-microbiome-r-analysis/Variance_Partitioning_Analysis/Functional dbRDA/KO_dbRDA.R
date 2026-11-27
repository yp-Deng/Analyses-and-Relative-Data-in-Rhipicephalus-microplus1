library(phyloseq)
library(vegan)
library(readxl)
library(ggplot2)
library(ggVennDiagram)
library(ggrepel)
library(ggforce)
library(dplyr)
library(tibble)

##Bray-Curtis matrix
meta <- read_excel("metadata.xlsx")
ko_raw <- read_excel("KO.xlsx")

ko_mat <- ko_raw %>%
  column_to_rownames(var = "SampleID") 

head(rownames(KO))

#Filter KO lower than 1e-6
ko_filt <- ko_mat[, colSums(ko_mat) > 1e-4]


# Hellinger transformation
ko_hel <- decostand(
  ko_filt,
  method="hellinger"
)


# Bray-Curtis
bc.dist <- vegdist(
  ko_hel,
  method="bray"
)

#PERMANOVA marginal test
ko_adonis_margin <- adonis2(
  bc.dist ~ 
    Clade+
    Stage+
    Temperature+
    Humidity+
    Altitude+
    Latitude+
    Longitude,
  data=meta,
  permutations=999,
  by="margin"
)

ko_adonis_margin
summary(ko_adonis_margin)
ko_adonis_margin <- as.data.frame(ko_adonis_margin)
write.csv(ko_adonis_margin, file = "ko_adonis_margin.csv", row.names = TRUE)

#dbRDA
ko_dbRDA <- capscale(
  bc.dist~
    Clade+
    Stage+
    Temperature+
    Humidity+
    Altitude+
    Latitude+
    Longitude,
  data=meta,
  add=TRUE
)

anova(ko_dbRDA,
      permutations=999)


ko_db <- anova(
  ko_dbRDA,
  by="terms",
  permutations=999
)

ko_db

summary(ko_db)
ko_db <- as.data.frame(ko_db)
write.csv(ko_db, file = "ko_dbRDA_Anova.csv", row.names = TRUE)


clade_colors <- c(
  "A"="#0072B2",
  "B"="#E69F00",
  "C"="#009E73"
)

ko_site_scores <- scores(
  ko_dbRDA,
  display="sites"
)

ko_site_scores <- as.data.frame(ko_site_scores)

head(rownames(ko_site_scores))
# rownames(ko_site_scores)<-metadata$SampleID
colnames(ko_site_scores)


ko_plot_data <- cbind(
  ko_site_scores,
  meta[match(
    rownames(ko_site_scores),
    meta$SampleID
  ),]
)

colnames(ko_plot_data)
table(ko_plot_data$Clade)

#提取环境箭头
env_scores <- scores(
  ko_dbRDA,
  display="bp"
)


env_scores <- as.data.frame(env_scores)

env_scores$Variable <- rownames(env_scores)

env_scores <- env_scores[
  !grepl(
    "Clade|Stage",
    env_scores$Variable
  ),
]

#Computing interprete rate
ko_eig <- summary(ko_dbRDA)$cont$importance

xlab <- paste0(
  "dbRDA1 (",
  round(ko_eig[2,1]*100,1),
  "%)"
)


ylab <- paste0(
  "dbRDA2 (",
  round(ko_eig[2,2]*100,1),
  "%)"
)


#plot
pdf("p_dbRDA.pdf", width = 10, height = 8)
p_dbRDA_ko <- ggplot(
  ko_plot_data,
  aes(
    x=CAP1,
    y=CAP2,
    color=Clade
  )
)+
  
  
  # 0轴
  
  geom_hline(
    yintercept=0,
    linetype="dashed",
    color="grey75"
  )+
  
  geom_vline(
    xintercept=0,
    linetype="dashed",
    color="grey75"
  )+
  
  
  # 95% ellipse
  
  stat_ellipse(
    aes(fill=Clade),
    type="norm",
    geom="polygon",
    alpha=0.15,
    color=NA
  )+
  
  
  # points
  
  geom_point(
    size=3.5,
    alpha=0.9
  )+
  
  
  # 环境箭头
  
  geom_segment(
    data=env_scores,
    aes(
      x=0,
      y=0,
      xend=CAP1*1.2,
      yend=CAP2*1.2
    ),
    inherit.aes=FALSE,
    arrow=
      arrow(
        length=unit(0.25,"cm")
      ),
    linewidth=0.8,
    color="black"
  )+
  
  
  # 环境变量名称
  
  geom_text_repel(
    data=env_scores,
    aes(
      x=CAP1*1.35,
      y=CAP2*1.35,
      label=Variable
    ),
    inherit.aes=FALSE,
    size=5,
    fontface="italic",
    color="black"
  )+
  
  
  scale_color_manual(
    values=clade_colors
  )+
  
  
  scale_fill_manual(
    values=clade_colors
  )+
  
  
  labs(
    x=xlab,
    y=ylab,
    color="Clade"
  )+
  
  
  theme_classic(
    base_size=16
  )+
  
  
  theme(
    
    legend.position="right",
    
    legend.title=
      element_text(
        face="bold"
      ),
    
    legend.text=
      element_text(
        size=14
      ),
    
    axis.title=
      element_text(
        face="bold"
      ),
    
    panel.border=
      element_rect(
        color="black",
        fill=NA,
        linewidth=0.8
      )
    
  )
p_dbRDA_ko

dev.off()

