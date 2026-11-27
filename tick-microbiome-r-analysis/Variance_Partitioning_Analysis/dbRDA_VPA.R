library(phyloseq)
library(vegan)
library(readxl)
library(ggplot2)
library(ggVennDiagram)
library(ggrepel)
library(ggforce)
library(dplyr)

##Bray-Curtis matrix
family <- read_excel("Family.xlsx")
meta <- read_excel("metadata.xlsx")

sample_id <- family$SampleID
rownames(family) <- family$SampleID

head(rownames(family))


# Hellinger transformation
family_hel <- decostand(
  family,
  method="hellinger"
)


# Bray-Curtis
bc.dist <- vegdist(
  family_hel,
  method="bray"
)

#PERMANOVA marginal test
adonis_margin <- adonis2(
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

adonis_margin
summary(adonis_margin)
res_df <- as.data.frame(adonis_margin)
write.csv(res_df, file = "adonis_margin.csv", row.names = TRUE)

#dbRDA
dbRDA <- capscale(
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

anova(dbRDA,
      permutations=999)


db <- anova(
  dbRDA,
  by="terms",
  permutations=999
)

db

summary(db)
db <- as.data.frame(db)
write.csv(db, file = "dbRDA_Anova.csv", row.names = TRUE)

#Factor session
clade_colors <- c(
  "A"="#0072B2",
  "B"="#E69F00",
  "C"="#009E73"
)

site_scores <- scores(
  dbRDA,
  display="sites"
)

site_scores <- as.data.frame(site_scores)

rownames(site_scores)<-metadata$SampleID
colnames(site_scores)


plot_data <- cbind(
  site_scores,
  meta[match(
    rownames(site_scores),
    meta$SampleID
  ),]
)

colnames(plot_data)
table(plot_data$Clade)
#提取环境箭头
env_scores <- scores(
  dbRDA,
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
eig <- summary(dbRDA)$cont$importance

xlab <- paste0(
  "dbRDA1 (",
  round(eig[2,1]*100,1),
  "%)"
)


ylab <- paste0(
  "dbRDA2 (",
  round(eig[2,2]*100,1),
  "%)"
)


#plot
pdf("p_dbRDA.pdf", width = 10, height = 8)
p_dbRDA <- ggplot(
  plot_data,
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
p_dbRDA

dev.off()

##Variation Partitioning Analysis
#three matrixes
host <- meta[,c(
  "Clade",
  "Stage"
)]

env <- meta[,c(
  "Temperature",
  "Humidity",
  "Altitude"
)]

space <- meta[,c(
  "Longitude",
  "Latitude"
)]

vpa <- varpart(
  family_hel,
  host,
  env,
  space
)


vpa

#VPA data management
vpa_df <- data.frame(
  
  Component=c(
    "Clade (pure)",
    "Environment (pure)",
    "Spatial (pure)",
    "Shared",
    "Residual"
  ),
  
  
  Variation=c(
    25.6,
    6.8,
    4.5,
    12.3,
    50.8
  )
  
)


vpa_df$Component <- factor(
  vpa_df$Component,
  levels=
    rev(vpa_df$Component)
)


#plot
vpa_colors <- c(
  "Clade (pure)"="#0072B2",
  "Environment (pure)"="#E69F00",
  "Spatial (pure)"="#009E73",
  "Shared"="#999999",
  "Residual"="#D9D9D9"
)

pdf("VPA_plot.pdf", width = 10, height = 8)
p_vpa <- ggplot(
  vpa_df,
  aes(
    x=Variation,
    y=Component,
    fill=Component
  )
)+
  
  
  geom_col(
    width=0.65
  )+
  
  
  geom_text(
    aes(
      label=paste0(
        Variation,
        "%"
      )
    ),
    hjust=-0.15,
    size=5
  )+
  
  
  scale_fill_manual(
    values=vpa_colors
  )+
  
  
  scale_x_continuous(
    expand=
      expansion(
        mult=c(0,0.1)
      )
  )+
  
  
  labs(
    x="Explained variation (%)",
    y=NULL
  )+
  
  
  theme_classic(
    base_size=16
  )+
  
  
  theme(
    
    legend.position="none",
    
    axis.text=
      element_text(
        color="black"
      ),
    
    axis.title.x=
      element_text(
        face="bold"
      )
    
  )

p_vpa

dev.off()
#plot combination
library(patchwork)
final_fig <- 
  p_dbRDA /
  p_vpa +
  
  plot_layout(
    heights=c(2,1)
  )

final_fig

ggsave(
  "dbRDA_VPA.pdf",
  final_fig,
  width=9,
  height=20,
  units="in"
)

#Distance-decay analysis
library(geosphere)

geo.dist <- distm(
  meta[,c(
    "Longitude",
    "Latitude"
  )],
  fun=distHaversine
)

geo.dist <- as.dist(
  geo.dist
)
summary(geo.dist)
geo.dist <- geo.dist/1000

#check datasets
labels(bc.dist)
rownames(meta)


#Mantel distance-decay test
mantel_decay <- mantel(
  bc.dist,
  geo.dist,
  method="spearman",
  permutations=9999
)

mantel_decay
mantel_all_results <- data.frame(
  
  Test="Overall",
  
  Mantel_r=
    as.numeric(mantel_decay$statistic),
  
  P_value=
    mantel_decay$signif
  
)

write.csv(mantel_all_results, "Mantel_Distance_Decay_Test.csv", row.names=FALSE)

#plot distance-decay curve
decay_df <- data.frame(
  
  Microbial_distance =
    as.vector(bc.dist),
  
  Geographic_distance =
    as.vector(geo.dist)
  
)


p_decay <- ggplot(
  decay_df,
  aes(
    x=Geographic_distance,
    y=Microbial_distance
  )
)+
  
  
  geom_point(
    size=3,
    alpha=0.5,
    color="#0072B2"
  )+
  
  
  geom_smooth(
    method="lm",
    formula=y~x,
    color="#D55E00",
    linewidth=1.2,
    fill="#D55E00",
    alpha=0.2
  )+
  
  
  theme_classic(
    base_size=16
  )+
  
  
  labs(
    
    x="Geographic distance (km)",
    
    y="Bray–Curtis dissimilarity"
    
  )

p_decay

ggsave(
  "Mantel_Distance_Decay_Curve.pdf",
  p_decay,
  width=9,
  height=10,
  units="in"
)

#adding clades information
mantel_clade_results <- data.frame(
  Clade = character(),
  Mantel_r = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for(i in unique(meta$Clade)){
  idx <- which(meta$Clade == i)
      # microbiome distance
  bc_sub <- as.dist(
    as.matrix(bc.dist)[idx,idx]
  )
    # geographic distance
  geo_sub <- as.dist(
    as.matrix(geo.dist)[idx,idx]
  )
    # Mantel test
  mantel_result <- mantel(
    bc_sub,
    geo_sub,
    method="spearman",
    permutations=9999
  )
    # 保存结果
  mantel_clade_results <- rbind(
    mantel_clade_results,
    data.frame(
      Clade=i,
      Mantel_r=as.numeric(mantel_result$statistic),
      P_value=mantel_result$signif
    )
  )
  
}


mantel_clade_results
write.csv(mantel_clade_results, "Clade_specific_distance_decay_Mantel.csv",row.names=FALSE)

#Partial Mantel Analysis
clade.dist <- dist(
  as.numeric(
    factor(meta$Clade)
  )
)

partial_result <- mantel.partial(
  bc.dist,
  geo.dist,
  clade.dist,
  method="spearman",
  permutations=9999
)

partial_mantel_results <- data.frame(
  Test="Geographic distance controlled by Clade",
  Mantel_r=
    as.numeric(partial_result$statistic),
  P_value=
    partial_result$signif
)

partial_mantel_results
write.csv(partial_mantel_results, "Partial_Mantel_control_Clade.csv", row.names=FALSE)

#plot overall and partial mantel distance curve
#data combination 
overall_df <- data.frame(
  Geographic_distance = as.vector(geo.dist),
  Microbial_distance = as.vector(bc.dist),
  Clade = "Overall"
)

clade_df <- data.frame()

for(i in unique(meta$Clade)){
  idx <- which(meta$Clade == i)

  # microbiome distance
  bc_sub <- as.matrix(bc.dist)[idx,idx]
  
  # geographic distance
  geo_sub <- as.matrix(geo.dist)[idx,idx]
  
  # 去除对角线0
  keep <- lower.tri(
    bc_sub
  )
  
  temp <- data.frame(
    Geographic_distance =
      geo_sub[keep],
    Microbial_distance =
      bc_sub[keep],
    Clade=i
  )
  
  clade_df <- rbind(
    clade_df,
    temp
  )
  
}

decay_df <- rbind(
  overall_df,
  clade_df
)

table(decay_df$Clade)

#Calculating linear regression parameters
library(broom)
lm_results <- decay_df %>%
  group_by(Clade) %>%
  do(
    tidy(
      lm(
        Microbial_distance~
          Geographic_distance,
        data=.
      )
    )
  )

lm_results


all_mantel_results <- data.frame(
  Clade=character(),
  Mantel_r=numeric(),
  P_value=numeric()
)


# Overall

m_all <- mantel(
  bc.dist,
  geo.dist,
  method="spearman",
  permutations=9999
)


all_mantel_results <- rbind(
  all_mantel_results,
  data.frame(
    Clade="Overall",
    Mantel_r=m_all$statistic,
    P_value=m_all$signif
  )
)



# Clade-specific

for(i in unique(meta$Clade)){
  idx <- which(meta$Clade==i)
  bc_sub <- as.dist(
    as.matrix(bc.dist)[idx,idx]
  )
  
  geo_sub <- as.dist(
    as.matrix(geo.dist)[idx,idx]
  )
  
  m <- mantel(
    bc_sub,
    geo_sub,
    method="spearman",
    permutations=9999
  )
  
  all_mantel_results <- rbind(
    all_mantel_results,
    data.frame(
      Clade=i,
      Mantel_r=m$statistic,
      P_value=m$signif
    )
    
  )
  
}


all_mantel_results

#Plot curves
colors <- c(
  "Overall"="black",
  "A"="#0072B2",
  "B"="#D55E00",
  "C"="#009E73"
)

p_decay_all <- ggplot(
  decay_df,
  aes(
    x=Geographic_distance,
    y=Microbial_distance,
    color=Clade
  )
)+
  
  # 点
  geom_point(
    alpha=0.25,
    size=1.8
  )+
  
  
  # regression
  geom_smooth(
    method="lm",
    se=TRUE,
    linewidth=1.2
  )+
  
  scale_color_manual(
    values=colors
  )+
  
  theme_classic(
    base_size=16
  )+
  
  
  labs(
    x="Geographic distance (km)",
    y="Bray–Curtis dissimilarity",
    color=NULL
    
  )+
  
  
  theme(
    legend.position="right",
    axis.title=
      element_text(
        face="bold"
      ),
    panel.border=
      element_rect(
        fill=NA,
        color="black",
        linewidth=0.8
      )
    
  )


p_decay_all

label_df <- data.frame(
  
  Clade=c(
    "Overall",
    "A",
    "B",
    "C"
  ),
  
  label=c(
    "Overall: r=0.25, P=0.001",
    "A: r=0.10, P=0.12",
    "B: r=0.32, P=0.005",
    "C: r=0.08, P=0.45"
  )
  
)

p_decay_all <- p_decay_all +
  geom_text(
    data=label_df,
    aes(
      x=600,
      y=c(
        0.66,
        0.64,
        0.62,
        0.6
      ),
      label=label
    ),
    inherit.aes=FALSE,
    hjust=0,
    size=4
  )

ggsave(filename = "Combination_Distance_Decay_Curves.pdf", plot = p_decay_all, width=10,height=8, units = "in")
