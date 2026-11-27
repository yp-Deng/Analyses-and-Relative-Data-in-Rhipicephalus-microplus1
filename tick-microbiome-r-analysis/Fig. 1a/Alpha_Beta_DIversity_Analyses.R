##Alpha多样性
library(amplicon)
library(agricolae)
library(dplyr)
library(maptools)
library(ggplot2)
library(patchwork) #组图

setwd("D:/BaiduNetdiskWorkspace/2023-待办事情/Data_Collection_Figs/alpha_beta_div/alpha")

gene_data <- read.table("gene_alpha.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)
COG_data <- read.table("COG_alpha.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)
family_data <- read.table("family_alpha.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)
KO_data <- read.table("KO_alpha.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)
species_data <- read.table("species_alpha.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)
design <- read.table("meta_clade.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)



if (TRUE){
  sub_design = subset(design, group %in% c("A","B","C"))
  sub_design$group  = factor(sub_design$group, levels=c("A","B","C"))
}else{
  sub_design = design
}

idx = rownames(sub_design) %in% rownames(gene_data)
sub_design=sub_design[idx,]
sub_gene=gene_data[rownames(sub_design),]

index = cbind(sub_gene, sub_design) 

m = "shannon_e"

model = aov(index[[m]] ~ group, data=index)
Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_table = as.data.frame(Tukey_HSD$group) 

write.table(paste(m, "\n\t", sep=""), file=paste(1, "alpha_",m,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(Tukey_HSD_table, file=paste("alpha_",m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

out = LSD.test(model,"group", p.adj="none") # alternative fdr
stat = out$groups
index$stat=stat[as.character(index$group),]$groups
max=max(index[,c(m)])
min=min(index[,c(m)])
x = index[,c("group",m)]
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
y=as.data.frame(y)
rownames(y)=y$group
index$y=y[as.character(index$group),]$Max + (max-min)*0.05

p = ggplot(index, aes(x=group, y=index[[m]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(m, "index")) + theme_classic() + main_theme +
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)

if (length(unique(sub_design$group))>3){
  p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
}

color = c("#00BFC4", "#F9766E", "#9E7C61", "#A9B7C6", "#7AA1D2", "#718096","#6D98BA")
p = p + scale_color_manual(values = color)

ggsave(paste(1,"alpha_", m, ".pdf", sep=""), p, width = 5, height = 3)





##Beta多样性
library(vegan)
library(ggplot2)
library(patchwork)


## Plotting settings ------------------------

groups <- read.table('group_info.txt', header = FALSE, colClasses=c("character","character"))
colnames(groups) <- c("sample", "group")

length  <- length(unique(as.character(groups$sample)))
length1 <- length(unique(as.character(groups$group)))
times1  <- length%/%8
res1    <- length%%8
times2  <- length%/%5
res2    <- length%%5

# dot color and number
mycol = c("#B2182B","#E69F00","#56B4E9","#009E73",
          "#F0E442","#0072B2","#D55E00","#CC79A7",
          "#CC6666","#9999CC","#66CC99","#99999","#ADD1E5")
col1=rep(mycol,times1)
col=c(col1,mycol[1:res1])

# dot shape and number
pich1 <- rep(c(15:18,20,7:14,0:6),times2)
pich  <- c(pich1,15:(15+res2))

# legend column number
if(length1>30){
  n = 2
}else{
  n = 1
}

## PCoA analysis with vegan ------------------------
data <- read.table('input_data.xls', header = TRUE, 
                   row.names = 1, sep = "\t", comment.char = "", check.names = FALSE)
data <- t(data)

distance <- vegdist(data, method = 'bray') #calculate distance
write.csv(as.matrix(distance), 'results\\samples_distance.csv', quote = F)

pcoa <- cmdscale(distance, k = 5, eig = TRUE)
points <- data.frame(scores(pcoa))
point  <- data.frame(sample = row.names(points), points) 
write.csv(point, file = "results\\PCoA_coordinates.csv", row.names = FALSE)

pc1 <- round((pcoa$eig/sum(pcoa$eig))*100,2)[1]
pc2 <- round((pcoa$eig/sum(pcoa$eig))*100,2)[2]


## PCoA plot using ggplot2 ------------------------
colnames(points)[1:2] <- c('dim1', 'dim2')
plotdata = data.frame(rownames(points), points$dim1, points$dim2, groups$group)
colnames(plotdata)=c("sample","dim1","dim2","group")


## p value ------------------------
adonis_result_dis <- adonis2(distance ~ group, data = groups, permutations = 999)

R2 = adonis_result_dis$R2[1]
pvalue = adonis_result_dis$`Pr(>F)`[1]

adonis <- paste("PERMANOVA:\nR2 = ", round(R2,4), "\nP-value = ", pvalue)


## PCoA plot + boxplot using ggplot2 ------------------------
plot <- ggplot(plotdata, aes(dim1, dim2)) + 
  geom_point(aes(colour = group, shape = group), size = 6) + 
  #geom_text(aes(label = sample), size = 4, family = "Arial", hjust = 0.5, vjust = -1)+ # display sample names
  scale_shape_manual(values = pich) + 
  scale_colour_manual(values = col) +
  #labs(title = "PCoA Plot") + 
  xlab(paste("PCoA axis1 ( ",pc1,"%"," )", sep = "")) + 
  ylab(paste("PCoA axis2 ( ",pc2,"%"," )", sep = "")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, colour = "black", face = "bold")) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        panel.grid = element_blank(), 
        axis.title = element_text(color = 'black', family = "Arial", size = 20),
        axis.text = element_text(colour = 'black', size = 16, margin = unit(0.6, "lines")),
        axis.ticks = element_line(color = 'black'), 
        axis.ticks.length = unit(0.4, "lines"), 
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 14),
        legend.key = element_blank(),
        legend.position = c(0.9, 0.1),
        legend.background = element_rect(colour = "black")) +
  guides(col = guide_legend(ncol = n), shape = guide_legend(ncol = n)) + 
  stat_ellipse(aes(x = dim1,y = dim2, color = group), data = plotdata) # group cluster

box1 <- ggplot(plotdata, 
               aes(x = groups$group, y = dim1, fill = groups$group)) +
  geom_boxplot(show.legend = FALSE) + 
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
  geom_jitter(show.legend = FALSE) + 
  scale_fill_manual(values = col) +
  coord_flip() + 
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color='black'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 16))

box2 <- ggplot(plotdata,
               aes(x = groups$group, y = dim2, fill = groups$group)) + 
  geom_boxplot(show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
  geom_jitter(show.legend = FALSE) + 
  scale_fill_manual(values = col) +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color='black'),
        axis.text.x = element_text(colour ='black', size = 16, angle = 45,
                                   vjust = 1, hjust = 1), 
        axis.text.y=element_blank())

box3 <- ggplot(plotdata, aes(dim1, dim2)) +
  geom_text(aes(x = -0.5,y = 0.6,
                label = adonis), 
            size = 4) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p <- box1 + box3 + plot + box2 +
  plot_layout(heights = c(1,4), widths = c(4,1), ncol = 2, nrow = 2)

cairo_pdf(filename = "results\\PCoA_box_dot_combine.pdf", height = 8, width = 8)
p
dev.off()

png(filename = "results\\PCoA_box_dot_combine.png",res = 600, height = 5400, width = 5400, type = "cairo")
p
dev.off()