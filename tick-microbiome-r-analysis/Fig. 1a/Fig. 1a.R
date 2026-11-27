#Fig. 1a
# 加载必要的包
library(tidyverse)
library(gghalves)
library(circlize)

# 读取数据
subject_info  = readxl::read_xlsx("data.xlsx")

# 设置颜色
age_color <- '#631879FF'
bmi_color <- '#FF410DFF'
sex_color <- c("Female" = '#F2300F', "Male" = '#E1BD6D')
ethnicity_color <- c("Caucasian" = '#ECCBAE', "Asian" = '#0B775E', #046C9A 文章配色
                     "Hispanics" = '#D69C4E', "Black" = '#81A88D',
                     "Unknown" = 'grey') #ABDDDE 文章配色
iris_color <- c("IS" = '#F21A00', "IR" = '#74A089', "Unknown" = 'grey') #3B9AB2 文章配色


# 数据处理
df <-
  data.frame(
    factors = subject_info$subject_id_random,
    x = 1,
    y = 1,
    subject_info,
    stringsAsFactors = TRUE
  ) %>%
  arrange(adjusted_age) %>%
  mutate(factors = factor(factors, levels = factors))


# 定义 circos 画布
pdf("plot_circos.pdf", width = 10, height = 10)
circos.clear()
circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  circle.margin = c(0.1, 0.1, 0.1, 0.1),
  cell.padding = c(0, 0, 0, 0)
)
circos.initialize(factors = df$factors, x = df$x, xlim = c(0.5, 1.5))



# 定义绘制函数
plot_circos_track <- function(temp_value, color, ylab, ylim_mult = c(0.8, 1.1), track_height = 0.2) {
  circos.track(
    factors = df$factors,
    y = temp_value,
    ylim = range(temp_value, na.rm = TRUE) * ylim_mult,
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      
      circos.yaxis(
        side = "left",
        at = c(ceiling(0.8*min(temp_value, na.rm = TRUE)),
               round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)) / 2, 0),
               round(max(temp_value, na.rm = TRUE), 0)),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.8,
        labels.niceFacing = FALSE
      )
      
      circos.lines(
        x = mean(xlim, na.rm = TRUE),
        y =  temp_value[i],
        type = "h",
        col = color,
        lwd = 2
      )
      
      if (ylab == "Age") {
        circos.text(
          x = 1,
          y = 110,
          labels = name,
          facing = "clockwise",
          niceFacing = TRUE,
          cex = 0.75
        )
      }
      
      circos.points(
        x = mean(xlim),
        y =  temp_value[i],
        pch = 16,
        cex = 0.8,
        col = color
      )
    }
  )
}

# 绘制 tracks
plot_circos_track(df$adjusted_age, age_color, "Age")
plot_circos_track(df$BMI, bmi_color, "BMI")

# 定义绘制分类变量的函数
plot_categorical_track <- function(temp_var, colors, track_height = 0.1) {
  temp_var[is.na(temp_var)] <- "grey"
  mapped_colors <- colors[temp_var]
  circos.track(
    factors = df$factors,
    y = df$y,
    ylim = c(0, 1),
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(x, y) {
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      circos.rect(
        xleft = xlim[1],
        ybottom = ylim[1],
        xright = xlim[2],
        ytop = ylim[2],
        col = mapped_colors[i],
        bg.border = "black"
      )
    }
  )
}

# 绘制分类 tracks
plot_categorical_track(df$Gender, sex_color)
plot_categorical_track(df$Ethnicity, ethnicity_color)
plot_categorical_track(df$IRIS, iris_color)

dev.off()


# 定义绘制箱线图和点图的函数
plot_box_dot <- function(data, color, binwidth) {
  ggplot(data.frame(class = "class", value = data), aes(x = class, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(binaxis = "y", color = color, fill = color, 
                 shape = 16, binwidth = binwidth, stackdir = "center") +
    theme_bw() +
    labs(x = "", y = "") +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
}

# 绘制年龄和BMI的箱线图和点图
plot_age <- plot_box_dot(df$adjusted_age, age_color, 1)
plot_bmi <- plot_box_dot(df$BMI, bmi_color, 0.6)

# 打印绘图
print(plot_age)
print(plot_bmi)

ggsave(plot_age, filename = "plot_age.pdf", width = 3, height = 10)
ggsave(plot_bmi, filename = "plot_bmi.pdf", width = 3, height = 10)


# 定义绘制条形图的函数
plot_bar <- function(data, colors, levels_order) {
  ggplot(data.frame(class = "class", value = data), aes(x = class)) +
    geom_bar(
      aes(fill = factor(value, levels = levels_order)),
      color = "black", 
      position = "stack", 
      show.legend = FALSE, 
      width = 2
    ) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    labs(x = "", y = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()
    )
}

# 绘制性别、种族和IRIS的条形图
plot_sex <- plot_bar(df$Gender, sex_color, c("Female", "Male"))
plot_ethnicity <- plot_bar(df$Ethnicity, ethnicity_color, c("Asian", "Black", "Caucasian", "Hispanics"))
plot_iris <- plot_bar(df$IRIS, iris_color, c("IR", "IS", "Unknown"))

# 打印绘图
print(plot_sex)
print(plot_ethnicity)
print(plot_iris)

ggsave(plot_sex, filename = "plot_sex.pdf", width = 1.5, height = 10)
ggsave(plot_ethnicity, filename = "plot_ethnicity.pdf", width = 1.5, height = 10)
ggsave(plot_iris, filename = "plot_iris.pdf", width = 1.5, height = 10)



物种稀释曲线
#设置目录
setwd("C:/Users/ghost/Desktop/2023-待办事情/Data_Collection_Figs/RareCurve")

#载入R包
library(devtools)
library(amplicon)
library(vegan)
library("BiocManager")
library(dplyr)

#读入分析数据和分组信息等
data <- read.table("feature_rare.txt", sep = "\t", header =  T, stringsAsFactors = FALSE, row.names = 1)
design <- read.table("metadata.txt",  sep = "\t", header =  T, stringsAsFactors = FALSE, row.names = 1)

#通过amplicon中的alpha_rare_curve函数绘制稀释曲线
p = alpha_rare_curve(data, design, groupID = "group")
p

#将输出文件保存
ggsave(paste0("p1.rare_curve.pdf"), p, width=89*1.5, height=80, units="mm")
ggsave(paste0("p1.rare_curve.png"), p, width=89*1.5, height=80, units="mm")

#其他可以参照https://shengxin.ren/article/255
#分析可以参照

#species accumulation curve
A <- data.frame(t(read.table("A.txt",sep = "\t", header =  T, stringsAsFactors = FALSE, row.names = 1)))
B <- data.frame(t(read.table("B.txt",sep = "\t", header =  T, stringsAsFactors = FALSE, row.names = 1)))
CC <- data.frame(t(read.table("C.txt",sep = "\t", header =  T, stringsAsFactors = FALSE, row.names = 1)))

sp1 <- specaccum(A)
sp11 <- specaccum(A,"random") #计算标准差
summary(sp11)
spB <- specaccum(B)
spB1 <- specaccum(B,"random") #计算标准差
summary(spB1)
spCC <- specaccum(CC)
spCC1 <- specaccum(CC,"random") #计算标准差
summary(spCC1)

p1 <- plot(spB, col = "grey", lwd = 2, lty = 2, xlab = "Number of Samples", ylab = "Number of Species", ci = 0, ylim = c(1800,6000))+
  #boxplot(spC1, add = TRUE, col = "grey")+
  lines(sp1, col = "steelblue", lwd = 2, lty = 2, ci = 0)+
  #boxplot(spY1, add = TRUE, col = "steelblue")+
  lines(spC, col = "coral2", lwd = 2, lty = 2, ci = 0)+
  #boxplot(spR1, add = TRUE, col = "coral2")+
  legend(21,6000, legend = c("clade A","clade B","clade C"), col = c("steelblue","grey","coral2"), lwd = 2:2)

p1 <- plot(spB, col = "#d53f50", lwd = 3, lty = 2, xlab = "Number of Samples", ylab = "Number of Species", ci = 0, ylim = c(1800,6000))+
  #boxplot(spB1, col = "#d53f50")+
  lines(sp1, col = "#f46e44", lwd = 3, lty = 2, ci = 0)+
  #boxplot(sp11, col = "#f46e44")+
  lines(spCC, col = "#0f6bac", lwd = 3, lty = 2, ci = 0)+
  #boxplot(spCC1, col = "#0f6bac")+
  lines(spC, col = "#86b6fe", lwd = 3, lty = 2, ci = 0)+
  #boxplot(spC1, col = "#86b6fe")+
  lines(spR, col = "#035f74", lwd = 3, lty = 2, ci = 0)+
  #boxplot(spR1, col = "#035f74")+
  lines(spY, col = "#fee08b", lwd = 3, lty = 2, ci = 0)+
  #boxplot(spY1, col = "#fee08b")+
  legend(0,6000, legend = c("clade B","clade A","clade C","adult","nymph","larva"), col = c("#d53f50","#f46e44","#0f6bac","#86b6fe","#035f74","#fee08b"), lwd = 2:2)


ggsave(paste0("p1.accu_curve.pdf"), p, width=89*1.5, height=80, units="mm")
ggsave(paste0("p1.accu_curve.png"), p, width=89*1.5, height=80, units="mm")