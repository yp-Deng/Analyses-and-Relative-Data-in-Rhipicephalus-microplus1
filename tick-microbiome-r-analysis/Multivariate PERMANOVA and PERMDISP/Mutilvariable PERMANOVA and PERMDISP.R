#加载R包
library(vegan)
library(readxl)
library(dplyr)

#Changing "Clade" to "Stage" when analyzing "stage" 

metadata <- read_excel(
  "metadata.xlsx"
)


otu <- read_excel(
  "OTU.xlsx"
)


# 第一列改名

colnames(otu)[1]="SampleID"



# 匹配SampleID

otu <- otu[
  match(
    metadata$SampleID,
    otu$SampleID
  ),
]


# OTU丰度矩阵

otu_matrix <- otu %>%
  select(-SampleID)



# 全部转换numeric

otu_matrix <- as.data.frame(
  lapply(
    otu_matrix,
    as.numeric
  )
)



# NA处理

otu_matrix[
  is.na(otu_matrix)
] <- 0



# matrix

otu_matrix <- as.matrix(
  otu_matrix
)



# Bray-Curtis

bray_dist <- vegdist(
  otu_matrix,
  method="bray"
)


#PERMANOVA analysis
permanova_result <- adonis2(
  
  bray_dist ~
    
    Clade +
    
    `Humidity(%)` +
    
    `Temperature(℃)` +
    
    `Altitude (m)` +
    
    `Longtitude (E)` +
    
    `Latitude (N)`+
    
    #`Stage` ,
  
  data = metadata2,
  
  permutations = 999,
  
  method = "bray"
  
)


print(permanova_result)

#Marginal effects of environmental variables
permanova_marginal <- adonis2(
  
  bray_dist ~
    
    Clade +
    
    `Humidity(%)` +
    
    `Temperature(℃)` +
    
    `Altitude (m)` +
    
    `Longtitude (E)` +
    
    `Latitude (N)`+
    
    #`Stage` ,
  
  data=metadata2,
  
  permutations=999,
  
  by="margin"
  
)


print(permanova_marginal)

#PERMDISP test
dispersion <- betadisper(
  
  bray_dist,
  
  metadata2$Clade
  
)

permdisp_result <- permutest(
  
  dispersion,
  
  permutations = 999
  
)


print(permdisp_result)

#Plot
boxplot(
  
  dispersion,
  
  xlab="Clade",
  
  ylab="Distance to centroid",
  
  main="PERMDISP among Clades"
  
)

nmds <- metaMDS(
  
  otu_matrix,
  
  distance="bray",
  
  k=2,
  
  trymax=100
  
)


plot(
  nmds,
  type="n"
)


points(
  
  nmds,
  
  display="sites",
  
  col=as.factor(metadata2$Clade),
  
  pch=19
  
)


legend(
  
  "topright",
  
  legend=levels(
    as.factor(metadata2$Clade)
  ),
  
  col=1:3,
  
  pch=19
  
)

#Partial Mantel Test
library(ecodist)
library(geosphere)

bray_dist

bray_dist <- vegdist(
  otu_matrix,
  method="bray"
)

#using longtitude and latitude as variables
env_data <- metadata2 %>%
  select(
    `Humidity(%)`,
    `Temperature(℃)`,
    `Altitude (m)`
  )

#标准化单位
env_scaled <- scale(
  env_data
)
#计算Euclidean距离
env_dist <- dist(
  env_scaled,
  method="euclidean"
)
#转换为矩阵
env_dist <- as.dist(
  env_dist
)

#建立坐标
coordinates <- metadata2 %>%
  select(
    `Longtitude (E)`,
    `Latitude (N)`
  )

# POSIXct转换回Excel数值

coordinates <- coordinates %>%
  mutate(
    
    longitude = as.numeric(
      as.Date(`Longtitude (E)`) -
        as.Date("1899-12-30")
    ),
    
    
    latitude = as.numeric(
      as.Date(`Latitude (N)`) -
        as.Date("1899-12-30")
    )
    
  ) %>%
  
  select(
    longitude,
    latitude
  )

library(geosphere)


#计算地理距离
geo_dist <- distm(
  
  coordinates,
  
  fun = distHaversine
  
)

# 转换为km

geo_dist <- geo_dist / 1000


geo_dist <- as.dist(
  geo_dist
)

#检查矩阵
dim(as.matrix(bray_dist))

dim(as.matrix(env_dist))

dim(as.matrix(geo_dist))

summary(coordinates)

mantel.partial(
  bray_dist,
  env_dist,
  geo_dist,
  permutations=999
)
#species差异是否因环境距离控制空间距离
partial_mantel <- mantel.partial(
  
  xdis = bray_dist,
  
  ydis = env_dist,
  
  zdis = geo_dist,
  
  method="pearson",
  
  permutations=999
  
)


print(partial_mantel)

#环境变量湿度的贡献
humidity_dist <- dist(
  scale(
    metadata2$`Humidity(%)`
  )
)


mantel.partial(
  
  bray_dist,
  
  humidity_dist,
  
  geo_dist,
  
  permutations=999
  
)
#环境变量温度的贡献
temperature_dist <- dist(
  scale(
    metadata2$`Temperature(℃)`
  )
)


mantel.partial(
  
  bray_dist,
  
  temperature_dist,
  
  geo_dist,
  
  permutations=999
  
)

#环境变量海拔的贡献
altitude_dist <- dist(
  scale(
    metadata2$`Altitude (m)`
  )
)


mantel.partial(
  
  bray_dist,
  
  altitude_dist,
  
  geo_dist,
  
  permutations=999
  
)
