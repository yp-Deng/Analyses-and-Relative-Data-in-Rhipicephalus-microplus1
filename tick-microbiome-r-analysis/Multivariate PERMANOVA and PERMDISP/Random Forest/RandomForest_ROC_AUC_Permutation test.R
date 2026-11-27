library(readxl)
library(dplyr)
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)

#inport data
family <- read_excel(
  "Family.xlsx"
)


metadata <- read_excel(
  "../metadata.xlsx"
)


head(family)
head(metadata)

#combining Family data with Clade data
data <- family %>%
  
  inner_join(
    metadata %>% 
      select(
        SampleID,
        Clade
      ),
    by="SampleID"
  )

#Check
table(data$Clade)

#Extract relative abundance of Family
X <- data %>%
  
  select(
    -SampleID,
    -Clade
  )

#set as numeric
X <- data.frame(
  lapply(
    X,
    as.numeric
  )
)

#lebeled
Y <- factor(
  data$Clade
)

#RandomForest model
set.seed(123)

rf_model <- randomForest(
  
  x=X,
  
  y=Y,
  
  ntree=1000,
  
  importance=TRUE,
  
  mtry=sqrt(ncol(X))
  
)

print(rf_model)

#select top 20 family with impotance
importance_table <- data.frame(
  
  Family =
    rownames(
      importance(rf_model)
    ),
  
  importance(rf_model)
  
)

importance_table <- importance_table %>%
  
  arrange(
    desc(
      MeanDecreaseAccuracy
    )
  )

top20_family <- importance_table %>%
  
  slice(
    1:20
  )


top20_family

#build model using top 20 family
selected_family <- top20_family$Family
X20 <- X[,selected_family]

#Build model using top 20 family
set.seed(123)


rf_final <- randomForest(
  
  x=X20,
  
  y=Y,
  
  ntree=1000,
  
  importance=TRUE
  
)


print(rf_final)

#10-fold cross validation
set.seed(123)


control <- trainControl(
  
  method="cv",
  
  number=10,
  
  classProbs=TRUE,
  
  savePredictions=TRUE
  
)

rf_cv <- train(
  
  x=X20,
  
  y=Y,
  
  method="rf",
  
  trControl=control,
  
  ntree=1000,
  
  importance=TRUE
  
)

rf_cv

#mix matrix
confusionMatrix(
  
  rf_cv$pred$pred,
  
  rf_cv$pred$obs
  
)

#ROC curves and AUC
pred <- rf_cv$pred
colnames(pred)

#multitypes of AUC
roc_multi <- multiclass.roc(
  
  response = pred$obs,
  
  predictor =
    pred[,c(
      "A",
      "B",
      "C"
    )]
  
)


auc_multi <- auc(
  roc_multi
)


auc_multi


#plot ROC of each Clade
roc_A <- roc(
  
  ifelse(
    pred$obs=="A",
    "A",
    "Other"
  ),
  
  pred$A
  
)


plot(
  roc_A,
  main="ROC - Clade A"
)


auc(roc_A)

roc_B <- roc(
  
  ifelse(
    pred$obs=="B",
    "B",
    "Other"
  ),
  
  pred$B
  
)


plot(
  roc_B,
  main="ROC - Clade B"
)


auc(roc_B)

roc_C <- roc(
  
  ifelse(
    pred$obs=="C",
    "C",
    "Other"
  ),
  
  pred$C
  
)


plot(
  roc_C,
  main="ROC - Clade C"
)


auc(roc_C)

plot(roc_A, col = 'red', print.auc = TRUE)
lines(roc_B, col = 'blue',print.auc = TRUE, print.auc.y = 0.4 )
lines(roc_C, col = 'grey',print.auc = TRUE, print.auc.y = 0.3 )

#compute real accuracy
set.seed(123)


rf_observed <- randomForest(
  
  x=X20,
  
  y=Y,
  
  ntree=1000
  
)


observed_accuracy <-
  
  mean(
    rf_observed$predicted == Y
  )


observed_accuracy

#999 permutation
set.seed(123)


n_perm <- 999


perm_accuracy <- numeric(
  n_perm
)



for(i in 1:n_perm){
  
  
  # 随机打乱标签
  
  Y_perm <- sample(Y)
  
  
  
  rf_perm <- randomForest(
    
    x=X20,
    
    y=Y_perm,
    
    ntree=1000
    
  )
  
  
  
  perm_accuracy[i] <-
    
    mean(
      rf_perm$predicted == Y_perm
    )
  
  
}

#COmpute permutation P value
p_value <- (
  
  sum(
    perm_accuracy >= observed_accuracy
  )+1
  
)/(
  
  n_perm+1
  
)


p_value

#plot permutation test
hist(
  
  perm_accuracy,
  
  breaks=30,
  
  main="Random Forest permutation test",
  
  xlab="Accuracy"
  
)


abline(
  
  v=observed_accuracy,
  
  col="red",
  
  lwd=3
  
)

#compute P value of each Family
real_importance <- importance(
  
  rf_observed,
  
  type=1
  
)[,"MeanDecreaseAccuracy"]

head(real_importance)

#Permutation importance matrix
perm_importance <- matrix(
  
  NA,
  
  nrow=length(real_importance),
  
  ncol=n_perm
  
)


rownames(perm_importance) <-
  
  names(real_importance)

#permutation cycle
for(i in 1:n_perm){
  
  
  Y_perm <- sample(Y)
  
  
  rf_perm <- randomForest(
    
    x=X20,
    
    y=Y_perm,
    
    ntree=1000,
    
    importance=TRUE
    
  )
  
  
  perm_importance[,i] <-
    
    importance(
      rf_perm,
      type=1
    )[,"MeanDecreaseAccuracy"]
  
  
}

#compute family P value
family_p <- apply(
  perm_importance,
  1,
  function(x){
    
    p <- (sum(x >= real_importance) + 1) / (n_perm + 1)
    
    return(p)
    
  }
)


family_p_result <- data.frame(
  
  Family = names(family_p),
  
  P_value = family_p
  
)


family_p_result
