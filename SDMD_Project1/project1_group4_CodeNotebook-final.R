## Load Packages
suppressPackageStartupMessages({
    library(Boruta)
    library(caret)
    library(class)
    library(corrplot)
    library(data.table)
    library(DataExplorer)
    library(dslabs)
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(gtsummary)
    library(kernlab)
    library(kknn)
    library(MASS)
    library(mlbench)
    library(mlr)
    library(mRMRe)
    library(naivebayes)
    library(naniar)
    library(PRROC)
    library(purrr)
    library(readr)
    library(reshape2)
    library(ROSE)
    library(skimr)
    library(splitstackshape)
    library(superml)
    library(tidyverse)
    library(VIM)
    })

## Read original dataset
df <- read.csv("dataset_original.csv", header = TRUE)

## Replace all missing values (?) with NA
df[df == "?"] <- NA

## Change all cols type to numeric
df[] <- lapply(df, as.numeric)

## Create dataset only with the relevant INPUT (Except response variable) cols for our problem
### Input cols: 2-112 except 95, 102, 105 
df_inputs <- df[-c(1, 95, 102, 105, 113:124)]

head(df_inputs)

## % of Missing Values
df.nas <- sum(is.na(df_inputs))

df.size <- nrow(df_inputs) * ncol(df_inputs)

pct_nas <- (df.nas * 100) / df.size
cat("Total percetange of missing values: ", pct_nas)

## Which cols have NA?
NAs = names(df_inputs)[sapply(df_inputs, anyNA)]
cat("Number of cols with NAs:", length(NAs), "\n\n")

cat("Names of cols with missing values: ", "\n")
print(NAs)

## Create dataset with only the vars with more than 25% of NAs
pct_nas <- 0.25
na_cols_drop <- c()

for (i in 1:length(df_inputs)){
    if (sum(is.na(df_inputs[, i])) < (1700 * pct_nas)) {
        na_cols_drop <- append(na_cols_drop, i)
    }
}


df_nas_25 <- df_inputs[-na_cols_drop]

## Plot % of NAs per var
plot_w <- 10
plot_h <- 4
options(repr.plot.width=plot_w, repr.plot.height=plot_h, unit="cm")

gg_miss_var(df_nas_25,  show_pct=TRUE) + labs(y="Missing Values [%]", x=" ") + 
    theme(text=element_text(size=14))

## Names of columns with more 80% or more of missing values
na_over_80 <- function(df) {
    names(df)[((sapply(df, function(x) mean(is.na(x)))) * 100) >= 80]
}

                       
drop <- na_over_80(df_inputs)
cat("Cols with more tha 80% of missing values: ", "\n")
print(drop)

## Remove columns with over 80% of missing values
df_inputs <- df_inputs[ , !(names(df_inputs) %in% drop)]
head(df_inputs)

imputed <- kNN(df_inputs, k=9)
df_imputed <- imputed[1:106]
head(df_imputed)

## Transform ordinal variables into binary variables, by doing dummies
ordinalvars = c("INF_ANAM","STENOK_AN","FK_STENOK","IBS_POST","GB","DLIT_AG", "ZSN_A",
             "ant_im","lat_im","inf_im","post_im","TIME_B_S","R_AB_1_n","R_AB_2_n",
             "NA_R_1_n","NA_R_2_n","NOT_NA_1_n","NOT_NA_2_n")

df_new <- df_imputed
for (i in ordinalvars){
  df_new[,i] = as.factor(df_imputed[,i])
}
  
df_dummies <- createDummyFeatures(df_new, method = "reference")

## Standardize numeric variables
df_dummies$AGE = scale(df_dummies$AGE)
df_dummies$S_AD_KBRIG = scale(df_dummies$S_AD_KBRIG)
df_dummies$D_AD_KBRIG = scale(df_dummies$D_AD_KBRIG)
df_dummies$S_AD_ORIT = scale(df_dummies$S_AD_ORIT)
df_dummies$D_AD_ORIT = scale(df_dummies$D_AD_ORIT)
df_dummies$Na_BLOOD = scale(df_dummies$Na_BLOOD)
df_dummies$ALT_BLOOD = scale(df_dummies$ALT_BLOOD)
df_dummies$AST_BLOOD = scale(df_dummies$AST_BLOOD)
df_dummies$K_BLOOD = scale(df_dummies$K_BLOOD)
df_dummies$L_BLOOD = scale(df_dummies$L_BLOOD)
df_dummies$ROE = scale(df_dummies$ROE)

head(df_dummies)

cormatrix = round(cor(df_dummies, method="pearson"),2)
head(cormatrix)

highlyCorrelated <- findCorrelation(cormatrix, cutoff=0.5)

cat("Number of featutes selected by the Pearson Correlation Method: ", length(colnames(df_dummies)[-highlyCorrelated]), "\n\n")
cat("Names of featutes selected by the Pearson Correlation Method: ", "\n")
print(colnames(df_dummies)[-highlyCorrelated])

for (i in 1:106){
   df_dummies[,i] <- as.numeric(as.character(df_dummies[,i]))
   }
dd = mRMR.data(data = df_dummies)

#Choose to retain 25% of the features
mrmr = mRMR.classic(data = dd, target_indices = c(1), feature_count = 40)
feature_index = mrmr@filters$'1'

mrmr_selectVars = colnames(df_dummies)[feature_index]
df_mrmr = df_dummies[,mrmr_selectVars]

cat("Number of featutes selected by the mRMR Method: ", length(mrmr_selectVars), "\n\n")
cat("Names of featutes selected by the mRMR Method: ", "\n")
print(mrmr_selectVars)

set.seed(852456)
boruta_output <- Boruta(x=df_dummies, y=df$FIBR_JELU)
# names(boruta_output)

roughFixMod <- TentativeRoughFix(boruta_output)
boruta_selectedVars <- getSelectedAttributes(roughFixMod)

imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2[order(-imps2$meanImp), ]  # descending sort

boruta_data = df_dummies[,boruta_selectedVars]

plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")

cat("Number of featutes selected by the Boruta Method: ", length(boruta_selectedVars), "\n\n")
cat("Names of featutes selected by the Boruta Method: ", "\n")
print(boruta_selectedVars)

defaultW <- getOption("warn")
options(warn = -1)

subsets <- c(10, 20, 30, 40, 50, 75)
ctrl <- rfeControl(functions = lmFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(df_dummies, df$FIBR_JELU,
                sizes = subsets,
                rfeControl = ctrl)
lmProfile

RFE_selectVars = predictors(lmProfile)
rfe_data = df_dummies[,RFE_selectVars]

plot(lmProfile, type = c("g", "o"))

cat("Number of featutes selected by the RFE Method: ", length(RFE_selectVars), "\n\n")
cat("Names of featutes selected by the RFE Method: ", "\n")
print(RFE_selectVars)

## Add response col to the dataset
df_final <- df_imputed
df_final$FIBR_JELU <- df$FIBR_JELU

df_final["FIBR_JELU"][df_final["FIBR_JELU"] == 0] <- "0"
df_final["FIBR_JELU"][df_final["FIBR_JELU"] == 1] <- "1"


df_final_dummies <- df_dummies
df_final_dummies$FIBR_JELU <- df$FIBR_JELU

## Create train (80%) and test (20%) sets
set.seed(12345)
df_split <- stratified(df_final, c('FIBR_JELU'), 0.8, bothSets=TRUE)

train <- df_split$SAMP1
train <- as.data.frame(train) # convert data.table to data.frame

test <- df_split$SAMP2
test <- as.data.frame(test) # convert data.table to data.frame

## Create train (80%) and test (20%) sets
set.seed(12345)
df_split <- stratified(df_final_dummies, c('FIBR_JELU'), 0.8, bothSets=TRUE)

train_dummies <- df_split$SAMP1
train_dummies <- as.data.frame(train_dummies) # convert data.table to data.frame

test_dummies <- df_split$SAMP2
test_dummies <- as.data.frame(test_dummies) # convert data.table to data.frame

## Sanity Check!! Must be the same in the test and train vs test_dummies and train_dummies

### test and train --------------------------------
cat("test and train:")
#### Number of elements per class in train set
table(train$FIBR_JELU)
#### Number of elements per class in test set
table(test$FIBR_JELU)

cat("\n\n")

### test_dummies and train_dummies ----------------
cat("test_dummies and train_dummies:")
#### Number of elements per class in train set
table(train_dummies$FIBR_JELU)
#### Number of elements per class in test set
table(test_dummies$FIBR_JELU)

## Create a new balanced train set, with the same size of the original train set, using both under- and over sampling
train.balanced <- ovun.sample(FIBR_JELU~., train, method="both", p=0.30, seed=54321)$data

train_dummies.balanced <- ovun.sample(FIBR_JELU~., train_dummies, method="both", p=0.30, seed=54321)$data

## Sanity Check!! Must be the same in the train.balanced vs train_dummies.balanced

### train.balanced --------------------------------
cat("train.balanced:")
#### Number of elements per class in train set
table(train.balanced$FIBR_JELU)
#### % of the minority class in new train set
(length(train.balanced["FIBR_JELU"][train.balanced["FIBR_JELU"] == 1]) / nrow(train.balanced)) * 100

cat("\n\n")

### train_dummies.balanced ------------------------
cat("train_dummies.balanced:")
#### Number of elements per class in train set
table(train_dummies.balanced$FIBR_JELU)
#### % of the minority class in new train set
(length(train_dummies.balanced["FIBR_JELU"][train_dummies.balanced["FIBR_JELU"] == 1]) / nrow(train_dummies.balanced)) * 100

y.train <- train.balanced$FIBR_JELU
y.test <- test$FIBR_JELU

y.train.factor <- factor(train.balanced$FIBR_JELU)
y.test.factor <- factor(test$FIBR_JELU)

x.train.all <- train.balanced[ , !(names(train.balanced) %in% c("FIBR_JELU"))]
x.test.all <- test[ , !(names(test) %in% c("FIBR_JELU"))]

mrmr_signif <- c("NITR_S", "GT_POST", "zab_leg_01", "INF_ANAM", "nr04", "DLIT_AG", "lat_im", "TIKL_S_n", 
                 "IBS_POST", "inf_im", "n_p_ecg_p_10", "ZSN_A", "n_p_ecg_p_09", "fibr_ter_07", "NA_KB", 
                 "fibr_ter_06", "n_p_ecg_p_03", "zab_leg_03", "ritm_ecg_p_01", "GB", "O_L_POST", "TIME_B_S", 
                 "zab_leg_02", "S_AD_KBRIG", "endocr_01", "GEPAR_S_n", "n_r_ecg_p_01", "STENOK_AN", "ALT_BLOOD", 
                 "MP_TP_POST", "fibr_ter_03", "ROE", "B_BLOK_S_n", "FK_STENOK", "SEX")

x.train.mrmr <- train.balanced[mrmr_signif]
x.test.mrmr <- test[mrmr_signif]

mrmr_signif_dummies <- c("NITR_S", "GT_POST", "zab_leg_01", "INF_ANAM.2" , "nr04", "DLIT_AG.1", "lat_im.4", 
                         "TIKL_S_n", "IBS_POST.1", "inf_im.4", "n_p_ecg_p_10", "ZSN_A.1", "n_p_ecg_p_09", 
                         "ZSN_A.3", "fibr_ter_07", "NA_KB", "GB.2", "lat_im.3", "fibr_ter_06", "n_p_ecg_p_03", 
                         "zab_leg_03", "ritm_ecg_p_01", "GB.1", "O_L_POST", "TIME_B_S.5", "zab_leg_02", "S_AD_KBRIG", 
                         "endocr_01", "GEPAR_S_n", "n_r_ecg_p_01", "STENOK_AN.6", "ALT_BLOOD", "GB.3", "MP_TP_POST", 
                         "fibr_ter_03", "ROE", "DLIT_AG.7", "B_BLOK_S_n", "FK_STENOK.2", "SEX")

x.train.mrmr.dummies <- train_dummies.balanced[mrmr_signif_dummies]
x.test.mrmr.dummies <- test_dummies[mrmr_signif_dummies]

boruta_signif <- c("S_AD_KBRIG", "D_AD_KBRIG", "S_AD_ORIT", "D_AD_ORIT", "K_SH_POST", "n_r_ecg_p_10", 
                   "K_BLOOD", "LID_S_n", "ASP_S_n", "GB", "ZSN_A", "R_AB_2_n", "NA_R_2_n")

x.train.boruta <- train.balanced[boruta_signif]
x.test.boruta <- test[boruta_signif]

boruta_signif_dummies <- c("S_AD_KBRIG", "D_AD_KBRIG", "S_AD_ORIT", "D_AD_ORIT", "K_SH_POST", "n_r_ecg_p_10", 
                           "K_BLOOD", "LID_S_n", "ASP_S_n", "GB.3", "ZSN_A.2", "R_AB_2_n.2", "NA_R_2_n.2")

x.train.boruta.dummies <- train_dummies.balanced[boruta_signif_dummies]
x.test.boruta.dummies <- test_dummies[boruta_signif_dummies]

rfe_signif <- c("n_r_ecg_p_10","R_AB_2_n", "n_r_ecg_p_08", "np01", "NOT_NA_2_n", "nr07", "SVT_POST",
                "fibr_ter_05", "GB", "np09", "NA_R_2_n", "n_p_ecg_p_08", "fibr_ter_06", "n_p_ecg_p_03",
                "ritm_ecg_p_04", "K_SH_POST", "fibr_ter_08", "ZSN_A", "fibr_ter_07", "FIB_G_POST",
                "n_r_ecg_p_09", "nr02", "np07", "np05", "STENOK_AN", "ritm_ecg_p_08", "n_r_ecg_p_06",
                "TIKL_S_n", "n_p_ecg_p_04", "LID_S_n", "FK_STENOK", "GT_POST", "n_p_ecg_p_01", "fibr_ter_02")

x.train.rfe <- train.balanced[rfe_signif]
x.test.rfe <- test[rfe_signif]

rfe_signif_dummies <- c("n_r_ecg_p_10", "R_AB_2_n.3", "n_r_ecg_p_08", "np01", "NOT_NA_2_n.3", "nr07" , 
                        "SVT_POST", "fibr_ter_05", "GB.1", "np09", "NA_R_2_n.2", "n_p_ecg_p_08", "fibr_ter_06", 
                        "n_p_ecg_p_03", "ritm_ecg_p_04", "K_SH_POST", "fibr_ter_08", "ZSN_A.4","fibr_ter_07",
                        "ZSN_A.2", "FIB_G_POST","n_r_ecg_p_09", "nr02", "np07", "np05", "STENOK_AN.2", 
                        "ritm_ecg_p_08", "n_r_ecg_p_06", "R_AB_2_n.2", "TIKL_S_n", "n_p_ecg_p_04", "STENOK_AN.4", 
                        "STENOK_AN.1", "LID_S_n" , "STENOK_AN.3", "FK_STENOK.1", "STENOK_AN.6", "GT_POST", 
                        "n_p_ecg_p_01", "fibr_ter_02")

x.train.rfe.dummies <- train_dummies.balanced[rfe_signif_dummies]
x.test.rfe.dummies <- test_dummies[rfe_signif_dummies]

## Cross Validation with 5 folds
fitControl_knn_mrmr <- trainControl(method = "cv", number = 5)

## Train
model_knn_mrmr <- caret::train(x = x.train.mrmr.dummies, y = y.train.factor, 
                               method = "knn", 
                               trControl = fitControl_knn_mrmr,
                               tuneLength = 10)
## Classify
y.pred.knn_mrmr <- predict(model_knn_mrmr, newdata = x.test.mrmr.dummies)

## Plot Cross Validation
plot(model_knn_mrmr)

## Performance
knn_mrmr_perf <- confusionMatrix(data=y.pred.knn_mrmr, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(knn_mrmr_perf)

### Precision
knn_mrmr_perf$byClass["Precision"]

### Sensitivity (or Recall)
knn_mrmr_perf$byClass["Sensitivity"]

### Specificity
knn_mrmr_perf$byClass["Specificity"]

### Accurancy
knn_mrmr_perf$overall["Accuracy"]

### Balanced Accurancy
knn_mrmr_perf$byClass["Balanced Accuracy"]

## Cross Validation with 5 folds
fitControl_knn_boruta <- trainControl(method = "cv", number = 5)

## Train
model_knn_boruta <- caret::train(x = x.train.boruta.dummies, y = y.train.factor, 
                               method = "knn", 
                               trControl = fitControl_knn_boruta,
                               tuneLength = 10)
## Classify
y.pred.knn_boruta <- predict(model_knn_boruta, newdata = x.test.boruta.dummies)

## Plot Cross Validation
plot(model_knn_boruta)

## Performance
knn_boruta_perf <- confusionMatrix(data=y.pred.knn_boruta, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(knn_boruta_perf)

### Precision
knn_boruta_perf$byClass["Precision"]

### Sensitivity (or Recall)
knn_boruta_perf$byClass["Sensitivity"]

### Specificity
knn_boruta_perf$byClass["Specificity"]

### Accurancy
knn_boruta_perf$overall["Accuracy"]

### Balanced Accurancy
knn_boruta_perf$byClass["Balanced Accuracy"]

## Cross Validation with 5 folds
fitControl_knn_rfe <- trainControl(method = "cv", number = 5)

## Train
model_knn_rfe <- caret::train(x = x.train.rfe.dummies, y = y.train.factor, 
                               method = "knn", 
                               trControl = fitControl_knn_rfe,
                               tuneLength = 10)
## Classify
y.pred.knn_rfe <- predict(model_knn_rfe, newdata = x.test.rfe.dummies)

## Plot Cross Validation
plot(model_knn_rfe)

## Performance
knn_rfe_perf <- confusionMatrix(data=y.pred.knn_rfe, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(knn_rfe_perf)

### Precision
knn_rfe_perf$byClass["Precision"]

### Sensitivity (or Recall)
knn_rfe_perf$byClass["Sensitivity"]

### Specificity
knn_rfe_perf$byClass["Specificity"]

### Accurancy
knn_rfe_perf$overall["Accuracy"]

### Balanced Accurancy
knn_rfe_perf$byClass["Balanced Accuracy"]

## Train
NBay_all <- gaussian_naive_bayes(data.matrix(x.train.all), y.train)

## Classify
y.pred.NB_all <- predict(NBay_all, data.matrix(x.test.all))

## Performance
NB_all_perf <- confusionMatrix(data=y.pred.NB_all, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(NB_all_perf)

### Precision
NB_all_perf$byClass["Precision"]

### Sensitivity (or Recall)
NB_all_perf$byClass["Sensitivity"]
NB_all_perf$byClass["Recall"]

### Specificity
NB_all_perf$byClass["Specificity"]

### Accurancy
NB_all_perf$overall["Accuracy"]

### Balanced Accurancy
NB_all_perf$byClass["Balanced Accuracy"]

## Train
NBay_mrmr <- gaussian_naive_bayes(data.matrix(x.train.mrmr), y.train)

## Classify
y.pred.NB_mrmr <- predict(NBay_mrmr, data.matrix(x.test.mrmr))

## Performance
NB_mrmr_perf <- confusionMatrix(data=y.pred.NB_mrmr, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(NB_mrmr_perf)

### Precision
NB_mrmr_perf$byClass["Precision"]

### Sensitivity (or Recall)
NB_mrmr_perf$byClass["Sensitivity"]

### Specificity
NB_mrmr_perf$byClass["Specificity"]

### Accurancy
NB_mrmr_perf$overall["Accuracy"]

### Balanced Accurancy
NB_mrmr_perf$byClass["Balanced Accuracy"]

## Train
NBay_boruta <- gaussian_naive_bayes(data.matrix(x.train.boruta), y.train)

## Classify
y.pred.NB_boruta <- predict(NBay_boruta, data.matrix(x.test.boruta))

## Performance
NB_boruta_perf <- confusionMatrix(data=y.pred.NB_boruta, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(NB_boruta_perf)

### Precision
NB_boruta_perf$byClass["Precision"]

### Sensitivity (or Recall)
NB_boruta_perf$byClass["Sensitivity"]

### Specificity
NB_boruta_perf$byClass["Specificity"]

### Accurancy
NB_boruta_perf$overall["Accuracy"]

### Balanced Accurancy
NB_boruta_perf$byClass["Balanced Accuracy"]

## Train
NBay_rfe <- gaussian_naive_bayes(data.matrix(x.train.rfe), y.train)

## Classify
y.pred.NB_rfe <- predict(NBay_rfe, data.matrix(x.test.rfe))

## Performance
NB_rfe_perf <- confusionMatrix(data=y.pred.NB_rfe, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(NB_rfe_perf)

### Precision
NB_rfe_perf$byClass["Precision"]

### Sensitivity (or Recall)
NB_rfe_perf$byClass["Sensitivity"]

### Specificity
NB_rfe_perf$byClass["Specificity"]

### Accurancy
NB_rfe_perf$overall["Accuracy"]

### Balanced Accurancy
NB_rfe_perf$byClass["Balanced Accuracy"]

## Cross Validation with 5 folds
fitControl_svm_mrmr <- trainControl(method = "cv", number = 5)

## Train
model_svm_mrmr <- caret::train(x = x.train.mrmr.dummies, y = y.train.factor, 
                               method = 'svmLinear', 
                               trControl = fitControl_svm_mrmr,
                               tuneLength = 10)
## Classify
y.pred.svm_mrmr <- predict(model_svm_mrmr, newdata = x.test.mrmr.dummies)

## Performance
svm_mrmr_perf <- confusionMatrix(data=y.pred.svm_mrmr, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(svm_mrmr_perf)

### Precision
svm_mrmr_perf$byClass["Precision"]

### Sensitivity (or Recall)
svm_mrmr_perf$byClass["Sensitivity"]

### Specificity
svm_mrmr_perf$byClass["Specificity"]

### Accurancy
svm_mrmr_perf$overall["Accuracy"]

### Balanced Accurancy
svm_mrmr_perf$byClass["Balanced Accuracy"]

## Cross Validation with 5 folds
fitControl_svm_boruta <- trainControl(method = "cv", number = 5)

## Train
model_svm_boruta <- caret::train(x = x.train.boruta.dummies, y = y.train.factor, 
                               method = 'svmLinear', 
                               trControl = fitControl_svm_boruta,
                               tuneLength = 10)
## Classify
y.pred.svm_boruta <- predict(model_svm_boruta, newdata = x.test.boruta.dummies)

## Performance
svm_boruta_perf <- confusionMatrix(data=y.pred.svm_boruta, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(svm_boruta_perf)

### Precision
svm_boruta_perf$byClass["Precision"]

### Sensitivity (or Recall)
svm_boruta_perf$byClass["Sensitivity"]

### Specificity
svm_boruta_perf$byClass["Specificity"]

### Accurancy
svm_boruta_perf$overall["Accuracy"]

### Balanced Accurancy
svm_boruta_perf$byClass["Balanced Accuracy"]

## Cross Validation with 5 folds
fitControl_svm_rfe <- trainControl(method = "cv", number = 5)

## Train
model_svm_rfe <- caret::train(x = x.train.rfe.dummies, y = y.train.factor, 
                              method = 'svmLinear', 
                              trControl = fitControl_svm_rfe,
                              tuneLength = 10,
                              fitBest = FALSE,
                              returnData = TRUE)

## Classify
y.pred.svm_rfe <- predict(model_svm_rfe, newdata = x.test.rfe.dummies)

## Performance
svm_rfe_perf <- confusionMatrix(data=y.pred.svm_rfe, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(svm_rfe_perf)

### Precision
svm_rfe_perf$byClass["Precision"]

### Sensitivity (or Recall)
svm_rfe_perf$byClass["Sensitivity"]

### Specificity
svm_rfe_perf$byClass["Specificity"]

### Accurancy
svm_rfe_perf$overall["Accuracy"]

### Balanced Accurancy
svm_rfe_perf$byClass["Balanced Accuracy"]

## Train 
lda_mrmr <- lda(x = x.train.mrmr.dummies, grouping = y.train) 

## Classify
y.pred.lda_mrmr <- predict(lda_mrmr, x.test.mrmr.dummies)$class

## Performance
lda_mrmr_perf <- confusionMatrix(data=y.pred.lda_mrmr, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(lda_mrmr_perf)

### Precision
lda_mrmr_perf$byClass["Precision"]

### Sensitivity (or Recall)
lda_mrmr_perf$byClass["Sensitivity"]

### Specificity
lda_mrmr_perf$byClass["Specificity"]

### Accurancy
lda_mrmr_perf$overall["Accuracy"]

### Balanced Accurancy
lda_mrmr_perf$byClass["Balanced Accuracy"]

## Train
lda_boruta <- lda(x = x.train.boruta.dummies, grouping = y.train) 

## Classify
y.pred.lda_boruta <- predict(lda_boruta, x.test.boruta.dummies)$class

## Performance
lda_boruta_perf <- confusionMatrix(data=y.pred.lda_boruta, reference=y.test.factor, positive = "1")

### Confusion Matrix
cat("Confusion Matrix: ", "\n")
as.table(lda_boruta_perf)

### Precision
lda_boruta_perf$byClass["Precision"]

### Sensitivity (or Recall)
lda_boruta_perf$byClass["Sensitivity"]

### Specificity
lda_boruta_perf$byClass["Specificity"]

### Accurancy
lda_boruta_perf$overall["Accuracy"]

### Balanced Accurancy
lda_boruta_perf$byClass["Balanced Accuracy"]

## Train 
lda_rfe <- lda(x = x.train.rfe.dummies, grouping = y.train) 

# ## Classify
y.pred.lda_rfe <- predict(lda_rfe, x.test.rfe.dummies)$class
