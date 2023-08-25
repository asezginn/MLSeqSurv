t_harmonised <- read.csv("R/ACC_ProteinCoding.csv",header=TRUE,sep=",")
# remove 0's and NA's from time
# remove NA's from status
t_harmonised <- t_harmonised[!is.na(t_harmonised$time),] # removes NA rows from time
t_harmonised <- t_harmonised[!is.na(t_harmonised$status),] # removes NA rows from status
t_harmonised <- t_harmonised[t_harmonised$time != 0,] # removes rows where time is 0


nrow(t_harmonised)
table(t_harmonised$status)

rownames(t_harmonised) <- t_harmonised[,1]
t_harmonised2 <- t_harmonised[,-1]

# t_harmonised2 <- t_harmonised # Bazı verilerde X sütunu yok

names(t_harmonised2)[names(t_harmonised2) == "t_harmonised$time"] <- "time"
names(t_harmonised2)[names(t_harmonised2) == "t_harmonised$status"] <- "status"
mydata <- t_harmonised2
mydata2 <- mydata %>% select(-time, -status)

satir <- cbind(rownames(mydata),c(1:nrow(mydata)), mydata$time, mydata$status)


tdatam <- TaskSurv$new("mydata", mydata, time = "time", event = "status")
task <- tdatam


split = partition(task,ratio = 0.70)
# roughly same distribution of the target label
prop.table(table(task$truth()))
prop.table(table(task$truth(split$train)))
prop.table(table(task$truth(split$test)))



##### train data #####
train_ind <- as.vector(split$train)
train_data <- mydata[train_ind,]
fwrite(train_data, file = paste("R/trainset_ACC",1,".csv",sep = ""))# replace 1 with the index of the for loop

##### test data ######
test_ind <- as.vector(split$test)
test_data <- mydata[test_ind,]
fwrite(test_data, file = paste("R/testset_ACC",1,".csv",sep = ""))# replace 1 with the index of the for loop


#### control #####
table(train_data$status)
table(test_data$status)


# Start from here
library(mlr3)
library(mlr3proba)
library(mlr3verse)
library(mlr3tuningspaces)
library(reticulate)
library(survivalmodels)
library(dplyr)
library(survival)
library(survminer)
library(plyr)
library(edgeR)
library(DESeq2)
library(data.table)
library(caret)
library(RSBID)
library(dplyr)
library(MASS)
library(ROCR)
library(limma)
library(prioritylasso)
library(ipflasso)
library(Hmisc)
library(ggplot2)
library(gridExtra)

train_s <- read.csv("data/trainset1.csv")
test_s <- read.csv("data/testset1.csv")

cols_t <- as.vector(colnames(train_data))
cols_train <- cols_t[! cols_t %in% c('time', 'status')]
cols_train <- as.vector(cols_train)
train_d <- subset(train_data, select=cols_train)

cols_te <- as.vector(colnames(test_data))
cols_test <- cols_te[! cols_te %in% c('time', 'status')]
cols_test <- as.vector(cols_test)
test_d <- subset(test_data, select=cols_test)

train_matrix <- as.matrix(t(train_d))
test_matrix <- as.matrix(t(test_d))

# call constructor
data_obj <- MLSeqSurv(train_s, test_s)

install_pycox(
  method = "auto",
  conda = "auto",
  pip = TRUE,
  install_torch = TRUE
)
#install_keras(pip = TRUE, install_tensorflow = TRUE)
#install.packages("devtools")
# devtools::install_github("rstudio/keras")
# devtools::install_github("rstudio/tensorflow")
library(keras)
library(tensorflow)
library(pseudo)

# test the functions
result <- survival(data = data_obj, method = "blackboost",
              fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
              trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), mstop = 100)

result <- survival(data = data_obj, method = "cforest",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), ntree = 100)

result <- survival(data = data_obj, method = "coxboost",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")))

result <- survival(data = data_obj, method = "ctree",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")))

result <- survival(data = data_obj, method = "gbm",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), bag.fraction = 0.9)

result <- survival(data = data_obj, method = "glmboost",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")))

result <- survival(data = data_obj, method = "lasso",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), alpha = 1,s = 0.01)

result <- survival(data = data_obj, method = "ridge",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), alpha = 0,s = 0.01)

result <- survival(data = data_obj, method = "elasticnet",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), alpha = 0.5,s = 0.01)

result <- survival(data = data_obj, method = "obliqueRSF",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")))

result <- survival(data = data_obj, method = "penalized",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), lambda1=10, lambda2=10)

result <- survival(data = data_obj, method = "ranger",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")))



result <- survival(data = data_obj, method = "rfsrc",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               tuneGrid = param_list)




result <- survival(data = data_obj, method = "rpart",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")))

result <- survival(data = data_obj, method = "svm",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), type = "hybrid", diff.meth = "makediff3", kernel = "lin_kernel", gamma.mu = c(100,1000))

result <- survival(data = data_obj, method = "xgboost_gbtree",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), booster = "gbtree")

result <- survival(data = data_obj, method = "xgboost_gblinear",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), booster = "gblinear")

result <- survival(data = data_obj, method = "xgboost_dart",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), booster = "dart")

result <- survival(data = data_obj, method = "gamboost",
                    fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
                    trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), center = TRUE, baselearner = "bols", family = "cindex")

result <- survival(data = data_obj, method = "mboost",
                    fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
                    trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), baselearner = "bols", center = TRUE, family = "cindex")





result <- survival(data = data_obj, method = "coxtime",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), optimizer = "sgd", frac = 0.3, epochs = 20)
result <- survival(data = data_obj, method = "deephit",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), optimizer = "sgd", frac = 0.3, epochs = 20)

result <- survival(data = data_obj, method = "deepsurv",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), optimizer = "sgd", frac = 0.3, epochs = 20)


result <- survival(data = data_obj, method = "loghaz",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), optimizer = "sgd", frac = 0.3, epochs = 20)

result <- survival(data = data_obj, method = "pchazard",
               fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")),
               trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 1), list("random_search")), optimizer = "sgd", frac = 0.3, epochs = 20)




result_ipf <- survival(data = data_obj, method = "ipflasso", balancer = "ROS", trainParams = list(TRUE, c(1,1), 5, 10, 0))

result_pl <- survival(data = data_obj, method = "prioritylasso", balancer = "ROS", trainParams = list(TRUE, "lambda.min", TRUE, 5, TRUE, 10) )

result_rfsrc <- survival(data = data_obj, method = "rfsrc",
                         fsParams = list(list("repeated_cv", repeats = 1, folds = 2), list("surv.cindex"), list("evals", n_evals = 3), list("random_search")),
                         trainParams = list(list("repeated_cv", repeats = 2, folds = 2), list("surv.cindex"), list("evals", n_evals = 3), list("random_search")),
                         tuneGrid = param_list)

resultIPF <- survPredict(result_ipf)

resultPL <- survPredict(result_pl)

resultMLR <- survPredict(result_rfsrc)


survPlot(resultIPF)

survPlot(resultPL)

survPlot(resultMLR)
# survivala test verisi gönderilemeyecek
# fonksiyonda train verisi üzerinde predict yapılacak (weight sütunu çıkarılcak)
# model c-indexi ile birleştirilip return edilecek.
# test verisi ile model yeni bir fonksiyona gönderilecek
# yeni fonksşyon


# might be better if we simply return the model instead of c-index
# this allows us to feed model into another predict function for survival hazard preds
# as well as a function that returns the c-index on train/test data
row1 <- c("numeric", "nu", 0, 1)
test <- data.frame(type = row1[1], paramid = row1[2], lower = row1[3], upper = row1[4])

test <- rbind(test, row1)
test <- rbind(test, row1)
test <- rbind(test, row1)

test$levels <- list(NA, c("coxph","weibull"), NA, NA)
test$type <- list("numeric", "character", "integer", "logical")
test$paramid <- list("mtry.ratio","splitrule","ntree","membership")
test$lower <- list(0, NA, 500, NA)
test$upper <- list(0.7, NA, 1500, NA)
test$levels <- list(NA, c("logrank","bs.gradient"), NA, NA)
# test$default <- list(NA, "logrank", 1000, TRUE)


param_list <- list()


i <- 1
for (i in 1:nrow(test)){
  row <- test[i,]
  if (row$type == "numeric"){
    param_list <- c(param_list, ParamDbl$new(row$paramid[[1]], lower = row$lower[[1]], upper = row$upper[[1]]))
  }
  else if (row$type == "character"){
    param_list <- c(param_list, ParamFct$new(row$paramid[[1]], levels = row$levels[[1]]))
  }
  else if (row$type == "integer"){
    param_list <- c(param_list, ParamInt$new(row$paramid[[1]], lower = row$lower[[1]], upper = row$upper[[1]]))
  }
  else if (row$type == "logical"){
    param_list <- c(param_list, ParamLgl$new(row$paramid[[1]]))
  }
  else{
    print("This parameter type is not supported")
  }

}

param_list

# ParamLgl$new(row$paramid[[1]], )

tune_test<- ParamSet$new(param_list)





get_surv_preds_IPF_bin <- function(t) {
  if (sum(trunc_time_grid <= t) != 0) {
    final_index <- max(which(trunc_time_grid <= t))
    haz <- as.matrix(hazard_preds_IPF_bin[, 1:final_index])
    anti_haz <- 1 - haz
    surv <- apply(anti_haz, MARGIN = 1, prod)
  }
  else {
    surv <- rep(1, nrow(hazard_preds_IPF_bin))
  }
  return(surv)
}

surv_preds_IPF_bin <- apply(X = matrix(newtimes), FUN = get_surv_preds_IPF_bin,
                            MARGIN = 1)

interval_75 <- findInterval(quantile(resultMLR@times)[4], resultMLR@times)

surv_preds_IPF_bin_spectime <- resultMLR@survHazards[,interval_75]

test_res <- rcorr.cens(surv_preds_IPF_bin_spectime, Surv(resultMLR@preprocessed_test$time, resultMLR@preprocessed_test$status))

test_res <- rcorr.cens(resultMLR@survHazards, Surv(resultMLR@preprocessed_test$time, resultMLR@preprocessed_test$status))




