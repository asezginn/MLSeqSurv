library(mlr3)
library(mlr3proba)
library(mlr3verse)
library(mlr3tuningspaces)
library(reticulate)
library(survivalmodels)
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(edgeR)
library(DESeq2)
library(foreach)
# library(doMC)
library(data.table)
library(future)
library(future.apply)
library(caret)

# options(future.globals.maxSize= 4294967296)
# options(expressions = 5e5)
# progressr::handlers(global = TRUE)

# path_to_python <- "/home/rstudio/.local/share/r-miniconda/envs/r-reticulate/bin/python"
# use_python(path_to_python)
# install_pycox(
#   method = "auto",
#   conda = "auto",
#   pip = TRUE,
#   install_torch = TRUE
# )
#install_keras(pip = TRUE, install_tensorflow = TRUE)
#install.packages("devtools")
# devtools::install_github("rstudio/keras")
# devtools::install_github("rstudio/tensorflow")
# library(keras)
# library(tensorflow)
# library(pseudo)


# mydata <- TH.data::wpbc
# mydata$status <- revalue(mydata$status, c("R"="1"))
# mydata$status <- revalue(mydata$status, c("N"="0"))
#
# mydata$time <- as.numeric(mydata$time)
# mydata$status <- as.numeric(as.character(mydata$status))
# mydata <- mydata[,-34]

dir.create(path = "Results_CompSurv/ACC_ProteinCoding") # Create necessary folders
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/cindex_time")
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/preFeatSel")
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/nearzero_testset")
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/nearzero_trainset")
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/selectedfeatures")
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/testset")
dir.create(path = "Results_CompSurv/ACC_ProteinCoding/trainset")


t_harmonised <- read.csv("cancerData/ACC_ProteinCoding.csv",header=TRUE,sep=",")
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

# plan(multisession(workers = 50))

for (index in 1:30){



  ######## SPLIT TRAIN/TEST #######
  ## define task
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
  fwrite(train_data, file = paste("Results_CompSurv/ACC_ProteinCoding/trainset/trainset",index,".csv",sep = ""))# replace 1 with the index of the for loop

  ##### test data ######
  test_ind <- as.vector(split$test)
  test_data <- mydata[test_ind,]
  fwrite(test_data, file = paste("Results_CompSurv/ACC_ProteinCoding/testset/testset",index,".csv",sep = ""))# replace 1 with the index of the for loop


  #### control #####
  table(train_data$status)
  table(test_data$status)

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



  ######### PRE-PROCESS ########
  ###### near-zero variance filtering ####
  dds_tr <- DESeqDataSetFromMatrix(countData = train_matrix,
                                   colData = as.data.frame(colnames(train_matrix)),
                                   design = ~ 1)
  keep_tr <- rowSums(counts(dds_tr)) >= 10
  dds_tr <- dds_tr[keep_tr,]
  dds_tr
  fwrite(as.data.frame(assay(dds_tr)), file = paste("Results_CompSurv/ACC_ProteinCoding/nearzero_trainset/nearzero_trainset",index,".csv",sep = ""))# replace 1 with the index of the for loop


  dds_ts <- DESeqDataSetFromMatrix(countData = test_matrix,
                                   colData = as.data.frame(colnames(test_matrix)),
                                   design = ~ 1)
  dds_ts <- dds_ts[keep_tr,]
  fwrite(as.data.frame(assay(dds_ts)), file = paste("Results_CompSurv/ACC_ProteinCoding/nearzero_testset/nearzero_testset",index,".csv",sep = ""))# replace 1 with the index of the for loop

  #### train set ####
  normalization <- "deseq"
  transformation <- "vst"

  dataSF_tr <- estimateSizeFactors(dds_tr)    # Estimate Size factors:
  dataDisp_tr <- estimateDispersions(dataSF_tr, fitType = "local")  # Estimate dispersions:
  transformedData_tr <- varianceStabilizingTransformation(dataDisp_tr, fitType = "local")

  dataVST_tr <- t(as.matrix(assay(transformedData_tr)))
  input_tr <- dataVST_tr   ## Input data from transformed expression data.

  trainParameters_tr <- list(disperFunc = dispersionFunction(dataDisp_tr),
                             sizeFactors = sizeFactors(dataDisp_tr))


  #### test set ####
  #Calculation of test set size factors using geometric means from training data
  #Genes in the row, samples in the column
  geomts = test_matrix / exp(rowMeans(log(train_matrix+0.5)))   ## Geometric mean of test data using estimators from train data
  sizeF.ts = apply(geomts, 2, function(x) median(x))

  test.dataSF <- estimateSizeFactors(dds_ts)    # Estimate Size factors:
  sizeFactors(test.dataSF) <- sizeF.ts             # Replace size factors with size factors which are estimates using train set parameters.

  ## Change dispersion function of test data with dispersion function of train data
  dispersionFunction(test.dataSF) <- trainParameters_tr$disperFunc

  transformedData <- varianceStabilizingTransformation(test.dataSF, fitType = "local", blind = FALSE)

  dataVST <- t(as.matrix(assay(transformedData)))
  input_ts <- dataVST   ## Input data from transformed expression data.

  input_tr_nzv <- nearZeroVar(input_tr, saveMetrics = FALSE)
  input_tr <- input_tr[,-as.vector(input_tr_nzv)]
  input_ts <- input_ts[,-as.vector(input_tr_nzv)]


  ############################### VARYANS-FILTERING INDEX  ############################################
  ss = apply(t(input_tr),1,sd)
  ort = apply(t(input_tr),1,mean)
  cv = ss/ort
  siralama = order(cv, decreasing = TRUE)
  cv_train <- t(input_tr)[siralama[1:1900],]
  cv_test <- t(input_ts)[siralama[1:1900],]

  tvsd_train2 <- as.data.frame(cbind(train_ind, satir[as.vector(train_ind),][,3:4], t(cv_train)))
  tvsd_test2 <- as.data.frame(cbind(test_ind, satir[as.vector(test_ind),][,3:4], t(cv_test)))
  #comp <- as.data.frame(rbind(tvsd_train2,tvsd_test2))
  names(tvsd_train2)[names(tvsd_train2) == "V2"] <- "time"
  names(tvsd_train2)[names(tvsd_train2) == "V3"] <- "status"
  names(tvsd_train2)[names(tvsd_train2) == "train_ind"] <- "sirano"
  tvsd_train2$time <- as.integer(tvsd_train2$time)
  tvsd_train2$status <- as.integer(tvsd_train2$status)
  tvsd_train2$sirano <- as.integer(tvsd_train2$sirano)
  # tvsd_train2$sirano <- as.integer(tvsd_train2$sirano)
  tvsd_train2[] <- apply(tvsd_train2, 2, function(x) as.double((x)))

  fwrite(tvsd_train2, file = paste("Results_CompSurv/ACC_ProteinCoding/preFeatSel/train_set", index, ".csv", sep = ""))


  names(tvsd_test2)[names(tvsd_test2) == "V2"] <- "time"
  names(tvsd_test2)[names(tvsd_test2) == "V3"] <- "status"
  names(tvsd_test2)[names(tvsd_test2) == "test_ind"] <- "sirano"
  tvsd_test2$time <- as.integer(tvsd_test2$time)
  tvsd_test2$status <- as.integer(tvsd_test2$status)
  tvsd_test2$sirano <- as.integer(tvsd_test2$sirano)
  tvsd_test2[] <- apply(tvsd_test2, 2, function(x) as.double((x)))

  fwrite(tvsd_test2, file = paste("Results_CompSurv/ACC_ProteinCoding/preFeatSel/test_set", index, ".csv", sep = "") )


}
