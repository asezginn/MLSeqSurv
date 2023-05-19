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

train_data <- read.csv("R/trainset_ACC1.csv")
test_data <- read.csv("R/testset_ACC1.csv")

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

data_obj <- MLSeqSurv(train_data, test_data)
# test the functions
surv(data_obj, method = "blackboost", preProcessing = "deseq-vst",
     fsParams = list(list("repeated_cv", repeats = 10, folds = 2), list("surv.cindex")),
     atParams = list(list("repeated_cv", repeats = 10, folds = 2), list("surv.cindex"), list("evals", n_evals = 5), list("random_search")),
     paramGrid = NULL)

