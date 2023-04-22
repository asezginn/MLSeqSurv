
library(dplyr)

method <- "blackboost"

origData <- read.csv("./R/ACC_ProteinCoding.csv")

t_harmonised <- origData
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

View(mydata)

data <- new()
if (method == "blackboost"){
  surv.blackboost(data = data, method = "blackboost", preProcessing = preProcessing, fsParams, atParams, paramGrid, ...) # Fix the parameters that are being sent once the generic functions are finalized.
}


surv.blackboost <- function(data = data, method = "blackboost", preProcessing = preProcessing, fsParams, atParams, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.blackboost", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  # measures = msrs(c("surv.cindex"))


  # fsParams is a list of parameters related to feature selection.
  # First element of this list is notable since it can/will contain multiple elements.
  # May need to look into unlisting them from a list format.
  # We can't simply dedicate more elements of the list because the list will contain varying amount of elements depending on the method chosen.
  instance = FSelectInstanceSingleCrit$new(
    task = task_fs,
    learner = learner,
    resampling = do.call(mlr3::rsmp, fsParams[[1]]),
    measure = do.call(mlr3::msr, fsParams[[2]]),
    terminator = do.call(mlr3::trm, fsParams[[3]])
  )

  fselector = do.call(mlr3::fs, fsParams[[4]])
  fselector$optimize(instance)

  features <- as.vector(instance$result_feature_set)

  data@train <- cbind(data@train[,1], subset(data@train, select=features), data@train$time, data@train$status) # this needs to be cleaned up and be more dynamic
  data@test <- cbind(data@test[,1], subset(data@test, select=features), data@test$time, data@test$status)
  # names of the column need to be fixed
  names(data@train)[names(data@train) == "data@train$time"] <- "time"
  names(data@train)[names(data@train) == "data@train$status"] <- "status"
  names(data@test)[names(data@test) == "data@test$time"] <- "time"
  names(data@test)[names(data@test) == "data@test$status"] <- "status"


  task_tune <- TaskSurv$new("task_tune", data@train, time = "time", event = "status")

  tune_ps = ParamSet$new(paramGrid) # It might be better to add the paramGrid directly into the atParams.

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp(atParams[[1]]),
                   measure=msr(atParams[[2]]),
                   terminator=trm(atParams[[3]]),
                   tuner=tnr(atParams[[4]]),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)

  # Create the desired S4 object by merging the results.



}

t=AutoTuner$new(learner=learners[[j]],
                resampling=do.call(mlr3::rsmp, list("repeated_cv", repeats = 10, folds = 2)),
                measure=msr("surv.cindex"),
                terminator=trm("evals", n_evals = 5),
                tuner=tnr("random_search"),
                search_space = tune_spaces[[j]]
)
rsmpParams <- list("repeated_cv", "repeats = 10", "folds = 2")
testList <- list(rsmpParams, "surv.cindex")

rsmp <- rsmpParams[[1]]

rsmpV <- unlist(rsmpParams) # Handling these parameters is fairly tricky.
rsmpV2 <- split(rsmpV)
testList[[1]]

randomFunc <- function(arg1, ...){

  args <- list(...)

  print(args$arg2)

}

randomFunc(1, arg2 = "a", arg3 = 2, list("a", 5, 10))

# An idea could be that we take the first parameter of resampling parameters as given and convert the rest of them. Not very clean
# Another idea is to take almost all of these parameters from ellipsis
# Or, we could assume that user will input these arguments in the correct order and just feed the elements of the list. Problem is, we don't know the
# length of the list. Thus, this isn't very clean either


# Here is another approach:
# Define your outer function
my_outer_function <- function(param1, param2, param_list) {
  # Check if the param_list is a list
  if (!is.list(param_list)) {
    stop("param_list must be a list")
  }

  # Use do.call() to pass the list as arguments to the inner function
  inner_result <- do.call(my_inner_function, as.list(param_list))

  # Return the result
  return(inner_result)
}

# Define your inner function
my_inner_function <- function(param3, param4) {
  # Use the parameters in your inner function logic
  # ...

  # Return the result
  result <- list(param3 = param3, param4 = param4)
  return(result)
}

# Call the outer function with example values
param_list <- list(param3 = "value3", param4 = "value4")
result <- my_outer_function("value1", "value2", param_list)

# This needs to be tested but it has potential
# There might be some issues calling the functions like rsmp with do.call since these are syntactic sugar.
# Following works:

fsParams <- list(list("repeated_cv", repeats = 10, folds = 2), list("surv.cindex"))
obj <-do.call(mlr3::Resampling$new, fsParams[[1]])
obj <- do.call(mlr3::msr, fsParams[[2]])
obj
msr("surv.cindex")


obj <- Resampling$new("repeated_cv", repeats = 10, folds = 2)
