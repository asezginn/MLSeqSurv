
# All of the generics defined for surv function.

surv.blackboost <- function(data = data, method = "blackboost", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.blackboost", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.cforest <- function(data = data, method = "cforest", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.cforest", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)




}

surv.coxboost <- function(data = data, method = "coxboost", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.coxboost", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.coxtime <- function(data = data, method = "coxtime", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.coxtime", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.ctree <- function(data = data, method = "ctree", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.ctree", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.deephit <- function(data = data, method = "deephit", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.deephit", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.deepsurv <- function(data = data, method = "deepsurv", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.deepsurv", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.elasticnet <- function(data = data, method = "elasticnet", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.elasticnet", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.gamboost <- function(data = data, method = "gamboost", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.gamboost", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.gbm <- function(data = data, method = "gbm", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.gbm", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.glmboost <- function(data = data, method = "glmboost", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.glmboost", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.lasso <- function(data = data, method = "lasso", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.lasso", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.loghaz <- function(data = data, method = "loghaz", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.loghaz", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.mboost <- function(data = data, method = "mboost", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.mboost", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.obliqueRSF <- function(data = data, method = "obliqueRSF", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.obliqueRSF", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.pchazard <- function(data = data, method = "pchazard", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.pchazard", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.penalized <- function(data = data, method = "penalized", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.penalized", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.ranger <- function(data = data, method = "ranger", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.ranger", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.rfsrc <- function(data = data, method = "rfsrc", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.rfsrc", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.ridge <- function(data = data, method = "ridge", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.ridge", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.rpart <- function(data = data, method = "rpart", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.rpart", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.svm <- function(data = data, method = "svm", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.svm", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.xgboost_dart <- function(data = data, method = "xgboost_dart", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.xgboost_dart", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.xgboost_gblinear <- function(data = data, method = "xgboost_gblinear", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.xgboost_gblinear", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}

surv.xgboost_gbtree <- function(data = data, method = "xgboost_gbtree", preProcessing = preProcessing, paramGrid, ...){

  # this function might return an S4 object as data
  # S4 object could contain: train data set, test data set, time column, event column, original data
  data <- preprocess(data, preProcessing)

  learner <- lrn("surv.xgboost_gbtree", id = method, mstop = 100, ...) # will take parameters dynamically

  task_fs = TaskSurv$new("task_fs", data@train, time = "time", event = "status")

  measures = msrs(c("surv.cindex"))

  instance = FSelectInstanceSingleCrit$new( # will take parameters dynamically
    task = task_fs,
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("surv.cindex"),
    terminator = trm("evals", n_evals = 5)
  )

  fselector = fs("random_search") # may take parameters
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

  tune_ps = ParamSet$new(paramGrid)

  at=AutoTuner$new(learner=learner,
                   resampling=rsmp("repeated_cv", repeats = 10, folds = 5),
                   measure=msr("surv.cindex"),
                   terminator=trm("evals", n_evals = 5),
                   tuner=tnr("random_search"),
                   search_space = tune_ps
  )

  at$train(task_tune)
  pred_result <- at$predict_newdata(newdata = data@test)$score(measures)





}



