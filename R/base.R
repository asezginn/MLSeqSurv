# Define the custom class "MLSeqSurv" with slots "train" and "test"
setClass("MLSeqSurv",
         slots = list(
           model = "ANY",
           cindex_train = "numeric",
           times = "ANY",
           cindex_test = "numeric",
           survHazards = "matrix",
           train = "data.frame",
           test = "data.frame",
           preprocessed_train = "data.frame",
           preprocessed_test = "data.frame"
         )
)




# Constructor for the MLSeqSurv class
# Add new fields as needed.
# This should be enough for now since we don't inherently need other fields such as event_column, time_column etc.
MLSeqSurv <- function(train, test) {
  obj <- new("MLSeqSurv")

  # Set the values of the "train" and "test" slots
  obj@train <- train
  obj@test <- test

  # Return the created object
  return(obj)
}

get_available_methods <- function(){

  all_methods_vector <- c("ipflasso", "prioritylasso", "blackboost", "cforest",
                          "coxboost", "coxtime", "ctree", "deephit", "deepsurv",
                          "elasticnet", "gamboost", "gbm", "glmboost", "lasso",
                          "loghaz", "mboost", "obliqueRSF", "pchazard", "penalized",
                          "ranger", "rfsrc", "ridge", "rpart", "svm",
                          "xgboost_dart", "xgboost_gblinear", "xgboost_gbtree")
  return(all_methods_vector)

}

survival <- function(data, method = c("ipflasso", "prioritylasso", "blackboost", "cforest",
                                  "coxboost", "coxtime", "ctree", "deephit", "deepsurv",
                                  "elasticnet", "gamboost", "gbm", "glmboost", "lasso",
                                  "loghaz", "mboost", "obliqueRSF", "pchazard", "penalized",
                                  "ranger", "rfsrc", "ridge", "rpart", "svm",
                                  "xgboost_dart", "xgboost_gblinear", "xgboost_gbtree"),
                 preProcessing = c("deseq-vst", "deseq-voom"),
                 fsParams, trainParams, balancer, tuneGrid, ...) {

  if (is.null(data)){
    stop("MLSeqSurv object is null")
  }

  if (is.null(method)){
    stop("Method is not specified.")
  }

  if (!(method %in% get_available_methods())){
    stop("Specified method is not available.")
  }

  # combine data into an S4 object
  # data <- MLSeqSurv(train_data, test_data)

  if (method == "blackboost"){
    survival.blackboost(data = data, method = "blackboost", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "cforest"){
    survival.cforest(data = data, method = "cforest", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "coxboost"){
    survival.coxboost(data = data, method = "coxboost", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "coxtime"){
    survival.coxtime(data = data, method = "coxtime", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "ctree"){
    survival.ctree(data = data, method = "ctree", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "deephit"){
    survival.deephit(data = data, method = "deephit", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "deepsurv"){
    survival.deepsurv(data = data, method = "deepsurv", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "gamboost"){
    survival.gamboost(data = data, method = "gamboost", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "gbm"){
    survival.gbm(data = data, method = "gbm", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "glmboost"){
    survival.glmboost(data = data, method = "glmboost", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "glmnet"){
    survival.glmnet(data = data, method = "glmnet", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "loghaz"){
    survival.loghaz(data = data, method = "loghaz", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "mboost"){
    survival.mboost(data = data, method = "mboost", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "obliqueRSF"){
    survival.obliqueRSF(data = data, method = "obliqueRSF", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "pchazard"){
    survival.pchazard(data = data, method = "pchazard", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "penalized"){
    survival.penalized(data = data, method = "penalized", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "ranger"){
    survival.ranger(data = data, method = "ranger", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "rfsrc"){
    survival.rfsrc(data = data, method = "rfsrc", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "rpart"){
    survival.rpart(data = data, method = "rpart", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "svm"){
    survival.svm(data = data, method = "svm", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "xgboost"){
    xgboost(data = data, method = "xgboost", fsParams, trainParams, tuneGrid, ...)
  }
  else if (method == "ipflasso"){
    survival.ipflasso(data = data, method = "ipflasso", balancer = balancer, trainParams = trainParams, ...)
  }
  else if (method == "prioritylasso"){
    survival.prioritylasso(data = data, method = "prioritylasso", balancer = balancer, trainParams = trainParams, ...)
  }

  ### Alternative approach is to call generic methods as soon as we take the input from the user via UseMethod.
  ### This allows this function to be relatively clean block of ifs and the actual work is done inside the methods which can call other generic functions


  ## Add relevant pre processing steps here.
  ## Might be cleaner to handle preprocessing in a different function.

  # data <- preprocess(data, preProcessing)
  # This method checks the parameter "preProcessing" and does preProcessing according to the value specified.
  # If the value is null, it shouldn't do anything and simply return the original data.


  # if (method == "prioritylasso"){
  #   model <- prioritylasso(X = data, Y = )
  # }
  # else if (method == "ipflasso"){
  #   model <- cvr.ipflasso(X = data, Y = )
  # }
  # else if (method == "blackboost"){
  #
  # }
  #
  # return(model)

}


survPredict <- function(model){

  if (is.null(model)){
    stop("Model is null")
  }

  if ("prioritylasso" %in% class(model@model)) {
    print("prioritylasso")
    survPredict.PL(model = model)
  } else if ("list" %in% class(model@model)) {
    print("IPFlasso")
    survPredict.IPF(model = model)
  } else if ("R6" %in% class(model@model)) {
    print("MLR3 model")
    survPredict.MLR(model = model)
  } else {
    stop("This model is not supported")
  }

}


paramListConstructor <- function(df){

  param_list <- list()

  for (i in 1:nrow(df)){
    row <- df[i,]
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
      stop(paste0("The parameter type ", row$type, " is not supported"))
    }

  }

  return(param_list)
}

exampleParamDf <- function(){

  row1 <- c("numeric", "nu", 0, 1)
  example_df <- data.frame(type = row1[1], paramid = row1[2], lower = row1[3], upper = row1[4])

  example_df <- rbind(example_df, row1)
  example_df <- rbind(example_df, row1)
  example_df <- rbind(example_df, row1)

  example_df$levels <- list(NA, c("coxph","weibull"), NA, NA)
  example_df$type <- list("numeric", "character", "integer", "logical")
  example_df$paramid <- list("mtry.ratio","splitrule","ntree","membership")
  example_df$lower <- list(0, NA, 500, NA)
  example_df$upper <- list(0.7, NA, 1500, NA)
  example_df$levels <- list(NA, c("logrank","bs.gradient"), NA, NA)

  return(example_df)
}
