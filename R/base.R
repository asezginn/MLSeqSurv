#' MLSeqSurv
#'
#' This package contains functions for training and predicting using machine learning models.
#'
#' @name MLSeqSurv

#' @docType package
#' @useDynLib MLSeqSurv, .registration = TRUE

setClass("MLSeqSurv", representation(model = "any"))

#' @rdname MLSeqSurv-class
#' @aliases MLSeqSurv-class
setMethod("predict", "MLSeqSurv", function(object, newdata) {
  # function content

  return(predict(object@model, newdata))
})

get_available_methods <- function(){

  all_methods_vector <- c("ipflasso", "prioritylasso")
  return(all_methods_vector)

}

#' @rdname MLSeqSurv-class
#' @aliases MLSeqSurv-classW
setMethod("surv", "MLSeqSurv", function(data, method = c("ipflasso", "prioritylasso"), preProcessing = c("deseq-vst", "deseq-voom"), paramGrid, ...) {

  # should data be an S4 object?
  if (is.null(data)){
    stop("Data is null")
  }
  if (!(is.data.frame(data))){
    data <- as.data.frame(data)
  }

  if (is.null(method)){
    stop("Method is not specified.")
  }

  if (!(method %in% get_available_methods())){
    stop("Specified method is not available.")
  }


  if (method == "blackboost"){
    surv.blackboost(data = data, method = "blackboost", preProcessing = preProcessing, paramGrid, ...) # Fix the parameters that are being sent once the generic functions are finalized.
  }
  else if (method == "cforest"){
    surv.cforest(data = data, method = "cforest", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "coxboost"){
    surv.coxboost(data = data, method = "coxboost", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "coxtime"){
    surv.coxtime(data = data, method = "coxtime", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "ctree"){
    surv.ctree(data = data, method = "ctree", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "deephit"){
    surv.deephit(data = data, method = "deephit", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "deepsurv"){
    surv.deepsurv(data = data, method = "deepsurv", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "elasticnet"){
    surv.elasticnet(data = data, method = "elasticnet", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "gamboost"){
    surv.gamboost(data = data, method = "gamboost", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "gbm"){
    surv.gbm(data = data, method = "gbm", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "glmboost"){
    surv.glmboost(data = data, method = "glmboost", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "lasso"){
    surv.lasso(data = data, method = "lasso", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "loghaz"){
    surv.loghaz(data = data, method = "loghaz", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "mboost"){
    surv.mboost(data = data, method = "mboost", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "obliqueRSF"){
    surv.obliqueRSF(data = data, method = "obliqueRSF", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "pchazard"){
    surv.pchazard(data = data, method = "pchazard", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "penalized"){
    surv.penalized(data = data, method = "penalized", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "ranger"){
    surv.ranger(data = data, method = "ranger", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "rfsrc"){
    surv.rfsrc(data = data, method = "rfsrc", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "ridge"){
    surv.ridge(data = data, method = "ridge", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "rpart"){
    surv.rpart(data = data, method = "rpart", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "svm"){
    surv.svm(data = data, method = "svm", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "xgboost_dart"){
    surv.xgboost_dart(data = data, method = "xgboost_dart", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "xgboost_gblinear"){
    surv.xgboost_gblinear(data = data, method = "xgboost_gblinear", preProcessing = preProcessing, paramGrid, ...)
  }
  else if (method == "xgboost_gbtree"){
    surv.xgboost_gbtree(data = data, method = "xgboost_gbtree", preProcessing = preProcessing, paramGrid, ...)
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

})
