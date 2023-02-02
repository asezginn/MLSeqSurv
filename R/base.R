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
setMethod("classify", "MLSeqSurv", function(data, method = c("ipflasso", "prioritylasso"), preProcessing = c("deseq", "voom")) {

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


  ## Add relevant pre processing steps here.
  ## Might be cleaner to handle preprocessing in a different function.

  data <- preprocess(data, preProcessing)
  # This method checks the parameter "preProcessing" and does preProcessing according to the value specified.
  # If the value is null, it shouldn't do anything and simply return the original data.


  if (method == "prioritylasso"){
    model <- prioritylasso(X = data, Y = )
  }
  else if (method == "ipflasso"){
    model <- cvr.ipflasso(X = data, Y = )
  }

  return(model)

})
