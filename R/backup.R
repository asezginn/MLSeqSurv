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

  all_methods_vector <- c()
  return(all_methods_vector)

}

#' @rdname MLSeqSurv-class
#' @aliases MLSeqSurv-classW
setMethod("classify", "MLSeqSurv", function(data, method = c("ipflasso", "prioritylasso"), preProcessing = c("deseq", "voom"), seeded = TRUE) {


  if (is.null(method)){
    stop("Method is not specified.")
  }


  ## Add relevant pre processing steps here.
  ## Might be cleaner to handle preprocessing in a different function.

  data <- preProcess(data, preProcessing)
  # This method checks the parameter "preProcessing" and does preProcessing according to the value specified.
  # If the value is null, it shouldn't do anything and simply return the value.

  task <- TaskSurv$new(data, time = "time", event = "status")
  # creating a task object from mlr3.
  # might be more reasonable to switch it to "event" instead of "status"

  learner <- lrn(method)
  # creating a learner with the specified method from mlr3.

  if (seeded){
    set.seed(1881)
  }

  splits <- partition(task) # ratio can be specified by user at a later stage as an additional parameter

  learner$train(task, splits$train)
  # training the learner on the task.




  return(data)
})
