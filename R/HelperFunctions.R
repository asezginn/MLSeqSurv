#### The functions of IPFLASSO ####

auc_n=function(y,prob,w){
  if(missing(w))
    concordance(y~prob)$concordance
  else concordance(y~prob,weights=w)$concordance
}

auc.mat_n=function(y,prob,weights=rep(1,nrow(y))){
  Weights=as.vector(weights*y)
  ny=nrow(y)
  Y=rep(c(0,1),c(ny,ny))
  Prob=c(prob,prob)
  auc_n(Y,Prob,Weights)
}

cvcompute <-
  function (cvstuff, foldid, nlams)
  {
    weights=cvstuff$weights
    mat=cvstuff$cvraw
    wisum = tapply(weights, foldid, sum)
    nfolds = max(foldid)
    outmat = matrix(NA, nfolds, ncol(mat))
    good = matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] = NA
    for (i in seq(nfolds)) {
      mati = mat[foldid == i, , drop = FALSE]
      wi = weights[foldid == i]
      outmat[i, ] = apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
      good[i, seq(nlams[i])] = 1
    }
    N = apply(good, 2, sum)
    list(cvraw = outmat, weights = wisum, N = N, type.measure=cvstuff$type.measure)
  }

cv.lognet <- function(predmat,y,type.measure,weights,foldid,grouped){
  prob_min = 1e-05
  prob_max = 1 - prob_min
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), , drop=FALSE]
  }
  N = nrow(y)
  nfolds = max(foldid)
  if ((N/nfolds < 10) && type.measure == "auc") {
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",
            call. = FALSE)
    type.measure = cvtype("deviance", "lognet")
  }
  predmat=1/(1+exp(-predmat))
  nlambda=ncol(predmat)
  nlams=rep(nlambda,nfolds)
  if (type.measure == "auc") {
    cvraw = matrix(NA, nfolds, nlambda)
    good = matrix(0, nfolds, nlambda)
    for (i in seq(nfolds)) {
      good[i, seq(nlams[i])] = 1
      which = foldid == i
      for (j in seq(nlams[i])) {
        cvraw[i, j] = auc.mat_n(y[which, ], predmat[which,
                                                    j], weights[which])
      }
    }
    N = apply(good, 2, sum)
    weights = tapply(weights, foldid, sum)
    grouped=FALSE
  }
  else {
    ywt = apply(y, 1, sum)
    y = y/ywt
    weights = weights * ywt
    N = nrow(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure,
                   mse = (y[, 1] - (1 - predmat))^2 + (y[, 2] - predmat)^2,
                   mae = abs(y[, 1] - (1 - predmat)) + abs(y[, 2] - predmat),
                   deviance = {
                     predmat = pmin(pmax(predmat, prob_min), prob_max)
                     lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                     ly = log(y)
                     ly[y == 0] = 0
                     ly = drop((y * ly) %*% c(1, 1))
                     2 * (ly - lp)
                   },
                   class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5)
    )
  }
  list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=grouped)
}

cvr.ipflasso_n <- function (X, Y, family, type.measure, standardize = TRUE, alpha = 1,
                            blocks, pf, nfolds, ncv, weights)
{
  M <- length(blocks)
  if (M != length(pf)) {
    stop("blocks and pf must have the same length.")
  }
  ulblocks <- as.numeric(unlist(blocks))
  if (!setequal(ulblocks, 1:ncol(X))) {
    stop("Each predictor should be included in exactly one block.")
  }
  if (family == "gaussian") {
    if (type.measure != "mse")
      warning("type.measure is set to mse.")
    type.measure <- "mse"
  }
  if (family == "cox") {
    if (type.measure != "deviance")
      warning("type.measure is set to partial likelihood.")
    type.measure <- "deviance"
  }
  if (family == "binomial" & !is.element(type.measure, c("auc",
                                                         "class"))) {
    warning("type.measure is set to class")
    type.measure <- "class"
  }
  if (any(pf < 0)) {
    stop("pf should have positive entries.")
  }
  pfvector <- numeric(ncol(X))
  for (m in 1:M) {
    pfvector[blocks[[m]]] <- pf[m]
  }
  a <- cvr.glmnet_n(X = X, Y = Y, weights = weights, family = family, standardize = standardize,
                    alpha = alpha, nfolds = nfolds, ncv = ncv, type.measure = type.measure,
                    penalty.factor = pfvector)
  coeff <- a$coeff
  if (family != "cox") {
    rownames(coeff)[1] <- "intercept"
  }
  if (type.measure == "auc") {
    ind.bestlambda <- which.max(a$cvm)
  }
  if (type.measure != "auc") {
    ind.bestlambda <- which.min(a$cvm)
  }
  nzero <- apply(coeff[-1, ], FUN = function(x) return(sum(x !=
                                                             0)), MARGIN = 2)
  return(list(coeff = coeff, ind.bestlambda = ind.bestlambda,
              lambda = a$lambda, cvm = a$cvm, nzero = nzero, family = family))
}




cvr.glmnet_n <- function (X, Y, weights, family, standardize = TRUE, alpha = 1, nfolds,
                          ncv, type.measure,  ...)
{
  set.seed(1)
  a1 <- cv.glmnet_n(x = X, y = Y, weights = weights, family = family, standardize = standardize,
                    nfolds = nfolds, type.measure = type.measure, alpha = alpha,
                    ...)
  lambda1 <- a1$lambda
  cvm <- a1$cvm
  if (ncv == 1) {
    return(list(coeff = rbind(as.numeric(a1$glmnet.fit$a0)[1:length(lambda1)],
                              as.matrix(a1$glmnet.fit$beta)[, 1:length(lambda1)]),
                lambda = lambda1, cvm = cvm))
  }
  else {
    lambda <- lambda1
    cvm <- matrix(cvm, nrow = 1)
    for (i in 2:ncv) {
      set.seed(i)
      a <- cv.glmnet_n(x = X, y = Y, weights = weights, family = family, standardize = standardize,
                       nfolds = nfolds, type.measure = type.measure,
                       lambda = lambda1, alpha = alpha, ...)
      newlambda <- intersect(lambda, a$lambda)
      cvm <- rbind(cvm[, is.element(lambda, newlambda)],
                   a$cvm[is.element(a$lambda, newlambda)])
      lambda <- newlambda
    }
    coeff <- rbind(as.numeric(a1$glmnet.fit$a0)[1:length(lambda)],
                   as.matrix(a1$glmnet.fit$beta)[, 1:length(lambda)])
    return(list(coeff = coeff, lambda = lambda, cvm = apply(cvm,
                                                            MARGIN = 2, FUN = mean)))
  }
}



cv.glmnet_n <- function (x, y, weights = NULL, offset = NULL, lambda = NULL,
                         type.measure = c("default", "mse", "deviance", "class", "auc",
                                          "mae", "C"), nfolds = 10, foldid = NULL, alignment = c("lambda",
                                                                                                 "fraction"), grouped = TRUE, keep = FALSE, parallel = FALSE,
                         gamma = c(0, 0.25, 0.5, 0.75, 1), relax = FALSE, trace.it = 0,
                         ...)
{
  type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv.glmnet")
  if (!is.null(lambda) && alignment == "fraction") {
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment = "lambda"
  }
  N = nrow(x)
  if (is.null(weights))
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  cv.call = glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped",
                  "keep"), names(glmnet.call), FALSE)
  if (any(which))
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  if (glmnet.control()$itrace)
    trace.it = 1
  else {
    if (trace.it) {
      glmnet.control(itrace = 1)
      on.exit(glmnet.control(itrace = 0))
    }
  }
  if (is.null(foldid))
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  if (relax)
    cv.relaxed.raw_n(x, y, weights, offset, lambda, type.measure,
                     nfolds, foldid, alignment, grouped, keep, parallel,
                     trace.it, glmnet.call, cv.call, gamma, ...)
  else cv.glmnet.raw_n(x, y, weights, offset, lambda, type.measure,
                       nfolds, foldid, alignment, grouped, keep, parallel, trace.it,
                       glmnet.call, cv.call, ...)
}

cvtype_n <- function(type.measure="mse",subclass="elnet"){
  type.measures = c("mse","deviance", "class", "auc", "mae","C")
  devname=switch(subclass,
                 elnet="Mean-squared Error",
                 lognet="Binomial Deviance",
                 fishnet="Poisson Deviance",
                 coxnet="Partial Likelihood Deviance",
                 multnet="Multinomial Deviance",
                 mrelnet="Mean-squared Error",
                 glmnetfit="GLM Deviance"
  )
  typenames = c(deviance = devname, mse = "Mean-Squared Error",
                mae = "Mean Absolute Error",auc = "AUC", class = "Misclassification Error",C="C-index")
  subclass.ch=switch(subclass,
                     elnet=c(1,2,5),
                     lognet=c(2,3,4,1,5),
                     fishnet=c(2,1,5),
                     coxnet=c(2,6),
                     multnet=c(2,3,1,5),
                     mrelnet=c(1,2,5),
                     glmnetfit=c(2,1,5)
  )
  subclass.type=type.measures[subclass.ch]
  if(type.measure=="default")type.measure=subclass.type[1]
  model.name=switch(subclass,
                    elnet="Gaussian",
                    lognet="Binomial",
                    fishnet="Poisson",
                    coxnet="Cox",
                    multnet="Multinomial",
                    mrelnet="Multi-response Gaussian",
                    glmnetfit="GLM"
  )
  if(!match(type.measure,subclass.type,FALSE)){
    type.measure=subclass.type[1]
    warning(paste("Only ",paste(subclass.type,collapse=", ")," available as type.measure for ",model.name," models; ", type.measure," used instead",sep=""),call.=FALSE)
  }
  names(type.measure)=typenames[type.measure]
  type.measure
}


cvstats_n=function(cvstuff,foldid,nfolds,lambda,nz,grouped,...){
  if (grouped){
    nlams=rep(dim(cvstuff$cvraw)[2],nfolds) ## All the same - this should go
    cvstuff= cvcompute(cvstuff, foldid, nlams)
  }
  cvm=with(cvstuff,apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE))
  cvsd=with(cvstuff, sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                                w = weights, na.rm = TRUE)/(N - 1)))
  nas=is.na(cvsd)
  if(any(nas)){
    lambda=lambda[!nas]
    cvm=cvm[!nas]
    cvsd=cvsd[!nas]
    nz=nz[!nas]
  }
  list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
         cvsd, cvlo = cvm - cvsd, nzero = nz)
}


getOptcv.glmnet_n <-
  function (lambda, cvm, cvsd,cvname)
  {
    if(match(cvname,c("AUC","C-index"),0))cvm=-cvm
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + cvsd)[idmin]
    id1se = cvm <= semin
    lambda.1se = max(lambda[id1se], na.rm = TRUE)
    id1se = match(lambda.1se, lambda)
    index=matrix(c(idmin,id1se),2,1,dimnames=list(c("min","1se"),"Lambda"))
    list(lambda.min = lambda.min, lambda.1se = lambda.1se, index = index)
  }

cv.glmnet.raw_n <-
  function (x, y, weights, offset, lambda, type.measure, nfolds, foldid, alignment,grouped, keep,
            parallel, trace.it, glmnet.call, cv.call, ...)
  {


    if (trace.it) cat("Training\n")
    glmnet.object = glmnet(x, y, weights = weights, offset = offset,
                           lambda = lambda, trace.it=trace.it,...)
    glmnet.object$call = glmnet.call
    subclass=class(glmnet.object)[[1]]
    type.measure=cvtype_n(type.measure,subclass)
    is.offset = glmnet.object$offset
    ###Next line is commented out so each call generates its own lambda sequence
    # lambda=glmnet.object$lambda
    if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
      nz = predict(glmnet.object, type = "nonzero")
      nz = sapply(nz, function(x) sapply(x, length))
      nz = ceiling(apply(nz, 1, median))
    }
    else nz = sapply(predict(glmnet.object, type = "nonzero"),
                     length)
    outlist = as.list(seq(nfolds))
    N=nrow(x)
    if (parallel) {
      #  if (parallel && require(foreach)) {
      outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
        {
          which = foldid == i
          #      if (is.matrix(y))
          if (length(dim(y))>1)
            y_sub = y[!which, ]
          else y_sub = y[!which]
          if (is.offset)
            offset_sub = as.matrix(offset)[!which, ]
          else offset_sub = NULL
          glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
                 offset = offset_sub, weights = weights[!which],
                 ...)
        }
    }
    else {
      for (i in seq(nfolds)) {
        if (trace.it) cat(sprintf("Fold: %d/%d\n", i, nfolds))
        which = foldid == i
        if (length(dim(y))>1)
          y_sub = y[!which, ]
        else y_sub = y[!which]
        if (is.offset)
          offset_sub = as.matrix(offset)[!which, ]
        else offset_sub = NULL
        outlist[[i]] = glmnet(x[!which, , drop = FALSE],
                              y_sub, lambda = lambda, offset = offset_sub,
                              weights = weights[!which],trace.it=trace.it, ...)
      }
    }
    lambda = glmnet.object$lambda
    class(outlist)=paste0(subclass,"list")
    predmat=buildPredmat(outlist,lambda,x,offset,foldid,alignment,y=y,weights=weights,
                         grouped=grouped,type.measure=type.measure,family=family(glmnet.object))
    ### we include type.measure for the special case of coxnet with the deviance vs C-index discrepancy
    ### family is included for the new GLM crowd
    ### Next we compute the measures
    #    if(subclass=="glmnetfit") attr(predmat,"family")=glmnet.object$family
    fun = paste("cv", subclass, sep = ".")
    cvstuff = do.call(fun, list(predmat,y,type.measure,weights,foldid,grouped))

    grouped=cvstuff$grouped
    if ((N/nfolds < 3) && grouped) {
      warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
              call. = FALSE)
      grouped = FALSE
    }

    out=cvstats_n(cvstuff,foldid,nfolds,lambda,nz,grouped)
    cvname = names(cvstuff$type.measure)
    names(cvname)=cvstuff$type.measure# to be compatible with earlier version; silly, I know
    out = c(out,list(call=cv.call,name = cvname, glmnet.fit = glmnet.object))
    if (keep)
      out = c(out, list(fit.preval = predmat, foldid = foldid))
    lamin=with(out,getOptcv.glmnet_n(lambda, cvm, cvsd, cvname))
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    obj
  }

#### The functions of the stacking ####
get_event_count <- function(data_df){
  event_count <- 0
  time_check <- 0

  data_df <- data_df[order(data_df[,c("time")],decreasing=FALSE),]

  for (i in 1:nrow(data_df)){ # calculates the total amount of columns

    if (data_df$status[i] == 1 && time_check != data_df$time[i]){
      time_check <- data_df$time[i]
      event_count <- event_count + 1
    }

  }
  return(event_count)
}
calculate_rows <- function(data_df){ # this function calculates the total amount of rows.

  row_count <- 0
  time_check <- 0

  for (i in 1:nrow(data_df)){ # Calculating total risk set length



    if (data_df$status[i] == 1 && data_df$time[i] != time_check){

      time_check <- data_df$time[i]

      for (j in 1:nrow(data_df)){ # a loop to find the length of the risk set at given time interval. We do this because we want to get the first instance for accurate results

        if (data_df$time[j] == data_df$time[i])  {

          risk_set_length <- nrow(data_df) - j + 1 # amount of people in the risk set at a given time interval
          row_count <- row_count + risk_set_length
          break
        }

      }

    }
  }
  return(row_count)

}
calculate_columns <- function(data_df, return_risk_ncol = FALSE){ # This function calculates the total amount of columns

  data_cov_df <- subset(data_df, select = -c(time, status)) # subsets the original dataset into the covariant set in order to calculate total number of columns
  cov_count <- ncol(data_cov_df)
  column_count <- ncol(data_cov_df)
  risk_matrix_ncol <- 0
  time_check <- 0

  for (i in 1:nrow(data_df)){ # calculates the total amount of columns

    if (data_df$status[i] == 1 && time_check != data_df$time[i]){
      time_check <- data_df$time[i]
      column_count <- column_count+1
      risk_matrix_ncol <- risk_matrix_ncol + 1
    }

  }
  if (return_risk_ncol){
    return(risk_matrix_ncol)
  }
  else{
    return(column_count)
  }
}
# Fill left matrix
create_left_matrix <- function(data_df, data_cov_df){

  column_count <- calculate_columns(data_df)
  row_count <- calculate_rows(data_df)
  cov_count <- ncol(data_cov_df)

  left_matrix <- matrix(0, nrow = row_count, ncol = column_count)
  left_df <- as.data.frame(left_matrix)
  time_column <- numeric(row_count)

  # These indexes help us operate on the data inside the loop
  row_index <- 1
  col_index <- 1
  time_check <- 0

  for (i in 1:nrow(data_df)){


    if (data_df$status[i] == 1 && data_df$time[i] != time_check){

      time_check <- data_df$time[i]

      for (j in 1:nrow(data_df)){

        if (data_df$time[j] == data_df$time[i]){

          risk_set_length <- nrow(data_df) - j + 1
          left_df[row_index:(row_index+risk_set_length-1),col_index] <- 1
          time_column[row_index:(row_index+risk_set_length-1)] <- data_df$time[j]
          left_df[row_index:(row_index+risk_set_length-1),(1+column_count-cov_count):column_count] <- data_cov_df[j:nrow(data_cov_df),] # covariant matrix
          col_index <- col_index + 1
          row_index <- row_index + risk_set_length

          break
        }

      }
    }
  }
  # left_df <- left_df[,-(which(colSums(left_df)==0))]
  left_df <- cbind(left_df, time_column)
  return(left_df)

}
create_prediction_column <- function(data_df){

  prediction_column <- numeric(calculate_rows(data_df))
  row_index <- 1
  repeat_count <- 0
  risk_set_length <- 0
  time_check <- 0

  for (i in 1:nrow(data_df)) {

    if (data_df$status[i] == 1 && data_df$time[i] != time_check){ # status is 1 here

      time_check <- data_df$time[i] # this variable exists so that outer loop does not execute on duplicate rows

      for (j in 1:nrow(data_df)){ # starting another loop to find the first instance where that status equals 1

        if (data_df$time[j] == data_df$time[i] && data_df$status[j] == 1) { # this if iterates through all of the entries which have the same time and status = 1

          repeat_count <- repeat_count + 1
        }


      }
      for (j in 1:nrow(data_df)){ # a loop to find the length of the risk set

        if (data_df$time[j] == data_df$time[i])  {

          risk_set_length <- nrow(data_df) - j + 1

          break
        }

      }

      prediction_column[row_index:(row_index + repeat_count - 1)] <- 1
      row_index <- row_index + risk_set_length
      repeat_count <- 0

    }

  }
  return(prediction_column)

}
stack_df <- function(data_df){

  # To Do
  # make sure data_df is a dataframe. if not convert it to one

  # check for time and status columns, if these columns do not exist warn the user.

  # sort data_df by time

  # create the left matrix

  # create the prediction column

  # bind them and output the result as df

  # keep the names of the covariates if possible

  if (!is.data.frame(data_df)){
    data_df <- as.data.frame(data_df)
  }

  if (is.null(data_df$time) || is.null(data_df$status)){
    stop("Please enter a valid dataframe with time and status columns!")
  }
  # extra checks for time and status columns


  data_df <- data_df[order(data_df[,c("time")],decreasing=FALSE),] # sorting the dataframe by time
  data_cov_df <- subset(data_df, select = -c(time, status))
  left_matrix <- create_left_matrix(data_df = data_df, data_cov_df = data_cov_df)
  column_count <- calculate_columns(data_df)
  # colnames(left_matrix[,(column_count-ncol(data_cov_df)+1):column_count]) <- colnames(data_cov_df)
  prediction_column <- create_prediction_column(data_df)
  final_df <- cbind.data.frame(left_matrix, prediction_column)

  return(final_df)

}
stacked_matrix_summarize <- function(stacked_df){

  # prediction_column <- stacked_df$prediction_column
  # stacked_df <- subset(stacked_df, select = -c(prediction_column))

  summarised_df <- stacked_df %>%
    group_by(time_column) %>%
    summarise_all(mean)
  #summarize(across(everything(), mean))
  return(summarised_df)
}

#### The functions of the test datasets manipulation ####

adjust_test_row <- function(summarized_df, test_df_row, test_df_time, eventCount){ # added commas to row indexes because R thinks they are dataframes for some reason



  if (!is.null(summarized_df$prediction_column)){
    summarized_df <- subset(summarized_df, select = -c(prediction_column))
  }

  time_index <- sum(summarized_df$time_column <= test_df_time) # this could return 0. Need to create a special case for it, also if it is larger than the max value need to add a new row


  if (time_index == 0){ # smaller than min, add time set? if we add new time set do we also need to add new risk set?
    # tempRow <- c(numeric(eventCount), test_df_row)
    # names(tempRow) <- colnames(summarized_df)
    # summarized_df <- rbind(tempRow, summarized_df) # this would result in a row thats completely filled with 0s. We can remove this row later on
    # time_index = 1
    #return(numeric((length(test_df_row)+eventCount-1)))
    result_row <- as.numeric(test_df_row[,-1])
    result_row <- c(numeric(eventCount), result_row)
  }

  else if (time_index == nrow(summarized_df)){
    if (test_df_row[,1] > summarized_df$time_column[time_index]){ # bigger than max, add time set? if we add new time set do we also need to add new risk set?
      result_row <- as.numeric(test_df_row[,-1]) - as.numeric(summarized_df[nrow(summarized_df),(eventCount+2):ncol(summarized_df)])
      result_row <- c(numeric(eventCount), result_row)
      # print("here")
    }
  }
  else{
    result_row <- as.numeric(test_df_row[,-1]) - as.numeric(summarized_df[time_index,(eventCount+2):ncol(summarized_df)]) # this pulls EVERYTHING besides the risk matrix and time column. For this to work make sure that the summarized_df does not contain the prediction_column

    result_row <- c(numeric(eventCount), result_row)
  }
  return(result_row)

}

integrate_test_data <- function(summarized_df, test_df, eventCount){

  # Need to add appopriate checks for inputs

  test_df_result <- data_frame()

  for (i in 1:nrow(test_df)){

    test_row_result <- adjust_test_row(summarized_df, test_df[i,], test_df[i,1], eventCount)
    # compare_colnames <<- c(compare_colnames, colnames(test_df_result))
    colnames(test_df_result) <- names(test_row_result)
    test_df_result <- rbind(test_df_result, as.vector(test_row_result))

  }
  test_df_result <- test_df_result[apply(test_df_result[,-1], 1, function(x) !all(x==0)),] # removes rows that have only 0's in them
  return(test_df_result)
}

adjust_test_data <- function(test_data){

  test_time <- test_data$time
  processedTest <- test_data
  processedTest <- subset(processedTest, select = -c(time,status))
  final_test <- cbind(test_time, processedTest)
  names(final_test)[names(final_test) == 'test_time'] <- 'time'
  return(final_test)

}

adjust_stacked_df <- function(stacked_data){

  stacked_data_time <- stacked_data$time_column
  stacked_data <- stacked_data
  stacked_data <- subset(stacked_data, select = -c(time_column))
  stacked_data <- cbind(stacked_data_time, stacked_data)
  names(stacked_data)[names(stacked_data) == 'stacked_data_time'] <- 'time_column'
  return(stacked_data)
}


#### Helper Functions ####
get_matching_probabilities <- function(summarized_df, test_df){ # can get na's. should probably replace them with zero

  time_pred_df <- subset(summarized_df, select=c("time_column", "prediction_column"))
  # browser()
  print(time_pred_df$time_column)
  print(as.vector(test_df[,1]))
  probs <- list()

  for (i in 1:nrow(test_df)){

    test_df_time <- test_df[i,1]
    test_df_time <- as.numeric(test_df_time)
    # print(test_df_time)

    time_index <- sum(time_pred_df$time_column <= test_df_time) # this should determine the risk set each row belongs to. Currently it gives result-1
    # probs <- c(probs,time_pred_df[(time_index+1),1])
    # print(as.numeric(time_pred_df[(time_index+1),1]))
    numToAppend <- as.numeric(time_pred_df[(time_index+1),2])
    probs <- append(probs, numToAppend)
    # print(time_index)


  }
  return(probs)

}

getFeaturesFromPrioLasso <- function(inputModel){ # doesn't return index as of now

  coeffdf <- inputModel$coefficients
  coeffdf <- as.data.frame(coeffdf)
  selectedFeatures <- which(coeffdf != 0)
  return(selectedFeatures)
}

getFeaturesFromIPF <- function(inputModel){

  bestLambda <- inputModel$ind.bestlambda
  df <- inputModel$coeff
  bestFeatures <- df[,bestLambda]
  bestFeatures <- as.data.frame(bestFeatures)
  bestFeatures <- which(bestFeatures != 0)
  return(bestFeatures)

}
