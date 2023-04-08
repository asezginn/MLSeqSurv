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

# get_matching_probabilities(summ_stack, test_final)
# View(test_final)

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


#### Start of the main loop ####

cancerType <- "SARC"

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


# priolasso_selected_features <- getFeaturesFromPrioLasso(pl_bin1)
# ipf_selected_features <- getFeaturesFromIPF(model1)
for (index in 1:30) {
  
  results_df <- data.frame(matrix(ncol = 26, nrow = 0))
  results_time  <- data.frame(matrix(ncol = 26, nrow = 0))
  
  #### Reading train and test data ####
  print(paste0("Reading Data for ", index))
  tvsd_train3 <- fread(file = paste("MLSeqSurv_voomFeatures/",cancerType, "/train/train_", index, ".csv", sep = ""))
  tvsd_test3 <- fread(file = paste("MLSeqSurv_voomFeatures/",cancerType, "/test/test_", index, ".csv", sep = ""))
  
  #### Stacking ####
  print(paste0("Initializing stacking for " , index))
  stack_matris <- stack_df(tvsd_train3[,-1])
  
  summ_stack <- stacked_matrix_summarize(stack_matris)
  event_count <- get_event_count(tvsd_train3[,-1])
  
  stack_matris$prediction_column <- as.factor(stack_matris$prediction_column)
  
  ### Imbalancing problem ###
  over1 <- ROS(stack_matris, "prediction_column") 
  agirlik <- over1[,(ncol(over1)-2)]
  over <- over1[,-(ncol(over1)-2)]
  
  table(over$prediction_column)
  
  ##### Manipulation test data ####
  print("Adjusting test data!")
  summ_stack <- summ_stack[,-(ncol(summ_stack)-1)]
  over_final <- adjust_stacked_df(summ_stack)
  test_final <- adjust_test_data(tvsd_test3[,-1])
  transformed_test <- integrate_test_data(summarized_df = over_final, test_df = test_final, eventCount = event_count)
  
  probVector <- get_matching_probabilities(summ_stack, test_final)
  probVector <- probVector %>% replace(is.na(.), 0)
  probVector <- unlist(probVector)
  
  #### Priority-lasso ####
  print(paste0("Starting Priority Lasso for " , index))
  y3 <-as.numeric(over$prediction_column)
  y3[y3==1] <- 0
  y3[y3==2] <- 1
  # fixInNamespace("predict.glmnet", ns = "glmnet")
  # make newx numeric
  over_matrix <- as.matrix(over[,1:(ncol(over)-2)]) # results in a character array
  over_matrix <- apply(over_matrix, 2, as.numeric) # converts it into a numeric array
  # over_matrix <- over_matrix[,-40] # remove random features to test the result
  # over_matrix <- over_matrix[,-80]
  #### 1. model ####
  # pl_bin1 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
  #                          family = "binomial", type.measure = "class",
  #                          #blocks = list(block1=1:ncol(over_matrix)),
  #                          blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                          block1.penalization = TRUE, 
  #                          lambda.type = "lambda.min",
  #                          standardize = TRUE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  # newdata_bin1 <- as.matrix(transformed_test)
  # ypredscore1 <-predict(object = pl_bin1, newdata = newdata_bin1, type = "response")
  # # # pred_mat1 <- matrix(ypredscore1, nrow = nrow(ypredscore1), ncol = ncol(ypredscore1))
  # # pred_mat1 <- as.matrix(over_final$prediction_column)
  # 
  # 
  # pred1 <- prediction(ypredscore1,tvsd_test3$status)
  # pred1@predictions <- list(unlist(pred1@predictions) + probVector)
  # perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
  # plot(perf1,col='red')
  # abline(a=0,b=1)
  # auc1 <- performance(pred1,measure="auc")@y.values
  # auc1
  
  #print("Model 1 is done")
  #### 2. model ####
  # pl_bin2 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
  #                          family = "binomial", type.measure = "class",
  #                          blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                          block1.penalization = TRUE, lambda.type = "lambda.min",
  #                          standardize = FALSE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  # newdata_bin2 <- as.matrix(transformed_test)
  # ypredscore2 <-predict(object = pl_bin2, newdata = newdata_bin2, type = "response")
  # 
  # #ypredscore <-ypredscore+(1/2)
  # 
  # pred2 <- prediction(ypredscore2,tvsd_test3$status)
  # 
  # pred2@predictions <- list(unlist(pred2@predictions) + probVector)
  # perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
  # plot(perf2,col='red')
  # abline(a=0,b=1)
  # auc2 <- performance(pred2,measure="auc")@y.values
  # auc2
  
  #print("Model 2 is done")
  #### 3. model ####
  # pl_bin3 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
  #                          family = "binomial", type.measure = "class",
  #                          blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                          block1.penalization = TRUE, lambda.type = "lambda.1se",
  #                          standardize = TRUE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  # newdata_bin3 <- as.matrix(transformed_test)
  # ypredscore3 <-predict(object = pl_bin3, newdata = newdata_bin3, type = "response")
  # 
  # #ypredscore <-ypredscore+(1/2)
  # 
  # 
  # pred3 <- prediction(ypredscore3,tvsd_test3$status)
  # pred3@predictions <- list(unlist(pred3@predictions) + probVector)
  # perf3 <- performance(pred3, measure = "tpr", x.measure = "fpr")
  # plot(perf3,col='red')
  # abline(a=0,b=1)
  # auc3 <- performance(pred3,measure="auc")@y.values
  # auc3
  
  #print("Model 3 is done")
  #### 4. model ####
  # pl_bin4 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
  #                          family = "binomial", type.measure = "class",
  #                          blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                          block1.penalization = TRUE, lambda.type = "lambda.1se",
  #                          standardize = FALSE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  # newdata_bin4 <- as.matrix(transformed_test)
  # ypredscore4 <-predict(object = pl_bin4, newdata = newdata_bin4, type = "response")
  # 
  # #ypredscore <-ypredscore+(1/2)
  # 
  # 
  # pred4 <- prediction(ypredscore4,tvsd_test3$status)
  # pred4@predictions <- list(unlist(pred4@predictions) + probVector)
  # perf4 <- performance(pred4, measure = "tpr", x.measure = "fpr")
  # plot(perf4,col='red')
  # abline(a=0,b=1)
  # auc4 <- performance(pred4,measure="auc")@y.values
  # auc4
  
  #print("Model 4 is done")
  #### 5. model ####
  timer_start_tune <- Sys.time()
  pl_bin5 <- prioritylasso(X = over_matrix, Y = y3,weights = agirlik, 
                           family = "binomial", type.measure = "auc",
                           blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                           block1.penalization = TRUE, lambda.type = "lambda.min",
                           standardize = TRUE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin5 <- as.matrix(transformed_test)
  ypredscore5 <-predict(object = pl_bin5, newdata = newdata_bin5, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  
  pred5 <- prediction(ypredscore5,tvsd_test3$status)
  pred5@predictions <- list(unlist(pred5@predictions) + probVector)
  perf5 <- performance(pred5, measure = "tpr", x.measure = "fpr")
  plot(perf5,col='red')
  abline(a=0,b=1)
  auc5 <- performance(pred5,measure="auc")@y.values
  auc5
  timer_stop_tune <- Sys.time()
  auc5_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 5 is done")
  #### 6. model ####
  timer_start_tune <- Sys.time()
  pl_bin6 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
                           family = "binomial", type.measure = "auc",
                           # blocks = list(block1=1:ncol(over_matrix)),
                           blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                           block1.penalization = TRUE, lambda.type = "lambda.min",
                           standardize = FALSE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin6 <- as.matrix(transformed_test)
  ypredscore6 <-predict(object = pl_bin6, newdata = newdata_bin6, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  
  pred6 <- prediction(ypredscore6,tvsd_test3$status)
  pred6@predictions <- list(unlist(pred6@predictions) + probVector)
  perf6 <- performance(pred6, measure = "tpr", x.measure = "fpr")
  plot(perf6,col='red')
  abline(a=0,b=1)
  auc6 <- performance(pred6,measure="auc")@y.values
  auc6
  timer_stop_tune <- Sys.time()
  auc6_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 6 is done")
  
  #### 7. model ####
  timer_start_tune <- Sys.time()
  pl_bin7 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
                           family = "binomial", type.measure = "auc",
                           blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                           block1.penalization = TRUE, lambda.type = "lambda.1se",
                           standardize = TRUE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin7 <- as.matrix(transformed_test)
  ypredscore7 <-predict(object = pl_bin7, newdata = newdata_bin7, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  pred7 <- prediction(ypredscore7,tvsd_test3$status)
  pred7@predictions <- list(unlist(pred7@predictions) + probVector)
  perf7 <- performance(pred7, measure = "tpr", x.measure = "fpr")
  plot(perf7,col='red')
  abline(a=0,b=1)
  auc7 <- performance(pred7,measure="auc")@y.values
  auc7
  timer_stop_tune <- Sys.time()
  auc7_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 7 is done")
  
  #### 8. model ####
  timer_start_tune <- Sys.time()
  pl_bin8 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
                           family = "binomial", type.measure = "auc",
                           blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                           block1.penalization = TRUE, lambda.type = "lambda.1se",
                           standardize = FALSE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin8 <- as.matrix(transformed_test)
  ypredscore8 <-predict(object = pl_bin8, newdata = newdata_bin8, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  
  pred8 <- prediction(ypredscore8,tvsd_test3$status)
  pred8@predictions <- list(unlist(pred8@predictions) + probVector) # add mean of risk sets to the results
  perf8 <- performance(pred8, measure = "tpr", x.measure = "fpr")
  plot(perf8,col='red')
  abline(a=0,b=1)
  auc8 <- performance(pred8,measure="auc")@y.values
  auc8
  timer_stop_tune <- Sys.time()
  auc8_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 8 is done")
  
  
  #### 5.2 model ####
  timer_start_tune <- Sys.time()
  pl_bin5.2 <- prioritylasso(X = over_matrix, Y = y3,weights = agirlik, 
                             family = "binomial", type.measure = "auc",
                             blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                             block1.penalization = TRUE, lambda.type = "lambda.min",
                             standardize = TRUE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin5.2 <- as.matrix(transformed_test)
  ypredscore5.2 <-predict(object = pl_bin5.2, newdata = newdata_bin5.2, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  
  pred5.2 <- prediction(ypredscore5.2,tvsd_test3$status)
  pred5.2@predictions <- list(unlist(pred5.2@predictions) + probVector)
  perf5.2 <- performance(pred5.2, measure = "tpr", x.measure = "fpr")
  plot(perf5.2,col='red')
  abline(a=0,b=1)
  auc5.2 <- performance(pred5.2,measure="auc")@y.values
  auc5.2
  timer_stop_tune <- Sys.time()
  auc5.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 5.2 is done")
  
  #### 6.2 model ####
  timer_start_tune <- Sys.time()
  pl_bin6.2 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
                             family = "binomial", type.measure = "auc",
                             # blocks = list(block1=1:ncol(over_matrix)),
                             blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                             block1.penalization = TRUE, lambda.type = "lambda.min",
                             standardize = FALSE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin6.2 <- as.matrix(transformed_test)
  ypredscore6.2 <-predict(object = pl_bin6.2, newdata = newdata_bin6.2, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  
  pred6.2 <- prediction(ypredscore6.2,tvsd_test3$status)
  pred6.2@predictions <- list(unlist(pred6.2@predictions) + probVector)
  perf6.2 <- performance(pred6.2, measure = "tpr", x.measure = "fpr")
  plot(perf6.2,col='red')
  abline(a=0,b=1)
  auc6.2 <- performance(pred6.2,measure="auc")@y.values
  auc6.2
  timer_stop_tune <- Sys.time()
  auc6.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 6.2 is done")
  
  #### 7.2 model ####
  timer_start_tune <- Sys.time()
  pl_bin7.2 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
                             family = "binomial", type.measure = "auc",
                             blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                             block1.penalization = TRUE, lambda.type = "lambda.1se",
                             standardize = TRUE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin7.2 <- as.matrix(transformed_test)
  ypredscore7.2 <-predict(object = pl_bin7.2, newdata = newdata_bin7.2, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  pred7.2 <- prediction(ypredscore7.2,tvsd_test3$status)
  pred7.2@predictions <- list(unlist(pred7.2@predictions) + probVector)
  perf7.2 <- performance(pred7.2, measure = "tpr", x.measure = "fpr")
  plot(perf7.2,col='red')
  abline(a=0,b=1)
  auc7.2 <- performance(pred7.2,measure="auc")@y.values
  auc7.2
  timer_stop_tune <- Sys.time()
  auc7.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 7.2 is done")
  
  #### 8.2 model ####
  timer_start_tune <- Sys.time()
  pl_bin8.2 <- prioritylasso(X = over_matrix, Y = y3, weights = agirlik, 
                             family = "binomial", type.measure = "auc",
                             blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                             block1.penalization = TRUE, lambda.type = "lambda.1se",
                             standardize = FALSE, nfolds = 5, cvoffset = TRUE, cvoffsetnfolds = 10)
  newdata_bin8.2 <- as.matrix(transformed_test)
  ypredscore8.2 <-predict(object = pl_bin8.2, newdata = newdata_bin8.2, type = "response")
  
  #ypredscore <-ypredscore+(1/2)
  
  
  
  pred8.2 <- prediction(ypredscore8.2,tvsd_test3$status)
  pred8.2@predictions <- list(unlist(pred8.2@predictions) + probVector) # add mean of risk sets to the results
  perf8.2 <- performance(pred8.2, measure = "tpr", x.measure = "fpr")
  plot(perf8.2,col='red')
  abline(a=0,b=1)
  auc8.2 <- performance(pred8.2,measure="auc")@y.values
  auc8.2
  timer_stop_tune <- Sys.time()
  auc8.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 8.2 is done")
  
  print("Priority Lasso all done!")
  
  #### Ipf-lasso ####
  print("Starting IPF Lasso")
  y3 <-as.numeric(over$prediction_column)
  y3[y3==1] <- 0
  y3[y3==2] <- 1
  
  ##### model 1 #####
  timer_start_tune <- Sys.time()
  model1<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=FALSE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin1 <- as.matrix(transformed_test)
  ypredscore1 <- ipflasso.predict(object=model1,Xtest=newdata_bin1)
  ypredscore1$probabilitiestest <- ypredscore1$probabilitiestest + probVector
  auc_ipf1 <- my.auc(ypredscore1$linpredtest, tvsd_test3$status)
  auc_ipf1 
  timer_stop_tune <- Sys.time()
  auc_ipf1_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 1 is done")
  
  ##### model 1.2 #####
  timer_start_tune <- Sys.time()
  model1.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=FALSE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin1.2 <- as.matrix(transformed_test)
  ypredscore1.2 <- ipflasso.predict(object=model1.2,Xtest=newdata_bin1.2)
  ypredscore1.2$probabilitiestest <- ypredscore1.2$probabilitiestest + probVector
  auc_ipf1.2 <- my.auc(ypredscore1.2$linpredtest, tvsd_test3$status)
  auc_ipf1.2 
  timer_stop_tune <- Sys.time()
  auc_ipf1.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 1.2 is done")
  
  ##### model 2 #####
  timer_start_tune <- Sys.time()
  model2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin2 <- as.matrix(transformed_test)
  ypredscore2 <- ipflasso.predict(object=model2,Xtest=newdata_bin2)
  ypredscore2$probabilitiestest <- ypredscore2$probabilitiestest + probVector
  auc_ipf2 <- my.auc(ypredscore2$linpredtest, tvsd_test3$status)
  auc_ipf2 
  timer_stop_tune <- Sys.time()
  auc_ipf2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 2 is done")
  
  ##### model 2.2 #####
  timer_start_tune <- Sys.time()
  model2.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin2.2 <- as.matrix(transformed_test)
  ypredscore2.2 <- ipflasso.predict(object=model2.2,Xtest=newdata_bin2.2)
  ypredscore2.2$probabilitiestest <- ypredscore2.2$probabilitiestest + probVector
  auc_ipf2.2 <- my.auc(ypredscore2.2$linpredtest, tvsd_test3$status)
  auc_ipf2.2 
  timer_stop_tune <- Sys.time()
  auc_ipf2.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 2.2 is done")
  
  ##### model 3 #####
  timer_start_tune <- Sys.time()
  model3<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(1,1),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin3 <- as.matrix(transformed_test)
  ypredscore3 <- ipflasso.predict(object=model3,Xtest=newdata_bin3)
  ypredscore3$probabilitiestest <- ypredscore3$probabilitiestest + probVector
  auc_ipf3 <- my.auc(ypredscore3$linpredtest, tvsd_test3$status)
  auc_ipf3 
  timer_stop_tune <- Sys.time()
  auc_ipf3_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 3 is done")
  
  ##### model 3.2 #####
  timer_start_tune <- Sys.time()
  model3.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(1,1),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin3.2 <- as.matrix(transformed_test)
  ypredscore3.2 <- ipflasso.predict(object=model3.2,Xtest=newdata_bin3.2)
  ypredscore3.2$probabilitiestest <- ypredscore3.2$probabilitiestest + probVector
  auc_ipf3.2 <- my.auc(ypredscore3.2$linpredtest, tvsd_test3$status)
  auc_ipf3.2 
  timer_stop_tune <- Sys.time()
  auc_ipf3.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 3.2 is done")
  
  ##### model 4 #####
  timer_start_tune <- Sys.time()
  model4<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik,family="binomial",standardize=TRUE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(2,1),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin4 <- as.matrix(transformed_test)
  ypredscore4 <- ipflasso.predict(object=model4,Xtest=newdata_bin4)
  ypredscore4$probabilitiestest <- ypredscore4$probabilitiestest + probVector
  auc_ipf4 <- my.auc(ypredscore4$linpredtest, tvsd_test3$status)
  auc_ipf4 
  timer_stop_tune <- Sys.time()
  auc_ipf4_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 4 is done")
  
  ##### model 4.2 #####
  timer_start_tune <- Sys.time()
  model4.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik,family="binomial",standardize=TRUE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(2,1),nfolds=5,ncv=10,type.measure="auc", alpha=0)
  newdata_bin4.2 <- as.matrix(transformed_test)
  ypredscore4.2 <- ipflasso.predict(object=model4.2,Xtest=newdata_bin4.2)
  ypredscore4.2$probabilitiestest <- ypredscore4.2$probabilitiestest + probVector
  auc_ipf4.2 <- my.auc(ypredscore4.2$linpredtest, tvsd_test3$status)
  auc_ipf4.2 
  timer_stop_tune <- Sys.time()
  auc_ipf4.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 4.2 is done")
  
  ##### model 5 #####
  # model5<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
  #                        blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                        pf=c(1,2),nfolds=5,ncv=10,type.measure="class", alpha=0)
  # newdata_bin5 <- as.matrix(transformed_test)
  # ypredscore5 <- ipflasso.predict(object=model5,Xtest=newdata_bin5)
  # ypredscore5$probabilitiestest <- ypredscore5$probabilitiestest + probVector
  # auc_ipf5 <- my.auc(ypredscore5$linpredtest, tvsd_test3$status)
  # auc_ipf5 
  # 
  # 
  # print("Model 5 is done")
  ##### model 6 #####
  timer_start_tune <- Sys.time()
  model6<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=FALSE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin6 <- as.matrix(transformed_test)
  ypredscore6 <- ipflasso.predict(object=model6,Xtest=newdata_bin6)
  ypredscore6$probabilitiestest <- ypredscore6$probabilitiestest + probVector
  auc_ipf6 <- my.auc(ypredscore6$linpredtest, tvsd_test3$status)
  auc_ipf6 
  timer_stop_tune <- Sys.time()
  auc_ipf6_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 6 is done")
  
  ##### model 6.2 #####
  timer_start_tune <- Sys.time()
  model6.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=FALSE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin6.2 <- as.matrix(transformed_test)
  ypredscore6.2 <- ipflasso.predict(object=model6.2,Xtest=newdata_bin6.2)
  ypredscore6.2$probabilitiestest <- ypredscore6.2$probabilitiestest + probVector
  auc_ipf6.2 <- my.auc(ypredscore6.2$linpredtest, tvsd_test3$status)
  auc_ipf6.2 
  timer_stop_tune <- Sys.time()
  auc_ipf6.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 6.2 is done")
  
  ##### model 7 #####
  timer_start_tune <- Sys.time()
  model7<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin7 <- as.matrix(transformed_test)
  ypredscore7 <- ipflasso.predict(object=model7,Xtest=newdata_bin7)
  ypredscore7$probabilitiestest <- ypredscore7$probabilitiestest + probVector
  auc_ipf7 <- my.auc(ypredscore7$linpredtest, tvsd_test3$status)
  auc_ipf7 
  timer_stop_tune <- Sys.time()
  auc_ipf7_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 7 is done")
  
  ##### model 7.2 #####
  timer_start_tune <- Sys.time()
  model7.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin7.2 <- as.matrix(transformed_test)
  ypredscore7.2 <- ipflasso.predict(object=model7.2,Xtest=newdata_bin7.2)
  ypredscore7.2$probabilitiestest <- ypredscore7.2$probabilitiestest + probVector
  auc_ipf7.2 <- my.auc(ypredscore7.2$linpredtest, tvsd_test3$status)
  auc_ipf7.2 
  timer_stop_tune <- Sys.time()
  auc_ipf7.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 7.2 is done")
  
  ##### model 8 #####
  timer_start_tune <- Sys.time()
  model8<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(1,1),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin8 <- as.matrix(transformed_test)
  ypredscore8 <- ipflasso.predict(object=model8,Xtest=newdata_bin8)
  ypredscore8$probabilitiestest <- ypredscore8$probabilitiestest + probVector
  auc_ipf8 <- my.auc(ypredscore8$linpredtest, tvsd_test3$status)
  auc_ipf8 
  timer_stop_tune <- Sys.time()
  auc_ipf8_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 8 is done")
  
  ##### model 8.2 #####
  timer_start_tune <- Sys.time()
  model8.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(1,1),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin8.2 <- as.matrix(transformed_test)
  ypredscore8.2 <- ipflasso.predict(object=model8.2,Xtest=newdata_bin8.2)
  ypredscore8.2$probabilitiestest <- ypredscore8.2$probabilitiestest + probVector
  auc_ipf8.2 <- my.auc(ypredscore8.2$linpredtest, tvsd_test3$status)
  auc_ipf8.2 
  timer_stop_tune <- Sys.time()
  auc_ipf8.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 8.2 is done")
  
  ##### model 9 #####
  timer_start_tune <- Sys.time()
  model9<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                         pf=c(2,1),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin9 <- as.matrix(transformed_test)
  ypredscore9 <- ipflasso.predict(object=model9,Xtest=newdata_bin9)
  ypredscore9$probabilitiestest <- ypredscore9$probabilitiestest + probVector
  auc_ipf9 <- my.auc(ypredscore9$linpredtest, tvsd_test3$status)
  auc_ipf9 
  timer_stop_tune <- Sys.time()
  auc_ipf9_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 9 is done")
  
  ##### model 9.2 #####
  timer_start_tune <- Sys.time()
  model9.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                           blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                           pf=c(2,1),nfolds=5,ncv=10,type.measure="auc", alpha=1)
  newdata_bin9.2 <- as.matrix(transformed_test)
  ypredscore9.2 <- ipflasso.predict(object=model9.2,Xtest=newdata_bin9.2)
  ypredscore9.2$probabilitiestest <- ypredscore9.2$probabilitiestest + probVector
  auc_ipf9.2 <- my.auc(ypredscore9.2$linpredtest, tvsd_test3$status)
  auc_ipf9.2 
  timer_stop_tune <- Sys.time()
  auc_ipf9.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  print("Model 9.2 is done")
  
  ##### model 10 #####
  # model10<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
  #                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                         pf=c(1,2),nfolds=5,ncv=10,type.measure="class", alpha=1)
  # newdata_bin10 <- as.matrix(transformed_test)
  # ypredscore10 <- ipflasso.predict(object=model10,Xtest=newdata_bin10)
  # ypredscore10$probabilitiestest <- ypredscore10$probabilitiestest + probVector
  # auc_ipf10 <- my.auc(ypredscore10$linpredtest, tvsd_test3$status)
  # auc_ipf10 
  # 
  # print("Model 10 is done")
  ##### model 11 #####
  timer_start_tune <- Sys.time()
  model11<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                          blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
                          pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=0.5)
  newdata_bin11 <- as.matrix(transformed_test)
  ypredscore11 <- ipflasso.predict(object=model11,Xtest=newdata_bin11)
  ypredscore11$probabilitiestest <- ypredscore11$probabilitiestest + probVector
  auc_ipf11 <- my.auc(ypredscore11$linpredtest, tvsd_test3$status)
  auc_ipf11 
  timer_stop_tune <- Sys.time()
  auc_ipf11_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 11 is done")
  
  ##### model 11.2 #####
  timer_start_tune <- Sys.time()
  model11.2<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
                            blocks = list(block1=(event_count+1):ncol(over_matrix), block2=1:event_count),
                            pf=c(1,2),nfolds=5,ncv=10,type.measure="auc", alpha=0.5)
  newdata_bin11.2 <- as.matrix(transformed_test)
  ypredscore11.2 <- ipflasso.predict(object=model11.2,Xtest=newdata_bin11.2)
  ypredscore11.2$probabilitiestest <- ypredscore11.2$probabilitiestest + probVector
  auc_ipf11.2 <- my.auc(ypredscore11.2$linpredtest, tvsd_test3$status)
  auc_ipf11.2 
  timer_stop_tune <- Sys.time()
  auc_ipf11.2_time <- difftime(timer_stop_tune, timer_start_tune, units = "mins")
  
  print("Model 11.2 is done")
  
  ##### model 12 #####
  # model12<-cvr.ipflasso_n(X=over_matrix,Y=y3, weights = agirlik, family="binomial",standardize=TRUE,
  #                         blocks = list(block1=1:event_count, block2=(event_count+1):ncol(over_matrix)),
  #                         pf=c(1,2),nfolds=5,ncv=10,type.measure="class", alpha=0.5)
  # newdata_bin12 <- as.matrix(transformed_test)
  # ypredscore12 <- ipflasso.predict(object=model12,Xtest=newdata_bin12)
  # ypredscore12$probabilitiestest <- ypredscore12$probabilitiestest + probVector
  # auc_ipf12 <- my.auc(ypredscore12$linpredtest, tvsd_test3$status)
  # auc_ipf12 
  
  #### lasso 0.6764706- alpha=0 0.7352941 - alpha=1 0.745098
  #### svm 0.8235294 - alpha=0 0.7058824 - alpha=1 0.7843137
  #### 
  
  
  
  # print("Model 12 is done")
  print("All done")
  #### MERGING RESULTS ####
  iterRes_prioLasso <- c(auc5, auc6, auc7, auc8, auc5.2, auc6.2, auc7.2, auc8.2)
  iterRes_ipf <- c(auc_ipf1, auc_ipf2, auc_ipf3, auc_ipf4, auc_ipf6, auc_ipf7, auc_ipf8, auc_ipf9, auc_ipf11,
                   auc_ipf1.2, auc_ipf2.2, auc_ipf3.2, auc_ipf4.2, auc_ipf6.2, auc_ipf7.2, auc_ipf8.2, auc_ipf9.2, auc_ipf11.2)
  iterRes <- c(iterRes_prioLasso, iterRes_ipf)
  results_df <- rbind(results_df, iterRes)
  colnames(results_df) <- c("priolasso5", "priolasso6", "priolasso7", "priolasso8",
                            "priolasso5.2", "priolasso6.2", "priolasso7.2", "priolasso8.2",
                            "ipflasso1", "ipflasso2", "ipflasso3", "ipflasso4",
                            "ipflasso6", "ipflasso7", "ipflasso8",
                            "ipflasso9", "ipflasso11",
                            "ipflasso1.2", "ipflasso2.2", "ipflasso3.2", "ipflasso4.2",
                            "ipflasso6.2", "ipflasso7.2", "ipflasso8.2",
                            "ipflasso9.2", "ipflasso11.2")
  fwrite(results_df, paste("MLSeqSurv_voomResults/", cancerType, "/iter", index, ".csv", sep = ""))
  
  results_df <- rbind(results_df, iterRes)
  
  probVector <- get_matching_probabilities(summ_stack, test_final)
  #priolasso_selected_features <- getFeaturesFromPrioLasso(pl_bin1)
  ipf_selected_features <- getFeaturesFromIPF(model1)
  
  # pl_bin1_features <- getFeaturesFromPrioLasso(pl_bin1)
  # pl_bin2_features <- getFeaturesFromPrioLasso(pl_bin2)
  # pl_bin3_features <- getFeaturesFromPrioLasso(pl_bin3)
  # pl_bin4_features <- getFeaturesFromPrioLasso(pl_bin4)
  pl_bin5_features <- getFeaturesFromPrioLasso(pl_bin5)
  pl_bin6_features <- getFeaturesFromPrioLasso(pl_bin6)
  pl_bin7_features <- getFeaturesFromPrioLasso(pl_bin7)
  pl_bin8_features <- getFeaturesFromPrioLasso(pl_bin8)
  pl_bin5.2_features <- getFeaturesFromPrioLasso(pl_bin5.2)
  pl_bin6.2_features <- getFeaturesFromPrioLasso(pl_bin6.2)
  pl_bin7.2_features <- getFeaturesFromPrioLasso(pl_bin7.2)
  pl_bin8.2_features <- getFeaturesFromPrioLasso(pl_bin8.2)
  #
  pl_features_matrix <- cbind(#pl_bin1_features,pl_bin2_features,pl_bin3_features,pl_bin4_features,
    pl_bin5_features,pl_bin6_features,pl_bin7_features,pl_bin8_features,
    pl_bin5.2_features,pl_bin6.2_features,pl_bin7.2_features,pl_bin8.2_features)
  pl_features_df <- as.data.frame(pl_features_matrix)
  #
  ipf_model1_features <- getFeaturesFromIPF(model1)
  ipf_model2_features <- getFeaturesFromIPF(model2)
  ipf_model3_features <- getFeaturesFromIPF(model3)
  ipf_model4_features <- getFeaturesFromIPF(model4)
  # ipf_model5_features <- getFeaturesFromIPF(model5)
  ipf_model6_features <- getFeaturesFromIPF(model6)
  ipf_model7_features <- getFeaturesFromIPF(model7)
  ipf_model8_features <- getFeaturesFromIPF(model8)
  ipf_model9_features <- getFeaturesFromIPF(model9)
  # ipf_model10_features <- getFeaturesFromIPF(model10)
  ipf_model11_features <- getFeaturesFromIPF(model11)
  # ipf_model12_features <- getFeaturesFromIPF(model12)
  
  ipf_model1.2_features <- getFeaturesFromIPF(model1.2)
  ipf_model2.2_features <- getFeaturesFromIPF(model2.2)
  ipf_model3.2_features <- getFeaturesFromIPF(model3.2)
  ipf_model4.2_features <- getFeaturesFromIPF(model4.2)
  # ipf_model5_features <- getFeaturesFromIPF(model5)
  ipf_model6.2_features <- getFeaturesFromIPF(model6.2)
  ipf_model7.2_features <- getFeaturesFromIPF(model7.2)
  ipf_model8.2_features <- getFeaturesFromIPF(model8.2)
  ipf_model9.2_features <- getFeaturesFromIPF(model9.2)
  # ipf_model10_features <- getFeaturesFromIPF(model10)
  ipf_model11.2_features <- getFeaturesFromIPF(model11.2)
  # ipf_model12_features <- getFeaturesFromIPF(model12)
  
  #
  maxLength <- ncol(over_matrix)
  length(ipf_model1_features) <- maxLength
  length(ipf_model2_features) <- maxLength
  length(ipf_model3_features) <- maxLength
  length(ipf_model4_features) <- maxLength
  # length(ipf_model5_features) <- maxLength
  length(ipf_model6_features) <- maxLength
  length(ipf_model7_features) <- maxLength
  length(ipf_model8_features) <- maxLength
  length(ipf_model9_features) <- maxLength
  # length(ipf_model10_features) <- maxLength
  length(ipf_model11_features) <- maxLength
  # length(ipf_model12_features) <- maxLength
  
  length(ipf_model1.2_features) <- maxLength
  length(ipf_model2.2_features) <- maxLength
  length(ipf_model3.2_features) <- maxLength
  length(ipf_model4.2_features) <- maxLength
  # length(ipf_model5_features) <- maxLength
  length(ipf_model6.2_features) <- maxLength
  length(ipf_model7.2_features) <- maxLength
  length(ipf_model8.2_features) <- maxLength
  length(ipf_model9.2_features) <- maxLength
  # length(ipf_model10_features) <- maxLength
  length(ipf_model11.2_features) <- maxLength
  # length(ipf_model12_features) <- maxLength
  
  #
  #
  ipf_features_matrix <- cbind(ipf_model1_features,ipf_model2_features,ipf_model3_features,
                               ipf_model4_features,ipf_model6_features,
                               ipf_model7_features,ipf_model8_features,ipf_model9_features,
                               ipf_model11_features,
                               ipf_model1.2_features,ipf_model2.2_features,ipf_model3.2_features,
                               ipf_model4.2_features,ipf_model6.2_features,
                               ipf_model7.2_features,ipf_model8.2_features,ipf_model9.2_features,
                               ipf_model11.2_features)
  ipf_features_df <- as.data.frame(ipf_features_matrix)
  
  fwrite(pl_features_df, paste0("MLSeqSurv_voomSelectedFeatures/", cancerType, "/Priority_Lasso_Features",index,".csv"))
  fwrite(ipf_features_df, paste0("MLSeqSurv_voomSelectedFeatures/", cancerType, "/IPF_Lasso_Features",index,".csv"))
  
  #### TIME RESULTS ####
  time_prioLasso <- c(auc5_time, auc6_time, auc7_time, auc8_time, auc5.2_time, auc6.2_time, auc7.2_time, auc8.2_time)
  time_ipf <- c(auc_ipf1_time, auc_ipf2_time, auc_ipf3_time, auc_ipf4_time, auc_ipf6_time, auc_ipf7_time, auc_ipf8_time, auc_ipf9_time, auc_ipf11_time,
                auc_ipf1.2_time, auc_ipf2.2_time, auc_ipf3.2_time, auc_ipf4.2_time, auc_ipf6.2_time, auc_ipf7.2_time, auc_ipf8.2_time, 
                auc_ipf9.2_time, auc_ipf11.2_time)
  time_res <- c(time_prioLasso, time_ipf)
  results_time <- rbind(results_time, time_res)
  colnames(results_time) <- c("priolasso5_time", "priolasso6_time", "priolasso7_time", "priolasso8_time",
                              "priolasso5.2_time", "priolasso6.2_time", "priolasso7.2_time", "priolasso8.2_time",
                              "ipflasso1_time", "ipflasso2_time", "ipflasso3_time", "ipflasso4_time",
                              "ipflasso6_time", "ipflasso7_time", "ipflasso8_time",
                              "ipflasso9_time", "ipflasso11_time",
                              "ipflasso1.2_time", "ipflasso2.2_time", "ipflasso3.2_time", "ipflasso4.2_time",
                              "ipflasso6.2_time", "ipflasso7.2_time", "ipflasso8.2_time",
                              "ipflasso9.2_time", "ipflasso11.2_time")
  fwrite(results_time, paste("MLSeqSurv_voomResults_time/", cancerType, "/iter", index, ".csv", sep = ""))
  
  
  
  
}


#### Get feature counts and write them to csv files ####
selectedFeatureCount_IPF_df <- data_frame()
selectedFeatureCount_PrioLasso_df <- data_frame()
for (index in 1:30){
  
  selectedFeaturesIPF <- fread(file = paste("MLSeqSurv_voomSelectedFeatures/",cancerType, "/IPF_Lasso_Features", index, ".csv", sep = ""))
  selectedFeaturesPrioLasso <- fread(file = paste("MLSeqSurv_voomSelectedFeatures/",cancerType, "/Priority_Lasso_Features", index, ".csv", sep = ""))
  
  
  featureCountIPF <- colSums(!is.na(selectedFeaturesIPF)) # shows the amount of features selected by an algorithm for a specific iteration
  featureCountPrioLasso <- colSums(!is.na(selectedFeaturesPrioLasso))
  
  selectedFeatureCount_IPF_df <- rbind(selectedFeatureCount_IPF_df, featureCountIPF)
  selectedFeatureCount_PrioLasso_df <- rbind(selectedFeatureCount_PrioLasso_df, featureCountPrioLasso)
  
  
  
}

colnames(selectedFeatureCount_IPF_df) <- names(featureCountIPF)
colnames(selectedFeatureCount_PrioLasso_df) <- names(selectedFeaturesPrioLasso)

fwrite(selectedFeatureCount_IPF_df, paste("MLSeqSurv_voomResults/", cancerType, "/FeatureCountsIPF.csv", sep = ""))
fwrite(selectedFeatureCount_PrioLasso_df, paste("MLSeqSurv_voomResults/", cancerType, "/FeatureCountsPriorityLasso.csv", sep = ""))
