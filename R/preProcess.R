# These functions need to be documented and relabeled appropriately.
# Will there be variance filtering? If there is going to be variance filtering, what is the cutoff
#### Functions necessary for Preproecessing ####

arrayWeights <- function(object, design=NULL, weights=NULL, var.design=NULL, var.group=NULL, prior.n=10,
                         method="none", maxiter=50L, tol=1e-5, trace=FALSE)
  #	Estimate array quality weights.
  #
  #	Created by Matt Ritchie 7 Feb 2005.
  #	Gordon Smyth simplified argument checking to use getEAWP, 9 Mar 2008.
  #	Cynthia Liu added var.design argument 22 Sep 2014.
  #	Rewrite by Gordon Smyth 12 Feb 2019.
  #	Last modified 14 Feb 2019.
{
  #	Check object
  y <- getEAWP(object)
  E <- y$exprs
  ngenes <- nrow(E)
  narrays <- ncol(E)

  #	Initial values for array weights
  w <- rep_len(1,narrays)
  names(w) <- colnames(E)

  #	Require at least 2 rows for estimates to be useful
  if(ngenes < 2L) return(w)

  #	Check design
  if(is.null(design)) {
    design <- matrix(1,narrays,1)
    p <- 1L
  } else {
    design <- as.matrix(design)
    if(mode(design) != "numeric") stop("design must be a numeric matrix")
    QR <- qr(design)
    p <- QR$rank

    #		If not full rank, remove superfluous columns
    if(p < ncol(design)) design <- design[,QR$pivot[1:p],drop=FALSE]
  }

  #	Require at least 2 residual df.
  if(narrays - p < 2L) return(w)

  #	Check weights
  if(is.null(weights)) weights <- y$weights
  if(!is.null(weights)) {
    weights <- as.matrix(weights)
    if(nrow(weights)!=ngenes || ncol(weights)!=narrays) stop("dimensions of weights should match those of expression object")
    r <- range(weights)
    if(is.na(r[1])) stop("NA weights not allowed")
    if(r[1] < 0) stop("Negative weights not allowed")
    if(is.infinite(r[2])) stop("Infinite weights not allowed")
    if(r[1]==0) {
      E[weights==0] <- NA
      weights[weights==0] <- 1
    }
  }

  #	var.group takes precedence over var.design
  if(!is.null(var.group)) {
    var.group <- droplevels(as.factor(var.group))
    if(length(var.group) != narrays) stop("var.group has wrong length")
    if(nlevels(var.group) < 2L) stop("Need at least two variance groups")
    contrasts(var.group) <- contr.sum(levels(var.group))
    var.design <- model.matrix(~var.group)
    var.design <- var.design[,-1,drop=FALSE]
  }

  #	Setup variance design matrix
  #	First column must be an intercept and other columns must add to zero
  if(is.null(var.design)) {
    Z2 <- contr.sum(narrays)
  } else {
    Z2 <- var.design
    Z2 <- t(t(Z2) - colMeans(Z2))
    QR <- qr(Z2)
    Z2 <- Z2[,QR$pivot[1:QR$rank],drop=FALSE]
  }

  #	Detect NA and infinite values. Convert latter into NAs.
  r <- range(E)
  if(!all(is.finite(r))) {
    E[is.infinite(E)] <- NA
    HasNA <- TRUE
  } else {
    HasNA <- FALSE
  }

  #	Check method
  method <- match.arg(method,c("auto","genebygene","reml"))
  if(method=="auto")
    if(HasNA || !is.null(weights))
      method <- "genebygene"
  else
    method <- "reml"

  if(method=="genebygene")
    return(.arrayWeightsGeneByGene(E, design=design, weights=weights, var.design=Z2, prior.n=prior.n, trace=trace))

  if(method=="reml") {
    if(HasNA) {
      iNA <- is.na(rowMeans(E))
      message("removing ",sum(iNA)," rows with missing or infinite values")
      E <- E[!iNA,]
      if(!is.null(weights)) weights <- weights[!iNA,]
      if(nrow(E) < 2L) return(w)
    }
    if(is.null(weights)) {
      return(.arrayWeightsREML(E, design=design, var.design=Z2, prior.n=prior.n, maxiter=maxiter, tol=tol, trace=trace))
    } else {
      return(.arrayWeightsPrWtsREML(E, design=design, weights=weights, var.design=Z2, prior.n=prior.n, maxiter=maxiter, tol=tol, trace=trace))
    }
  }

}

.arrayWeightsGeneByGene <- function(E, design=NULL, weights=NULL, var.design=NULL, prior.n=10, trace=FALSE)
  #	Estimate array variances via gene-by-gene update algorithm
  #	Created by Matt Ritchie 7 Feb 2005.
  #	Cynthia Liu added var.design argument 22 Sep 2014.
  #	Fixes and speedups by Gordon Smyth 15 Feb 2019, 28 Apr 2020.
  #	Last modified 28 Feb 2020.
{
  ngenes <- nrow(E)
  narrays <- ncol(E)
  if(is.null(design)) design <- matrix(1,narrays,1)
  nparams <- ncol(design)

  #	Columns of var.design should sum to zero
  if(is.null(var.design)) {
    Z2 <- contr.sum(narrays)
  } else {
    Z2 <- var.design
  }
  Z <- cbind(1,Z2)

  #	Intialise array gammas to zero (with prior weight of prior.n genes having leverage=0)
  ngam <- ncol(Z2)
  gam <- rep_len(0, ngam)
  aw <- rep_len(1,narrays)
  info2 <- prior.n*crossprod(Z2)

  #	If requested, progess will be output 10 times at equal intervals
  if(trace) {
    cat("gene range(w)\n")
    ReportInterval <- pmax(as.integer(ngenes/10),1L)
  }

  #	Step progressive algorithm once through all genes
  Zero <- rep_len(0,narrays)
  One <- rep_len(1,narrays)
  for(i in 1:ngenes) {
    if(is.null(weights)) {
      w <- aw
    } else {
      w <- aw*weights[i,]
    }
    y <- E[i,]
    if(anyNA(y)) {
      obs <- is.finite(y)
      nobs <- sum(obs)
      if(nobs <= 2L) next
      X <- design[obs, , drop = FALSE]
      y <- y[obs]
      w <- w[obs]
      fit <- lm.wfit(X, y, w)
      if(fit$df.residual < 2L) next
      h1 <- d <- Zero
      d[obs] <- w*fit$residuals^2
      h1[obs] <- 1-hat(fit$qr)
    } else {
      fit <- lm.wfit(design, y, w)
      d <- w*fit$residuals^2
      h1 <- 1-hat(fit$qr)
    }
    s2 <- mean(fit$effects[-(1:fit$rank)]^2)
    if(s2 < 1e-15) next
    info <- crossprod(Z, h1*Z)
    info2 <- info2 + info[-1,-1,drop=FALSE] - (info[-1,1,drop=FALSE]/info[1,1]) %*% info[1,-1,drop=FALSE]
    z <- d/s2 - h1
    dl <- crossprod(Z2, z)
    gam <- gam + solve(info2, dl)
    aw <- drop(exp(Z2 %*% (-gam)))

    #		Progress output
    if(trace && (i %% ReportInterval==0L)) cat(i,range(aw),"\n")
  }

  aw
}

#Calculation of deseq normalization factors necessary for calcNormFactorsGSD function.
calcFactorRLEGSD <- function(data.train, data.test, lib.size, lib.size.test){
  gm <- exp(rowMeans(log(data.train)))
  f = apply(data.train, 2, function(u) median((u/gm)[gm > 0]))
  f.test = apply(data.test, 2, function(u) median((u/gm)[gm > 0]))
  f = f / lib.size
  f.test = f.test / lib.size.test
  deseqsizefactors = list(f, f.test)
  return(deseqsizefactors)
}

calcNormFactorsGSD <- function(data.train, data.test, lib.size = NULL, method = c("TMM", "deseq", "none"), refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05,
                               doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, ...){
  x <- as.matrix(data.train)
  xtest <- as.matrix(data.test)
  #if (any(is.na(x)||is.na(xtest)))
  #  stop("NAs not permitted")
  if (is.null(lib.size))
    lib.size <- colSums(x)

  lib.size.test <- colSums(xtest)
  method <- match.arg(method)
  allzero <- rowSums(x > 0) == 0
  if (any(allzero)){
    x <- x[!allzero, , drop = FALSE]
  }

  xtest <- xtest[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 1){
    method = "none"
  }

  if (method == "TMM"){
    f75 <- calcFactorQuantileGSD(data = x, lib.size = lib.size, p = 0.75)
    f75.test <- calcFactorQuantileGSD(data = xtest, lib.size = lib.size.test, p = 0.75)

    refColumn <- which.min(abs(f75 - mean(f75)))
    f <- rep(NA, ncol(x))
    f.test <- rep(NA, ncol(xtest))
    for (i in 1:ncol(x)) f[i] <- calcFactorWeightedGSD(obs = x[,i], ref = x[, refColumn], libsize.obs = lib.size[i],
                                                       libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                       sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
    for (i in 1:ncol(xtest)) f.test[i] <- calcFactorWeightedGSD(obs = xtest[,i], ref = x[, refColumn], libsize.obs = lib.size.test[i],
                                                                libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                                sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
    normf = list(f,f.test)
  } else if (method == "deseq"){
    normf = calcFactorRLEGSD(data.train = x, data.test = xtest, lib.size = lib.size, lib.size.test = lib.size.test)#/lib.size
  } else {
    normf = list(rep(1, ncol(x)), rep(1, ncol(xtest)))
  }

  names(normf) = c("train", "test")

  f = as.numeric(normf[[1]]) / (exp(mean(log(normf[[1]]))))
  f.test = as.numeric(normf[[2]]) / (exp(mean(log(normf[[1]]))))
  normf2 = list(f, f.test, lib.size, lib.size.test)
  names(normf2) = c("TrainNormFactor","TestNormFactor","TrainLibSize","TestLibSize")
  return(normf2)
}


voom_lasso <- function(data.train, data.test, design=NULL,lib.size=NULL,normalize.method="none",block=NULL,correlation=NULL,weights=NULL,
                     span=0.5,plot=FALSE,save.plot=FALSE)
  #	Linear modelling of count data with mean-variance modelling at the observation level.
  #	Creates an EList object for entry to lmFit() etc in the limma pipeline.
  #	Gordon Smyth and Charity Law
  #	Created 22 June 2011.  Last modified 23 January 2020.
{
  ##### Near-zero variance filtering ####
  tr_dge <- DGEList(counts=data.train)
  keep <- filterByExpr(tr_dge, design=NULL,group = NULL, lib.size = NULL)
  ts_dge <- DGEList(counts=data.test)



  out <- list()

  # #	Extract counts from known data objects
  # if(is(counts,"DGEList")) {
  #   out$genes <- counts$genes
  #   out$targets <- counts$samples
  #   if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
  #   if(is.null(lib.size)) lib.size <- counts$samples$lib.size*counts$samples$norm.factors
  #   counts <- counts$counts
  # } else {
  #   isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
  #   if(isExpressionSet) {
  #     if(length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
  #     if(length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
  #     counts <- Biobase::exprs(counts)
  #   } else {
  #     counts <- as.matrix(counts)
  #   }
  # }
  #
  # #	Check counts
  # n <- nrow(counts)
  # if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")
  # m <- min(counts)
  # if(is.na(m)) stop("NA counts not allowed")
  # if(m < 0) stop("Negative counts now allowed")
  #
  #	Check design
  # if(is.null(design)) {
  #   design <- matrix(1,ncol(counts),1)
  #   rownames(design) <- colnames(counts)
  #   colnames(design) <- "GrandMean"
  # }
  #
  # #	Check lib.size
  # if(is.null(lib.size)) lib.size <- colSums(counts)
  #
  NormFactors = calcNormFactorsGSD(data.train = data.train, data.test = data.test, method = "deseq")
  TrainNormFactor = NormFactors$TrainNormFactor
  TestNormFactor = NormFactors$TestNormFactor
  TrainLibSize = NormFactors$TrainLibSize
  TestLibSize = NormFactors$TestLibSize
  lib.size.tr = TrainNormFactor * TrainLibSize
  lib.size.ts = TestNormFactor * TestLibSize

  #design.tr = model.matrix(~group)
  design.tr = matrix(1, ncol(data.train), 1)
  rownames(design.tr) = colnames(data.train)

  design.ts <- matrix(1, ncol(data.test), 1)
  rownames(design.ts) <- colnames(data.test)
  colnames(design.ts) <- "GrandMean"

  #	Fit linear model to log2-counts-per-million
  y.tr <- t(log2(t(data.train + 0.5)/(lib.size.tr + 1) * 1e+06))
  y.ts <- t(log2(t(data.test + 0.5)/(lib.size.ts + 1) * 1e+06))

  y.tr <- y.tr[rownames(y.tr) %in% rownames(train_matrix), ]  # Extract rows from data
  y.ts <- y.ts[rownames(y.ts) %in% rownames(train_matrix), ]  # Extract rows from data

  #y <- normalizeBetweenArrays(y,method="none")

  train_dge <- y.tr[keep, ]
  test_dge <- y.ts[keep, ]
  tr_y <- as.matrix(train_dge)
  ts_y <- as.matrix(test_dge)


  #### Variance Filtering ####
  ss = apply(tr_y,1,sd)
  ort = apply(ts_y,1,mean)
  cv = ss/ort
  siralama = order(cv, decreasing = TRUE)


  fit.tr <- lmFit(tr_y, design.tr,block=NULL,correlation=NULL,weights=weights)
  fit.ts <- lmFit(ts_y, design.ts,block=NULL,correlation=NULL,weights=NULL)

  if (is.null(fit.tr$Amean))
    fit.tr$Amean <- rowMeans(tr_y, na.rm = TRUE)

  fit.ts$Amean = fit.tr$Amean
  fit.ts$sigma = fit.tr$sigma
  fit.ts$coefficients = fit.tr$coefficients[,1]


  # #	If no replication found, set all weight to 1
  # NWithReps <- sum(fit$df.residual > 0L)
  # if(NWithReps < 2L) {
  #   if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
  #   if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
  #   out$E <- y
  #   out$weights <- y
  #   out$weights[] <- 1
  #   out$design <- design
  #   if(is.null(out$targets))
  #     out$targets <- data.frame(lib.size=lib.size)
  #   else
  #     out$targets$lib.size <- lib.size
  #   return(new("EList",out))
  # }
  #

  #	Fit lowess trend to sqrt-standard-deviations by log-count-size
  sx <- fit.tr$Amean + mean(log2(lib.size.tr + 1)) - log2(1e+06)
  sy <- sqrt(fit.tr$sigma)


  # allzero <- rowSums(counts)==0
  # if(any(allzero)) {
  #   sx <- sx[!allzero]
  #   sy <- sy[!allzero]
  # }


  l <- lowess(sx,sy,f=span)
  # if(plot) {
  #   plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
  #   title("voom: Mean-variance trend")
  #   lines(l,col="red")
  # }

  #	Make interpolating rule
  #	Special treatment of zero counts is now removed;
  #	instead zero counts get same variance as smallest gene average.
  #	l$x <- c(0.5^0.25, l$x)
  #	l$x <- c(log2(0.5), l$x)
  #	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
  #	var0 <- max(var0,1e-6)
  #	l$y <- c(var0, l$y)

  f <- approxfun(l, rule = 2)

  #	Find individual quarter-root fitted counts
  fitted.values.tr <- fit.tr$coefficients %*% t(fit.tr$design)
  fitted.values.ts <- fit.ts$coefficients %*% t(fit.ts$design)

  # if(fit$rank < ncol(design)) {
  #   j <- fit$pivot[1:fit$rank]
  #   fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
  # } else {
  #   fitted.values <- fit$coefficients %*% t(fit$design)
  # }

  fitted.cpm.tr <- 2^fitted.values.tr
  fitted.cpm.ts <- 2^fitted.values.ts

  fitted.count.tr <- 1e-06 * t(t(fitted.cpm.tr) * (lib.size.tr + 1))
  fitted.count.ts <- 1e-06 * t(t(fitted.cpm.ts) * (lib.size.ts + 1))

  fitted.logcount.tr <- log2(fitted.count.tr)
  fitted.logcount.ts <- log2(fitted.count.ts)

  #	Apply trend to individual observations
  w.tr <- 1/f(fitted.logcount.tr)^4
  w.ts <- 1/f(fitted.logcount.ts)^4

  dim(w.tr) <- dim(fitted.logcount.tr)
  dim(w.ts) <- dim(fitted.logcount.ts)

  dimnames(w.tr) = dimnames(tr_y)
  dimnames(w.ts) = dimnames(ts_y)

  out$E <- tr_y
  out$TestExp <- ts_y
  out$weights <- w.tr
  out$TestWeights <- w.ts
  out$design <- design.tr
  out$TestDesign <- design.ts
  out$Traintargets <- data.frame(lib.size.tr=lib.size.tr)
  out$Testtargets <- data.frame(lib.size.ts=lib.size.ts)
  new("EList", out)

  # #	Output
  # out$E <- y
  # out$weights <- w
  # out$design <- design
  # if(is.null(out$targets))
  #   out$targets <- data.frame(lib.size=lib.size)
  # else
  #   out$targets$lib.size <- lib.size
  # if(save.plot) {
  #   out$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )")
  #   out$voom.line <- l
  # }
  #
  # new("EList",out)
}

voomWithQualityWeights_lasso <- function(data.train, data.test, design=NULL, lib.size=NULL, normalize.method="none", plot=FALSE, span=0.5,
                                       var.design=NULL, var.group=NULL, method="genebygene", maxiter=50, tol=1e-5, trace=FALSE,
                                       col=NULL, ...)
  #	Combine voom weights with sample-specific weights estimated by arrayWeights() function for RNA-seq data
  #	Matt Ritchie, Cynthia Liu, Gordon Smyth
  #	Created 22 Sept 2014. Last modified 7 June 2019.
{
  #	Setup side-by-side plots showing (1) the voom trend and (2) the array weights
  # if(plot) {
  #   oldpar <- par(mfrow=c(1,2))
  #   on.exit(par(oldpar))
  # }

  #	Voom without array weights
  v <- voom_lasso(data.train,data.test, design=design, lib.size=lib.size, normalize.method=normalize.method, plot=FALSE, span=span, ...)

  #	Estimate array weights on top of voom weights
  aw <- limma::arrayWeights(v, design=design, method=method, maxiter=maxiter, tol=tol, var.design=var.design, var.group=var.group)

  #	Update voom weights now using the array weights, plotting trend if requested
  v <- voom_lasso(data.train,data.test, design=design, weights=aw, lib.size=lib.size, normalize.method=normalize.method, plot=plot, span=span, ...)

  #	Update array weights
  aw <- limma::arrayWeights(v, design=design, method=method, maxiter=maxiter, tol=tol, trace=trace, var.design=var.design, var.group=var.group)

  #	Incorporate the array weights into the voom weights
  v$weights <- t(aw * t(v$weights))
  v$targets$sample.weights <- aw

  #	Plot array weights
  if(plot) {
    barplot(aw, names=1:length(aw), main="Sample-specific weights", ylab="Weight", xlab="Sample", col=col)
    abline(h=1, col=2, lty=2)
  }

  v
}



#### Preprocessing function ####

# data will be in the split form already

preprocess <- function(data, preProcessing = "deseq-vst", filterGene = FALSE, filterVariance = FALSE, varianceAmount = 1000, ...){

  if (is.null(preProcessing)){
    cat("Preprocess method not specified. Using the original data without any preprocessing")
    return(data)
  }

  if (preProcessing == "deseq-voom"){

    ### apply deseq-voom preprocessing steps here

    train_data <- data@train
    test_data <- data@test

    ### Extract columns of time and status from train dataset ###
    cols_t <- as.vector(colnames(train_data))
    cols_train <- cols_t[! cols_t %in% c('time', 'status')]
    cols_train <- as.vector(cols_train)
    train_d <- subset(train_data, select=cols_train)

    ### Extract columns of time and status from test dataset ###
    cols_te <- as.vector(colnames(test_data))
    cols_test <- cols_te[! cols_te %in% c('time', 'status')]
    cols_test <- as.vector(cols_test)
    test_d <- subset(test_data, select=cols_test)

    ### Transpose train and test datasets ###
    train_matrix <- as.matrix(t(train_d))
    test_matrix <- as.matrix(t(test_d))

    #### Normalization and transportation ####
    norm.method <- "deseq"
    group <- c(rep(1,ncol(train_matrix)))
    #design = ~1
    lib.size = NULL
    span = 0.5

    transformation <- voomWithQualityWeights_lasso(train_matrix, test_matrix, design=NULL, lib.size=NULL, normalize.method="deseq", plot=FALSE,
                                          span=0.5, var.design=NULL, var.group=NULL, maxiter=50, tol=1e-5,
                                          trace=FALSE, col=NULL, method = "genebygene")
    sirala <- transformation$siralama
    cv_train <- t(input_tr)[sirala[1:2000],]
    cv_test <- t(input_ts)[sirala[1:2000],]

    #### Creating final version of train and test datasets ####
    # This part needs to be cleaned up
    tvsd_train2 <- as.data.frame(cbind(train_ind, train_data$time, train_data$status, t(cv_train), train_W))
    tvsd_test2 <- as.data.frame(cbind(test_ind, test_data$time, test_data$status, t(cv_test)))
    names(tvsd_train2)[names(tvsd_train2) == "V2"] <- "time"
    names(tvsd_train2)[names(tvsd_train2) == "V3"] <- "status"
    names(tvsd_train2)[names(tvsd_train2) == "train_ind"] <- "sirano"
    names(tvsd_train2)[names(tvsd_train2) == "train_W2"] <- "weight"
    tvsd_train2$time <- as.integer(tvsd_train2$time)
    tvsd_train2$status <- as.integer(tvsd_train2$status)
    tvsd_train2$sirano <- as.integer(tvsd_train2$sirano)
    tvsd_train2[] <- apply(tvsd_train2, 2, function(x) as.double((x)))


    names(tvsd_test2)[names(tvsd_test2) == "V2"] <- "time"
    names(tvsd_test2)[names(tvsd_test2) == "V3"] <- "status"
    names(tvsd_test2)[names(tvsd_test2) == "test_ind"] <- "sirano"
    tvsd_test2$time <- as.integer(tvsd_test2$time)
    tvsd_test2$status <- as.integer(tvsd_test2$status)
    tvsd_test2$sirano <- as.integer(tvsd_test2$sirano)
    tvsd_test2[] <- apply(tvsd_test2, 2, function(x) as.double((x)))


    data@train <- tvsd_train2
    data@test <- tvsd_test2
    return(data)
  }

  else if (method == "deseq-vst"){

    ### apply deseq-vst preprocessing steps here

    # change function names to lasso


    train_data <- data@train
    test_data <- data@test

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

    #### near-zero variance filtering ####

    dds_tr <- DESeqDataSetFromMatrix(countData = train_matrix,
                                     colData = as.data.frame(colnames(train_matrix)),
                                     design = ~ 1)
    keep_tr <- rowSums(counts(dds_tr)) >= 10
    dds_tr <- dds_tr[keep_tr,]
    dds_ts <- DESeqDataSetFromMatrix(countData = test_matrix,
                                     colData = as.data.frame(colnames(test_matrix)),
                                     design = ~ 1)
    dds_ts <- dds_ts[keep_tr,]

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


    names(tvsd_test2)[names(tvsd_test2) == "V2"] <- "time"
    names(tvsd_test2)[names(tvsd_test2) == "V3"] <- "status"
    names(tvsd_test2)[names(tvsd_test2) == "test_ind"] <- "sirano"
    tvsd_test2$time <- as.integer(tvsd_test2$time)
    tvsd_test2$status <- as.integer(tvsd_test2$status)
    tvsd_test2$sirano <- as.integer(tvsd_test2$sirano)
    tvsd_test2[] <- apply(tvsd_test2, 2, function(x) as.double((x)))


    # Merge all the relevant data into a single S4 object
    data@train <- tvsd_train2
    data@test <- tvsd_test2
    return(data)
  }

  else {
    cat("Available preprocessing methods are: \"deseq-vst\" and \"deseq-voom\".")
    stop()
  }
}
