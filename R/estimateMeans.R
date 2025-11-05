#' @importFrom Rcpp sourceCpp
#' @importFrom methods new
#' @useDynLib TargetedShrinkage
NULL

#' Estimate the loss given estimated means and known population means
#' 
#' Estimate a quadratic loss function given estimated means and known population means. The estimated means
#' should come from an estimator such as \code{estimatedMeans}, \code{colMeans} or \code{JorionMeans} and
#' give a single estimate of the loss. Estimate the risk by taking the average of multiple estimates. This
#' function computes quadratic loss, which is symmetrical in the two parameters.
#' 
#' The first parameter can be given as a matrix wherein each row is a vector of mean values. Then the function
#' computes the mean quadratic loss, and estimate of risk, by taking the average quadratic loss over the rows.
#' @param estimatedMeans The estimated mean values
#' @param populationMeans The population mean values
#' @param covarianceMatrix A covariance matrix: default is NULL for an identity matrix
#' @return The quadratic loss
#' @export
quadraticLoss <- function(estimatedMeans,populationMeans,covarianceMatrix=NULL){
  if (is.null(covarianceMatrix)) {
    if (!is.matrix(estimatedMeans)) {
      estimatedMeans <- as.numeric(estimatedMeans)
      populationMeans <- as.numeric(populationMeans)
      if (length(estimatedMeans) != length(populationMeans)) {
        stop('estimatedMeans and populationMeans must at least have the same length.')
      }
      return(sum((estimatedMeans-populationMeans)^2))
    } else {
      # estimatedMeans is a matrix
      mean(apply(estimatedMeans, 1, function(row) quadraticLoss(row,populationMeans)))
    }
  } else {
    ## Better: a covariance matrix was supplied
    if (!is.matrix(estimatedMeans)) {
      # vector supplied
      differences <- estimatedMeans-populationMeans
      p <- length(differences)
      ## The following is identical with the start of targetedShrinkage
      eigen.solution <- eigen(covarianceMatrix)
      eigen.values <- eigen.solution$values
      eigen.indices <- which(eigen.values >= sqrt(.Machine$double.eps))
      eigen.values <- eigen.values[eigen.indices] # ignore very small values
      q <- length(eigen.values)
      one.p <- rep(1,p)
      if (q < p) {
        eigen.vectors <- eigen.solution$vectors[,eigen.indices] 
        generalised.inverse <- eigen.vectors %*% diag(1/eigen.values) %*% t(eigen.vectors)
        # compute the quadratic loss
        return(as.numeric(differences %*% generalised.inverse %*% differences))
      } else {
        # compute the quadratic loss directly
        return(as.numeric(differences %*% solve(covarianceMatrix,differences)))
      }
    } else {
      # estimatedMeans is a matrix
      mean(apply(estimatedMeans, 1, function(row) quadraticLoss(row,populationMeans,covarianceMatrix)))
    }
  }
}

#' Compute shrinkage means using a James-Stein estimator.
#' 
#' Compute a James–Stein estimate of shrinkage means. The shrinkage factor \eqn{a} is described
#' as \eqn{\hat w} in
#' Philippe Jorion, Bayes-Stein Estimation for Portfolio Analysis,
#' Journal of Financial and Quantitative Analysis, 21(3), 1986, equation. (11) and shrinkage is
#' \eqn{(1-a)X+at} where \eqn{t} is the parameter `shrinkageTarget`. Typical shrinkage targets are 0
#' (original target used by Stein), the mean of means (used by jorionShrinkage) and the result of
#' `targeted`. This version of James—Stein uses the covariance matrix to weight the risk estimate.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param shrinkageTarget (default is mean of means)
#' @param useLindleyAdjustment Set to TRUE to use adjustment (\eqn{p-3} rather than \eqn{p-2}) suggested by Lindley (1962)
#' when shrinking towards an estimated mean
#' @param debug Set to `TRUE` to obtain debugging information
#' @return The vector of means after shrinkage
#' @export
JamesSteinShrinkage <- function(means,covariances,shrinkageTarget=NULL,useLindleyAdjustment=FALSE,debug=FALSE){
  n <- length(as.vector(means))
  if (useLindleyAdjustment && n < 4){
    stop('James\u2013Stein shrinkage requires p > 4 with Lindley (1962) adjustment')
  } else if (n < 3) {
    stop('James\u2013Stein shrinkage requires p > 3')
  }
  ## sanity checks
  if(nrow(covariances) != n)
    stop('Covariance matrix has wrong number of rows')
  if(ncol(covariances) != n)
    stop('Covariance matrix has wrong number of columns')
  if (is.null(shrinkageTarget)) {
    shrinkageTarget <- mean(means)
  }
  if(debug)
    cat('shrinkage target = ',shrinkageTarget,'\n',sep='')
  r <- means-shrinkageTarget
  b <- solve(covariances,r)
  if (useLindleyAdjustment) {
    alpha <- min(1,(n-3)/as.numeric(r %*% b))
  } else {
    alpha <- min(1,(n-2)/as.numeric(r %*% b))
  }
  if(debug){
    #cat('raw shrinkage factor = ',(n-2)/as.numeric(r %*% b),'\n',sep='')
    cat('shrinkage factor = ',alpha,'\n',sep='')
  }
  result <- alpha*shrinkageTarget+(1-alpha)*means
  return(result)
}

#' Compute shrinkage means using a Jorion estimator.
#' 
#' Compute a Jorion estimate of shrinkage means. The method is taken from Philippe Jorion, Bayes-Stein
#' Estimation for Portfolio Analysis, Journal of Financial and Quantitative Analysis, 21(3), 1986,
#' equations (10) and (17). 
#' This centres on the mean of means and
#' uses a covariance matrix to adjust for covariance and unequal variance.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param T The length of the data generating means
#' @return The vector of means after shrinkage
#' @export
jorionMeans <- function(means,covariances,T){
  n <- length(as.vector(means))
  ## sanity checks
  if(nrow(covariances) != n)
    stop('Covariance matrix has wrong number of rows')
  if(ncol(covariances) != n)
    stop('Covariance matrix has wrong number of columns')
  m <- mean(means)
  r <- means-m
  b <- solve(covariances,r)
  phi <- (n+2)/as.numeric(r %*% b)
  alpha <- phi/(phi+T)
  result <- alpha*m+(1-alpha)*means
  return(result)
}

#' Synonym for jorionMeans.
#' 
#' Compute a Jorion estimate of shrinkage means. This centres on the mean of means and
#' uses a covariance matrix to adjust for covariance and unequal variance.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param T The length of the data generating means
#' @return The vector of means after shrinkage
#' @export
JorionShrinkage <- jorionMeans

#' Estimate the optimal target for weighted mean shrinkage.
#' 
#' Estimate the target that would be used by targetedShrinkage. This function is provided mainly
#' to allow testing of alternatives to the James–Stein estimator used by targetedShrinkage.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param debug Set to `TRUE` to obtain debugging information 
#' @return The vector of means after shrinkage
#' @export
target <- function(means,covariances,debug=FALSE){
  p <- length(as.vector(means))
  ## sanity checks
  if(nrow(covariances) != p)
    stop('Covariance matrix has wrong number of rows')
  if(ncol(covariances) != p)
    stop('Covariance matrix has wrong number of columns')
  ## The following is identical with the start of targetedShrinkage
  eigen.solution <- eigen(covariances)
  eigen.values <- eigen.solution$values
  eigen.indices <- which(eigen.values >= sqrt(.Machine$double.eps))
  eigen.values <- eigen.values[eigen.indices] # ignore very small values
  q <- length(eigen.values)
  one.p <- rep(1,p)
  if (q < p) {
    eigen.vectors <- eigen.solution$vectors[,eigen.indices] 
    generalised.inverse <- eigen.vectors %*% diag(1/eigen.values) %*% t(eigen.vectors)
    # compute the shrinkage target
    precision.one.p <- generalised.inverse %*% one.p
  } else {
    precision.one.p <- solve(covariances,one.p)
  }
  # compute the shrinkage target
  target.numerator <- means %*% precision.one.p
  shrinkageTarget <- as.numeric(target.numerator/(one.p %*% precision.one.p))
  if(debug)
    cat('Target = ',shrinkageTarget,'.\n',sep='')
  ## The preceding was identical with the start of targetedShrinkage
  return(shrinkageTarget)
}

#' Create a shrinkage estimate based on a covariance matrix and mean values.
#' 
#' Unlike Jorion shrinkage, this method does not shrink towards a mean of mean values (or indeed towards 0) but
#' to a good estimate of the optimal shrinkage target, which is a mean weighted by covariance matrix. Otherwise it
#' behaves like a conventional James–Stein estimator.
#' The parameter `covariances` will typically be a Ledoit–Wolf shrinkage estimate
#' of the covariance matrix of the data (which we use to estimate the covariance matrix of the means), but could, for
#' example be another covariance matrix estimate or the known covariance matrix.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param debug Set to `TRUE` to obtain debugging information
#' @param scale Set to a positive value to not use James–Stein scaling.
#' @return The vector of means after shrinkage
#' @export
targetedShrinkage <- function(means,covariances,debug=FALSE,scale=NA){
  p <- length(as.vector(means))
  ## sanity checks
  if(nrow(covariances) != p)
    stop('Covariance matrix has wrong number of rows')
  if(ncol(covariances) != p)
    stop('Covariance matrix has wrong number of columns')
  if (p < 4 && is.na(scale)) {
    stop('targeted shrinkage requires p > 3 unless a scale is given')
  }
  one.p <- rep(1,p)
  eigen.solution <- eigen(covariances)
  eigen.values <- eigen.solution$values
  eigen.indices <- which(eigen.values >= sqrt(.Machine$double.eps))
  eigen.values <- eigen.values[eigen.indices] # ignore very small values
  q <- length(eigen.values)
  one.p <- rep(1,p)
  if (q < p) {
    if (q < 4 && is.na(scale)) {
      stop(paste(sep='','targeted shrinkage requires q > 3 (q = ',q,', p = ',p,
                 ') unless a scale is given'))
    }
    eigen.vectors <- eigen.solution$vectors[,eigen.indices] 
    generalised.inverse <- eigen.vectors %*% diag(1/eigen.values) %*% t(eigen.vectors)
    # compute the shrinkage target
    precision.one.p <- generalised.inverse %*% one.p
  } else {
    precision.one.p <- solve(covariances,one.p)
  }
  # compute the shrinkage target
  target.numerator <- means %*% precision.one.p
  shrinkageTarget <- as.numeric(target.numerator/(one.p %*% precision.one.p))
  if(debug)
    cat('Shrinkage target = ',shrinkageTarget,'.\n',sep='')
  # compute the shrinkage factor
  if (q < p) {
    precision.x <- generalised.inverse %*% means
  } else {
    precision.x <- solve(covariances,means)
  }
  if (debug) {
    print(solve(covariances)) 
  }
  if (is.na(scale)) {
    numerator <- (q-3)
  } else {
    numerator <- p*scale
  }
  denominator <- means %*% precision.x - shrinkageTarget * target.numerator
  alpha <- as.numeric(numerator/denominator)
  alpha <- min(1,alpha)
  ## alpha <- max(0,alpha) # not needed.
  if(debug){
    cat('shrinkage factor',alpha,'\n')
    cat('target',shrinkageTarget,'\n')
  }
  result <- (1-alpha)*means+alpha*shrinkageTarget;
  return(result)
}

#' Estimate mean shrinkage under a whitening transformation.
#' 
#' This function shrinks under
#' a whitening transformation. The whitening transformation uses principal component analysis for whitening,
#' though other whitening transformations could easily be added. The effect is a nonuniform shrinkage or,
#' equivalently, uniform shrinkage but with different targets for different components of the vector of
#' means.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param useLindleyAdjustment Set to TRUE to use adjustment (\eqn{p-3} rather than \eqn{p-2}) suggested by Lindley (1962)
#' when shrinking towards an estimated mean
#' @param debug Set to `TRUE` to obtain debugging information 
#' @return The vector of means after shrinkage
#' @export
whiteningShrinkage <- function(means,covariances,useLindleyAdjustment=FALSE,debug=FALSE){
  # The following is close to the start of targetedShrinkage
  p <- length(as.vector(means))
  if (useLindleyAdjustment && p < 4){
    stop('Whitening shrinkage requires p > 4 with Lindley (1962) adjustment')
  } else if (p < 3) {
    stop('Whitening shrinkage requires p > 3')
  }
  ## sanity checks
  if(nrow(covariances) != p)
    stop('Covariance matrix has wrong number of rows')
  if(ncol(covariances) != p)
    stop('Covariance matrix has wrong number of columns')
  Sigma.inv <- LaplacesDemon::Cov2Prec(covariances)
  ed <- eigen(Sigma.inv,symmetric = TRUE)
  D <- diag(sqrt(ed$values))
  B <- ed$vectors %*% D
  # The lines above were identical with the start of targetedShrinkage
  ## Carry out whitening transformation
  Y <- B %*% means
  if (debug) {
    cat('Y: ',Y,'\n')
  }
  ## Find target (in Y)
  Y.target <- mean(Y)
  ## Find shrinkage factor
  denominator <- 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      denominator <- denominator + (Y[j]-Y[i])^2  
    }
  }
  if (debug) {
    cat('numerator:',p(p-2),'\n')
    cat('denominator:',denominator,'\n')
  }
  if (useLindleyAdjustment) {
    alpha <- min(1,p*(p-3)/denominator)
  } else {
    alpha <- min(1,p*(p-2)/denominator)
  }
  if (debug) {
    cat('alpha:',alpha,'\n')
    cat('target:',Y.target,'\n')
  }
  ## Carry out shrinkage in transformed data
  Y.hat <- (1-alpha)*Y+alpha*Y.target
  if (debug) {
    cat('Y.hat:',Y.hat,'\n')
  }
  ## Colouring transformation
  B.inv <- diag(1/diag(D)) %*% t(ed$vectors)
  return(as.numeric(B.inv %*% Y.hat))
}

#' Create a nonuniform shrinkage estimate based on a covariance matrix and mean values.
#' 
#' Unlike Jorion shrinkage, this method does not shrink towards a mean of mean values (or indeed towards 0) but
#' to a good estimate of the optimal shrinkage target, which is a mean weighted by covariance matrix. It then uses
#' shrinkage between consecutive pairs of means to try to get a non-uniform shrinkage. The shrinkage should match
#' that of `targetedShrinkage` in the special case where the covariance matrix is a muliple of the identity matrix
#' and the means are uniformly spaced.
#' The parameter `covariances` will typically be a Ledoit–Wolf shrinkage estimate
#' of the covariance matrix of the data (which we use to estimate the covariance matrix of the means), but could, for
#' example be another covariance matrix estimate or the known covariance matrix.
#' @param means The vector of means
#' @param covariances The matrix of covariances
#' @param T The length of the data generating means
#' @param debug Set to `TRUE` to obtain debugging information 
#' @return The vector of means after shrinkage
#' @export
nonuniformShrinkage <- function(means,covariances,T,debug=FALSE){
  ## The following is identical with the start of targetedShrinkage
  p <- length(as.vector(means))
  if (p < 3) {
    stop('Function only defined for p > 2.')
  }
  ## sanity checks
  if(nrow(covariances) != p)
    stop('Covariance matrix has wrong number of rows')
  if(ncol(covariances) != p)
    stop('Covariance matrix has wrong number of columns')
  Sigma.inv <- LaplacesDemon::Cov2Prec(covariances)
  ed <- eigen(Sigma.inv,symmetric = TRUE)
  D <- diag(sqrt(ed$values))
  B <- ed$vectors %*% D
  ## The preceding was identical with the start of targetedShrinkage
  ## Carry out whitening transformation
  Y <- B %*% means
  if (debug) {
    cat('Y: ',Y,'\n')
  }
  ## Find target (in Y)
  Y.target <- mean(Y)
  y.hat <- rep(0,p)
  ## Keep track of order of Y
  o <- order(Y)
  ## Shrinkage towards Y.target
  c <- 12*(p-2)/(p^2*(p-1))
  for (i in 2:p) {
    alpha <- c /(Y[o[i]]-Y[o[i-1]])^2
    y.hat[o[i]] <- y.hat[o[i-1]]+(1-alpha)*(Y[o[i]]-Y[o[i-1]])
  }
  y.hat.mean <- mean(y.hat)
  y.hat <- y.hat-y.hat.mean+Y.target
  ## Colouring transformation
  B.inv <- diag(1/diag(D)) %*% t(ed$vectors)
  return(as.numeric(B.inv %*% y.hat))
}

#' Create a bootstrap matrix suitable for `estimateMeansObj`.
#' 
#' Create a matrix suitable for `estimateMeansObj` with lots of rows and one column for each row of the data
#' it is used to resample.
#' @param T The number of rows in the data
#' @param B The number of bootstrap replications (default = 10000)
#' @return A `matrix` suitable for estimateMeansObj
#' @export
createBootstrapMatrix <- function(T,B=10000){
  return(matrix(nrow=B,sample(x = 1:T,size = B*T,replace = TRUE)))
}

#' Create a function object based on a bootstrap matrix. This matrix of positive integers 
#'
#' Create and return an object that can be used for the bootstrap shrinkage estimator. The matrix `bootstrapMatrix`
#' should have one column for each row of any data it is applied to and many (2000–10000 typically) rows. Each row
#' should have integers in the range 1,…,T, where T is the number of rows of data to which the bootstrap will be applied.
#' 
#' @param bootstrapMatrix The matrix used for bootstrap
#' @return the object: an object of class EstimateMeans
#' @seealso createBootstrapMatrix
#' @seealso estimateMeans
#' @export
estimateMeansObj <- function(bootstrapMatrix){
  # Check matrix makes some sense
  if (!is.matrix(bootstrapMatrix)) {
    stop('bootstrapMatrix must be a matrix.')
  }
  if (!is.integer(bootstrapMatrix)) {
    stop('bootstrapMatrix must have only integer entries.')
  }
  if (min(bootstrapMatrix) < 1) {
    stop('bootstrapMatrix must have only have positive integer entries.')
  }
  T <- ncol(bootstrapMatrix)
  if (max(bootstrapMatrix) > T) {
    stop(paste('bootstrapMatrix entries cannot exceed ',T,' in value.',sep=''))
  }
  env <- Rcpp::Module("OrderBasedShrinkageModule",PACKAGE = "OrderBasedShrinkage")
  return(new(env$EstimateMeans,bootstrapMatrix))
}
