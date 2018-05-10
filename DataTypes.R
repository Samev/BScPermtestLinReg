#Different data generation functions
#For simple data the first column returned is the explanatory X and the second is the response Y.
#For multivariate data the second last column is the variable of interest and the last is the response Y. The previous variables are nuicance variables.

#Heavy tails
simpleHeavyTail <- function(n,beta1,v){
  x <- rexp(n,1)
  return(cbind(x, beta1*x + rt(n,v)/sqrt(v/(v-2))))
}

#Skew errors with exponential X distribution. Errors are standardized.
simpleExpSkew_XExp <- function(n, beta1, m){
  x <- rexp(n,1)
  error <- rexp(n,1)^m
  error <- (error-mean(error))/sd(error)
  return(cbind(x, beta1*x+error))
}

#Skew errors with uniform X distributions. Errors are standardized
simpleExpSkew_XUnif <- function(n, beta1, m){
  x <- runif(n,min=0, max=3)
  error <- rexp(n,1)^m
  error <- (error-mean(error))/sd(error)
  return(cbind(x, beta1*x+error))
}

#Outliers
simpleOutliers <- function(n,beta1, outliers){
  x <- rexp(n,1)
  error <- rnorm(n,0,1)
  out <- runif(n,0,1)>(1-outliers)
  error[out] <- sign(error[out])*10
  return(cbind(x,beta1*x + error))
}

#Non-constant variance.
simpleNonConstVar <- function(n, beta1, v){
  x <- rexp(n,1)
  error <- rnorm(n,0,1+v*abs(x)/6)
  return(cbind(x,beta1*x + error))
}


#Multivariate heavy tails
multHeavyTail <- function(n, beta, rho, sig, v){
  sigMatrix <- diag(1, length(beta))
  sigMatrix[sigMatrix == 0] <- rho
  x <- mvrnorm(n = n, mu = rep(0, times=length(beta)), Sigma = sig*sigMatrix)
  error <- rt(n,v)/sqrt(v/(v-2))
  return(cbind(x, x %*% beta + error))
}

#Multivariate skew errors
multExpSkew <- function(n, beta, rho, sig, m){
  sigMatrix <- diag(1, length(beta))
  sigMatrix[sigMatrix == 0] <- rho
  x <- mvrnorm(n = n, mu = rep(0, times=length(beta)), Sigma = sig*sigMatrix)
  error <- rexp(n,1)^m
  error <- (error-mean(error))/sd(error)
  return(cbind(x, x %*% beta + error))
}

#Multivariate non-constant variance
multNonConstVar <- function(n, beta, rho, sig, v){
  sigMatrix <- diag(1, length(beta))
  sigMatrix[sigMatrix == 0] <- rho
  x <- mvrnorm(n = n, mu = rep(0, times=length(beta)), Sigma = sig*sigMatrix)
  error <- rnorm(n,0,1+v*abs(x[, length(beta)])/6)
  return(cbind(x, x %*% beta + error))
}

#Multivariate outliers
multOutliers <- function(n, beta, rho, sig, v){
  sigMatrix <- diag(1, length(beta))
  sigMatrix[sigMatrix == 0] <- rho
  x <- mvrnorm(n = n, mu = rep(0, times=length(beta)), Sigma = sig*sigMatrix)
  error <- rnorm(n,0,1)
  out <- runif(n,0,1)>(1-v)
  error[out] <- sign(error[out])*10
  return(cbind(x, x %*% beta + error))
}