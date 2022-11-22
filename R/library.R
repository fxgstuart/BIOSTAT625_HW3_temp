#'fit.coef
#'
#'Gets the estimated coefficients of a linear model
#'
#'@import bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return the estimated coefficients of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.coef(X,y)
#'
#'@export
#'
fit.coef <- function(X,y,intercept=T){
  if(intercept){
    X <- cbind(rep(1,nrow(X)),X)
  }
  return(as.vector(solve(crossprod(X))%*%crossprod(X,y)))
}


#'fit.values
#'
#'Gets the fitted values of a linear model
#'
#'@import bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return the fitted values of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.values(X,y)
#'
#'@export
#'
fit.values <- function(X,y,intercept=T){
  beta <- fit.coef(X,y,intercept)
  if(intercept){
    X <- cbind(rep(1,nrow(X)),X)
  }
  return(as.vector(X%*%beta))
}

#'fit.sigma.squared
#'
#'Gets the estimated residual of a linear model
#'
#'@import bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return the estimated residual of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.sigma.squared(X,y)
#'
#'@export
#'
fit.sigma.squared <- function(X,y,intercept=T){
  n <- nrow(X)
  p <- ncol(X)

  p <- p+intercept
  y.hat <- fit.values(X,y,intercept)
  epsilon.hat <- y-y.hat
  mse <- as.numeric(crossprod(epsilon.hat))/(n-p)

  return(as.numeric(mse))
}

#'fit.sd
#'
#'Gets the estimated standard deviation of each parameter(beta) of a linear model
#'
#'@import bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return The estimated standard deviation of each parameter of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.sd(X,y)
#'
#'@export
#'
fit.sd <- function(X,y,intercept=T){
  sigma.squared <- fit.sigma.squared(X,y,intercept)

  if(intercept){
    X <- cbind(rep(1,nrow(X)),X)
  }
  var.betahat <- diag(solve(crossprod(X)))*c(sigma.squared)
  sd.betahat <- sqrt(var.betahat)

  return(sd.betahat)
}

#'fit.t.test
#'
#'Conducting t test for each parameter(beta) of a linear model
#'
#'@import bench
#'
#'@importFrom stats pf pt qt
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return The test statistic and corresponding p-value of each parameter of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.t.test(X,y)
#'
#'@export
#'
fit.t.test <- function(X,y,intercept=T){
  n <- nrow(X)
  p <- ncol(X)+intercept
  betahat <- fit.coef(X,y,intercept)
  sd.betahat <- fit.sd(X,y,intercept)
  t.stat <- c(betahat/sd.betahat)
  p.value <- c(2*pmin(1-pt(t.stat,n-p),pt(t.stat,n-p)))

  return(list(t = t.stat, p.value = p.value))
}

#'fit.confint
#'
#'Construct the confidence intervals parameters of a linear model
#'
#'@import bench
#'
#'@importFrom stats pf pt qt
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@param level numerical value, stands for the significance of the confidence interval
#'
#'@return The upper and lower bounds of confidence interval for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.confint(X,y)
#'
#'@export
#'
fit.confint <- function(X,y,intercept=T,level=0.95){
  n <- nrow(X)
  p <- ncol(X)+intercept
  betahat <- fit.coef(X,y,intercept)
  sd.betahat <- fit.sd(X,y,intercept)
  t.value <- qt(1-(1-level)/2,n-p)

  beta.CI.upper <- betahat + t.value * sd.betahat
  beta.CI.lower <- betahat - t.value * sd.betahat
  CI.beta.matrix <- cbind(beta.CI.lower, beta.CI.upper)
  colnames(CI.beta.matrix) <- c("lower","upper")

  return(CI.beta.matrix)
}

#'fit.hat.matrix
#'
#'Calculate the hat matrix of a linear model
#'
#'@import bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return The upper and lower bounds of confidence interval for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.hat.matrix(X)
#'
#'@export
#'
fit.hat.matrix <- function(X,intercept=T){
  if(intercept){
    X <- cbind(rep(1,nrow(X)),X)
  }

  return(X%*%solve(crossprod(X))%*%t(X))
}

#'fit.overall.test
#'
#'Conduct the test of significance (overall test) of a linear model
#'
#'@import bench
#'
#'@importFrom stats pf pt qt
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return The test statistic and corresponding p-value of the overall test for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.overall.test(X,y)
#'
#'@export
#'
fit.overall.test <- function(X,y,intercept=T){
  n <- nrow(X)
  p <- ncol(X)+intercept
  yhat <- fit.values(X,y,intercept)
  SSE <- sum((y-yhat)^2)
  SSR <- sum((yhat-mean(y))^2)

  F.stat <- (SSR/(p-intercept))/(SSE/(n-p))
  p.value <- 1-pf(F.stat, p-intercept, n-p)
  return(list(F = F.stat, p.value = p.value, df1 = p-intercept, df2 = n-p))
}

#'fit.R.squared
#'
#'Gets the estimated R squared of a linear model
#'
#'@import bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@return The estimated R squared of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.R.squared(X,y)
#'
#'@export
#'
fit.R.squared <- function(X,y,intercept=T){
  n <- nrow(X)
  p <- ncol(X)+intercept
  yhat <- fit.values(X,y,intercept)
  SSE <- sum((y-yhat)^2)
  SSR <- sum((yhat-mean(y))^2)

  res1 <- SSR/(SSR+SSE)
  res2 <- 1-(SSE/(n-p))/((SSE+SSR)/(n-intercept))

  return(list(R.squared = res1, adj.R.squared = res2))
}

#'fit.partial.test
#'
#'Conduct the partial test for several given parameters of a linear model
#'
#'@import bench
#'
#'@importFrom stats pf pt qt
#'
#'@import car bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@param idx a vector of indices of tested parameters
#'
#'@return The test statistic and corresponding p-value of the partial test for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'fit.partial.test(X,y,idx=c(3:5))
#'
#'@export
#'
fit.partial.test <- function(X,y,intercept=T,idx){
  X.temp <- X[,-c(idx)]
  yhat.temp <- fit.values(X.temp,y,intercept)
  SSR.temp <- sum((yhat.temp-mean(y))^2)

  yhat <- fit.values(X,y,intercept)
  SSR <- sum((yhat-mean(y))^2)

  SS <- SSR-SSR.temp
  sigma.squared <- fit.sigma.squared(X,y,intercept)

  F.stat <- (SS/length(idx))/sigma.squared
  p.value <- 1-pf(F.stat, length(idx), nrow(X)-ncol(X)-intercept)

  return(list(F.stat = F.stat, p.value = p.value))
}


#'fit.GLH.test
#'
#'Conduct the general linear hypothesis test for several given parameters of a linear model
#'
#'@importFrom stats pf pt qt
#'
#'@import car bench
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param intercept logical, "TRUE" by default.
#'If "FALSE", the model will not include the intercept.
#'
#'@param T.matrix hypothesis matrix (or vector) giving linear combinations of coefficients by rows
#'
#'@param c right-hand-side vector for hypothesis, with as many entries as rows in the hypothesis matrix;
#'can be omitted, in which case it defaults to a vector of zeroes
#'
#'@return The test statistic and corresponding p-value of the general linear hypothesis test for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'y <- rnorm(12)
#'T.matrix <- matrix(0,nrow=2,ncol=11)
#'T.matrix[1,4]<-1
#'T.matrix[1,5]<--1
#'T.matrix[2,5]<-1
#'T.matrix[2,6]<--1
#'fit.GLH.test(X,y,T.matrix=T.matrix)
#'
#'@export
#'
fit.GLH.test <- function(X,y,intercept=T,T.matrix,c=NULL){
  if(is.null(c)) c<-rep(0,nrow(T.matrix))
  betahat <- fit.coef(X,y,intercept)
  sigma.squared <- fit.sigma.squared(X,y,intercept)
  T.rank <- qr(T.matrix)$rank
  if(intercept){
    X <- cbind(rep(1,nrow(X)),X)
  }

  vec.temp <- T.matrix%*%betahat-c
  mat.temp <- T.matrix%*%solve(crossprod(X))%*%t(T.matrix)

  F.stat <- t(vec.temp)%*%solve(mat.temp)%*%vec.temp/(T.rank)/sigma.squared
  p.value <- 1-pf(F.stat, T.rank, nrow(X)-ncol(X))

  return(list(F.stat = F.stat, p.value = p.value))
}
