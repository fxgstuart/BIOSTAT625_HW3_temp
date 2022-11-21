#'fit.coef
#'
#'Gets the estimated coefficients of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return the estimated coefficients of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
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
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return the fitted values of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.values(X,y)
#'
#'@export
#'
fit.values <- function(X,y){
  beta <- fit.coef(X,y)
  return(as.vector(X%*%beta))
}

#'fit.predict
#'
#'Gets the predict values of a linear model given some new data
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param Xnew covariate matrix of new data that needs to predict
#'
#'@return the predicted values of Xnew by linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'Xnew <- matrix(rnorm(100),nrow=10,ncol=10)
#'Xnew <- cbind(rep(1,10),Xnew)
#'y <- rnorm(12)
#'fit.predict(X,y,Xnew)
#'
#'@export
#'
fit.predict <- function(X,y,Xnew){
  beta <- fit.coef(X,y)
  return(as.vector(Xnew%*%beta))
}

#'fit.sigma.squared
#'
#'Gets the estimated residual of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return the estimated residual of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.sigma.squared(X,y)
#'
#'@export
#'
fit.sigma.squared <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  y.hat <- fit.values(X,y)
  epsilon.hat <- y-y.hat
  mse <- crossprod(epsilon.hat)/(n-p)

  return(mse)
}

#'fit.sd
#'
#'Gets the estimated standard deviation of each parameter(beta) of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return The estimated standard deviation of each parameter of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.sd(X,y)
#'
#'@export
#'
fit.sd <- function(X,y){
  sigma.squared <- fit.sigma.squared(X,y)
  var.betahat <- diag(crossprod(X))*c(sigma.squared)
  sd.betahat <- sqrt(var.betahat)

  return(sd.betahat)
}

#'fit.t.test
#'
#'Conducting t test for each parameter(beta) of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return The test statistic and corresponding p-value of each parameter of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.t.test(X,y)
#'
#'@export
#'
fit.t.test <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  betahat <- fit.coef(X,y)
  sd.betahat <- fit.sd(X,y)
  t.stat <- c(betahat/sd.betahat)
  p.value <- c(2*(1-pt(t.stat,n-p)))

  return(list(t = t.stat, p.value = p.value))
}

#'fit.overall.test
#'
#'Conduct the test of significance (overall test) of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return The test statistic and corresponding p-value of the overall test for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.overall.test(X,y)
#'
#'@export
#'
fit.overall.test <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  yhat <- fit.values(X,y)
  SSE <- sum((y-yhat)^2)
  SSR <- sum((yhat-mean(y))^2)

  F.stat <- (SSR/p)/(SSE/(n-p))
  p.value <- 1-pf(F.stat, p, n-p)
  return(list(F = F.stat, p.value = p.value))
}

#'fit.R.squared
#'
#'Gets the estimated R squared of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@return The estimated R squared of linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.R.squared(X,y)
#'
#'@export
#'
fit.R.squared <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  yhat <- fit.values(X,y)
  SSE <- sum((y-yhat)^2)
  SSR <- sum((yhat-mean(y))^2)

  return(SSR/(SSR+SSE))
}

#'fit.partial.test
#'
#'Conduct the partial test for several given parameters of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
#'
#'@param idx a vector of indices of tested parameters(cannot includes 1, which stands for intercept)
#'
#'@return The test statistic and corresponding p-value of the partial test for linear model constructing by X and y
#'
#'@examples
#'X <- matrix(rnorm(120),nrow=12,ncol=10)
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.overall.test(X,y)
#'
#'@export
#'
fit.partial.test <- function(X,y,idx){
  X.temp <- X[,-c(idx)]
  yhat.temp <- fit.values(X.temp,y)
  SSR.temp <- sum((yhat.temp-mean(y))^2)

  yhat <- fit.values(X,y)
  SSR <- sum((yhat-mean(y))^2)

  SS <- SSR-SSR.temp
  sigma.squared <- fit.sigma.squared(X,y)

  F.stat <- (SS/length(idx))/sigma.squared
  p.value <- 1-pf(F.stat, length(idx), nrow(X)-ncol(X))

  return(list(F.stat = F.stat, p.value = p.value))
}


#'fit.GLH.test
#'
#'Conduct the general linear hypothesis test for several given parameters of a linear model
#'
#'@param X input value of covariate matrix/vector
#'
#'@param y input value of responses
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
#'X <- cbind(rep(1,12),X)
#'y <- rnorm(12)
#'fit.overall.test(X,y)
#'
#'@export
#'
fit.GLH.test <- function(X,y,T.matrix,c=NULL){
  if(is.null(c)) c<-rep(0,nrow(T.matrix))
  betahat <- fit.coef(X,y)
  sigma.squared <- fit.sigma.squared(X,y)
  T.rank <- qr(T.matrix)$rank

  vec.temp <- T.matrix%*%betahat-c
  mat.temp <- T.matrix%*%solve(crossprod(X))%*%t(T.matrix)

  F.stat <- t(vec.temp)%*%solve(mat.temp)%*%vec.temp/(T.rank)/sigma.squared
  p.value <- 1-pf(F.stat, T.rank, nrow(X)-ncol(X))

  return(list(F.stat = F.stat, p.value = p.value))
}
