
fit.coef <- function(X,y){
  return(as.vector(solve(crossprod(X))%*%crossprod(X,y)))
}

fit.values <- function(X,y){
  beta <- fit.coef(X,y)
  return(as.vector(X%*%beta))
}

fit.predict <- function(X,y,Xnew){
  beta <- fit.coef(X,y)
  return(as.vector(Xnew%*%beta))
}

fit.sigma.squared <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  y.hat <- fit.values(X,y)
  epsilon.hat <- y-y.hat
  mse <- crossprod(epsilon.hat)/(n-p)
  
  return(mse)
}

fit.sd <- function(X,y){
  sigma.squared <- fit.sigma.squared(X,y)
  var.betahat <- diag(crossprod(X))*c(sigma.squared)
  sd.betahat <- sqrt(var.betahat)
  
  return(sd.betahat)
}

fit.t.test <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  betahat <- fit.coef(X,y)
  sd.betahat <- fit.sd(X,y)
  t.stat <- c(betahat/sd.betahat)
  p.value <- c(2*(1-pt(t.stat,n-p)))
  
  return(list(t = t.stat, p.value = p.value))
}

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

fit.R.squared <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  yhat <- fit.values(X,y)
  SSE <- sum((y-yhat)^2)
  SSR <- sum((yhat-mean(y))^2)
  
  return(SSR/(SSR+SSE))
}

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

fit.GLH.test <- function(X,y,T.matrix,c){
  betahat <- fit.coef(X,y)
  sigma.squared <- fit.sigma.squared(X,y)
  T.rank <- qr(T.matrix)$rank
  
  vec.temp <- T.matrix%*%betahat-c
  mat.temp <- T.matrix%*%solve(crossprod(X))%*%t(T.matrix)
  
  F.stat <- t(vec.temp)%*%solve(mat.temp)%*%vec.temp/(T.rank)/sigma.squared
  p.value <- 1-pf(F.stat, T.rank, nrow(X)-ncol(X))
  
  return(list(F.stat = F.stat, p.value = p.value))
}
