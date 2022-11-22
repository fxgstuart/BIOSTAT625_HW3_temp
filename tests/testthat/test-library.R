fit <- lm(hp~mpg+disp+drat+wt+qsec, data = mtcars)
tb <- summary(fit)
X <- as.matrix(mtcars[,c("mpg","disp","drat","wt","qsec")])
y <- as.vector(mtcars[,"hp"])
Xtrain <- X[1:24,]
ytrain <- y[1:24]
Xtest <- X[25:32,]
fit2 <- lm(ytrain~Xtrain)

test_that("fit.coef works", {
  expect_equal(fit.coef(X,y), unname(fit$coefficients))
})

test_that("fit.values works",{
  expect_equal(fit.values(X,y),unname(fit$fitted.values))
})

test_that("fit.predict works",{
  expect_equal(fit.predict(Xtrain,ytrain,Xtest),unname(predict.lm(fit2,Xtest)))
})

test_that("fit.sigma.squared works",{
  expect_equal(sqrt(fit.sigma.squared(X,y)),unname(tb$sigma))
})

test_that("fit.sd works",{
  expect_equal(unname(fit.sd(X,y)),unname(tb$coefficients[,2]))
})

test_that("fit.t.test works",{
  expect_equal(unname(cbind(fit.t.test(X,y)$t,fit.t.test(X,y)$p.value)),unname(tb$coefficients[,c(3,4)]))
})

test_that("fit.confint works",{
  expect_equal(unname(fit.confint(X,y)),unname(confint(fit)))
})

test_that("fit.hat.matrix works",{
  expect_equal(unname(diag(fit.hat.matrix(X))),unname(influence(fit)$hat))
})

test_that("fit.overall.test works",{
  res <- fit.overall.test(X,y)
  res_vec <- c(res$F,res$df1,res$df2)
  expect_equal(unname(res_vec),unname(tb$fstatistic))
})

test_that("fit.R.squared works",{
  res <- fit.R.squared(X,y)
  res_vec <- c(res$R.squared,res$adj.R.squared)
  expect_equal(unname(res_vec),unname(c(tb$r.squared,tb$adj.r.squared)))
})

test_that("fit.partial.test works",{
  res <- fit.partial.test(X,y,idx = c(2))
  res_vec <- c(res$F.stat,res$p.value)
  true_table <- car::Anova(fit,type="III")
  expect_equal(unname(res_vec),unname(c(true_table$`F value`[3],true_table$`Pr(>F)`[3])))
})

test_that("fit.GLH.test works",{
  T.matrix <- matrix(0,nrow=2,ncol=6)
  T.matrix[1,3]<-1
  T.matrix[1,4]<--1
  T.matrix[2,4]<-1
  T.matrix[2,5]<--1
  res <- fit.GLH.test(X,y,T.matrix = T.matrix)
  res_vec <- c(res$F.stat,res$p.value)
  true_table <- car::linearHypothesis(model=fit,hypothesis.matrix=T.matrix,rhs=c(0,0))
  expect_equal(unname(res_vec),unname(c(true_table$F[2],true_table$`Pr(>F)`[2])))
})
