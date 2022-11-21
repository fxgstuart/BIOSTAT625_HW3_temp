fit <- lm(hp~mpg+disp+drat+wt+qsec, data = mtcars)
X <- as.matrix(mtcars[,c("mpg","disp","drat","wt","qsec")])
y <- as.vector(mtcars[,"hp"])
Xtrain <- X[1:24,]
ytrain <- y[1:24]
Xtest <- X[25:32,]
fit2 <- lm(ytrain~Xtrain)

test_that("fit.coef works", {
  expect_equal(fit.coef(X,y), unname(fit$coefficients))
})

test_that("fit.values works"){
  expect_equal(fit.values(X,y),unname(fit$fitted.values))
}

test_that("fit.predict works"){
  expect_equal(fit.predict(Xtrain,ytrain,Xtest),unname(predict.lm(fit2,Xtest)))
}


