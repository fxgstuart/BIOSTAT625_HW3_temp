fit <- lm(hp~mpg+disp+drat+wt+qsec, data = mtcars)
X <- as.matrix(mtcars[,c("mpg","disp","drat","wt","qsec")])
y <- as.vector(mtcars[,"hp"])

test_that("fit.coef works", {
  expect_equal(fit.coef(X,y), unname(fit$coefficients))
})


