# FitLM

<!-- badges: start -->
[![R-CMD-check](https://github.com/fxgstuart/BIOSTAT625_HW3_temp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fxgstuart/BIOSTAT625_HW3_temp/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/fxgstuart/BIOSTAT625_HW3_temp/branch/main/graph/badge.svg)](https://app.codecov.io/gh/fxgstuart/BIOSTAT625_HW3_temp?branch=main)
<!-- badges: end -->

## Background

### Linear regression

In statistics, [linear regression](https://en.wikipedia.org/wiki/Linear_regression) is a linear approach for modelling the relationship between a scalar response and one or more explanatory variables (also known as dependent and independent variables). In linear regression, the relationships are modeled using linear predictor functions whose unknown model parameters are estimated from the data.

### "lm" function in R

The existing `stats::lm` in `R` is a commonly used function that can help fit linear regression model. It can provide us with estimation, related statistics, and hypothesis test results. 

However, the `stats::lm` function can not directly output some important results,including the confidence interval(CI) of $\hat{\beta}$, general linear hypothesis(GLH) test results, partial F test results and hat matrix.

Therefore, I develop this `FitLM` package to include most basic functions of `stats::lm`, and implement some new functions such as `FitLM::fit.partial.test` and `FitLM::fit.GLH.test` to improve the performance when analyzing linear regression models.

## Basic usage

### installation

* `FitLM` can be downloaded from github:

```
devtools::install_github("fxgstuart/BIOSTAT625_HW3_temp", build_vignettes = T)
```

### Loading packages

```
library(FitLM)
```


Use the following code to find more details about the "FitLM" package:

```
browseVignettes("FitLM")
```
