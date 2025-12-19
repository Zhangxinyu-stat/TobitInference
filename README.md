# TobitInference
R code for "Simultaneous inference for high-dimensional Tobit regression model".

## Usage
- `InferenceTheta(x, y, c_q, alpha, B)`
- `InferenceBeta(x, y, c_q, alpha, B)`

## Required Packages
- `MASS`
- `clime`
- `tobitnet`

## Inputs
- `x`: An `n × p` design matrix of predictors, where `n` is the sample size and `p` is the dimension.
- `y`: A numeric vector of length `n` observed responses under left-censoring.
- `c_q`: The censoring threshold. 
- `alpha`: Significance level for simultaneous inference.
- `B`: Number of multiplier bootstrap replications.

## Examples
```r
library(MASS)
library(clime)
library(tobitnet)

source("main.R")

## setting
n <- 200; p <- 100; alpha <- 0.05; B <- 500

## generate x, y, and c_q
Sigma <- toeplitz((1/2)^seq(0, p-1))
x <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
beta <- c(rep(sqrt(3), 3), rep(0, p-3))
y_star <- x %*% beta + rnorm(n, 0, 1)
var_y <- as.numeric(t(beta) %*% Sigma %*% beta) + 1
c_q <- sqrt(var_y) * qnorm(1/8)
y <- pmax(y_star, c_q)

# simultaneous inference
InferenceTheta(x, y, c_q, alpha, B)
InferenceBeta(x, y, c_q, alpha, B)

