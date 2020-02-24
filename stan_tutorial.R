#Set working directory
setwd("~/Documents/Princeton/PSY504 Statistics")
library(data.table)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)


#Hi @everyone, for the OLS stan practice, please use the mtcars data linked here: https://gist.github.com/seankross/a412dfbd88b3db70b74b  This dataset is included with R so there's no need to download it if you use R
#For the OLS regression, simply build a model that predicts "mpg" based on "hp" (horsepower) and "wt" (weight) of the vehicle.
#mpq ~ hp + wt

## Load mtcars data
mtcars <- mtcars
### Set up 
y <- mtcars$mpg
x <- matrix(c(mtcars$hp, mtcars$wt),ncol=2)
N <- as.numeric(nrow(mtcars))
K <- ncol(x)
### Parse into list for Stan
mtcars_data <- list(N=N, K=K, x=x, y=y)
# specify model
ols_text <- 
  "data {
int<lower=0> N;   // number of data items
int<lower=0> K;   // number of predictors
matrix[N, K] x;   // predictor matrix
vector[N] y;      // outcome vector
}
parameters {
real alpha;           // intercept
vector[K] beta;       // coefficients for predictors
real<lower=0> sigma;  // error scale
}
model {
y ~ normal(x * beta + alpha, sigma);  // likelihood
}
"
## Fit the model 
mtcars_fit <- stan(model_code=ols_text, data = mtcars_data,
                   verbose=TRUE)

summary(mtcars_fit)
