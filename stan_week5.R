### Load packages 
library(tidyr)
library(dplyr)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Set wd
setwd("~/myrepos/psyteam504")

## Load file
lotteries <- read.csv("lotteries.csv")
## Subset for fast
lotteries_subset <- sample_n(lotteries, 1000, replace = TRUE)


## Set up parameters
y <- lotteries_subset$R # DV
K <- 2 # No of groups
N <- as.numeric(nrow(lotteries_subset))

### Parse into list for Stan
lotteries_data <- list(K=K, N=N, y=y)


## Specify Model
multi_text <-
  "data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
}
parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] mu;             // locations of mixture components
  vector<lower=0>[K] sigma;  // scales of mixture components
}
model {
  vector[K] log_theta = log(theta);  // cache log calculation
  sigma ~ lognormal(0, 2);
  mu ~ normal(0, 10);
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K)
      lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    target += log_sum_exp(lps);
  }
}
"

# Fit Model
lotteries_multi_fit <- stan(model_code=multi_text, data = lotteries_data,
                      verbose=TRUE)

summary(lotteries_multi_fit)

theta1 <- extract(lotteries_multi_fit, 'theta[1]', permuted=FALSE)
