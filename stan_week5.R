### Load packages 
library(tidyr)
library(dplyr)
library(rstan)
library(parallel)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Set wd
setwd("~/myrepos/psyteam504")

## Load file
lotteries <- read.csv("lotteries.csv")
## Subset for fast
set.seed(08544)
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

tings_model <-
  'data {
  int N;
  int<lower=1> K; // the number of clusters to fit
  vector[N] y; // lotteries
}
parameters {
  real<lower=0> intercept[K];
  real y_coeff;
  real<lower=0> sigma[K]; // we need K sigmas because there are K components
  simplex[K] alpha; // mixture proportions aka weights
}
model {
  real mean_diff[K-1];
  real component_logprob[K];

  intercept ~ gamma(10,1);
  y_coeff ~ normal(0,20);
  
  sigma ~ cauchy(0, 2.5);
  alpha ~ dirichlet(rep_vector(2.0, K));
  
  for (k in 1:(K-1)) {
    mean_diff[k] = intercept[k+1] - intercept[k];
    mean_diff[k] ~ gamma(10, 1);
  }

  for (n in 1:N) {
    for (k in 1:K) {
      component_logprob[k] = 
          log(alpha[k]) + 
          normal_lpdf(y[n] | intercept[k], sigma[k]);
      }
      target += log_sum_exp(component_logprob);
    }
}'

# Fit Model
lotteries_multi_fit <- stan(model_code=multi_text, data = lotteries_data,
                      verbose=TRUE, chains = 1)

summary(lotteries_multi_fit)

lotteries_multi_fit2 <- stan(model_code=tings_model, data = lotteries_data,
                            verbose=TRUE, chains = 4, iter=2000, warmup=1000)
summary(lotteries_multi_fit2)
stan_trace(lotteries_multi_fit2)

ting_extract <- extract(lotteries_multi_fit2)
