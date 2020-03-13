### Load packages 
library(tidyr)
library(dplyr)
library(rstan)

## Set wd
setwd("~/myrepos/psyteam504")

###################################
##### Mixture model ###############
###################################

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
                            verbose=TRUE, chains = 1)

summary(lotteries_multi_fit)

##############################
##### Linear regression ######
##############################

# Gathering data
lotteries_theta <- as.data.frame(extract(lotteries_multi_fit, "theta"))
lotteries_subset_all <- merge(lotteries_subset, lotteries_theta)
lotteries_subset_all %>%
  select(partid, theta.1, theta.2) -> df.lotteries

bart_data <- read.csv("bart_pumps.csv", header = TRUE)
bart_subset <- sample_n(df.bart, 1000, replace = TRUE)

bart_subset %>%
  left_join(df.lotteries, by = "partid") %>%
  select(-theta.2) %>%
  rename(theta = theta.1) -> df.data

## Model
linear.mixed.effects.stan.prg = "
data {
  int<lower=0> N;   // number of data items
  int S; // unique participants
  vector[N] pumps; // outcome
  vector[N] explode;
  vector[N] theta;
  int partid[N]; // the participant id for each row
}
parameters {
  real intercept;
  real explode_coef;
  real theta_coef;

  // we don't really care about the random effects and it saves a lot of memory to put them here
  real intercept_adj[S];

  real<lower=0> sigma; // residual for the overall regression
  real<lower=0> intercept_rsigma; // random effect variance for intercept
}
model {
  // intercept ~ gamma(10, 1);
  // pumps ~ normal(0, 20); 
  // explode_coef ~ normal(0, 20);
  // theta_coef ~ lognormal(0, 2);
  // sigma ~ cauchy(0, 2.5);
  // intercept_rsigma ~ cauchy(0, 2.5);

  for (s in 1:S) {
    intercept_adj[s] ~ normal(0, intercept_rsigma);
  }

  for (n in 1:N) {
    pumps[n] ~ normal((intercept + intercept_adj[partid[n]]) + 
                       explode_coef * explode[n] + 
                       theta_coef * theta[n], sigma);
  }
}
"
df.data$partid = as.integer(factor(df.data$partid)) # it's important to use "factor" here to re-level the makes since we took a subset of the data

test <- na.omit(df.data)
fit.lmm <- stan(model_code = linear.mixed.effects.stan.prg,
                chains = 1,
                data = list(N = nrow(test), S = length(unique(test$partid)),
                            theta = test$theta, pumps = test$pumps, explode = test$explode,
                            partid = test$partid), verbose=TRUE)

summary(fit.lmm, pars=c('intercept', 'explode_coef', 'theta_coef', 'sigma', 'intercept_rsigma'))
