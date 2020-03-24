### PSY504 Psych Team
### Crystal Lee & Leon Mait
### Created: 03/13/20
### Modified: 03/24/20
### Last edited by: Leon 

## Load packages 
library(tidyr)
library(dplyr)
library(rstan)

## Set wd
setwd("~/myrepos/psyteam504")

###################################
##### Mixture model ###############
###################################

## Load file
lotto <- read.csv("lotteriesOvert.csv")
nrow(lotto) #1507

# R = Proportion of risky (= higher variance) choices
# H = Proportion of higher expected value choices (EV = probability of payoff (%) * payoff amount ($))
# CV = Proportion of higher coefficient of variance choices (CoV = the dispersion of data points in a data series around the mean)

## Set up parameters
y <- lotto$R # DV
K <- 2 # No of groups
N <- as.numeric(nrow(lotto))

## Parse into list for Stan
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

## Fit Model
lotteries_multi_fit <- stan(model_code=multi_text, data = lotteries_data,
                            verbose=TRUE, chains = 1)

summary.lottery <- summary(lotteries_multi_fit)

## Calculate posterior cluster probabilities per participant
mu1 <- as.data.frame(summary.lottery$summary)$mean[3]
sd1 <- as.data.frame(summary.lottery$summary)$mean[5]

lotto %>%
  group_by(partid) %>%
  mutate(cluster1 = pnorm(R, mu1, sd1)) -> lotteries_theta

##############################
##### Linear regression ######
##############################

# Gathering data
# lotteries_theta <- as.data.frame(extract(lotteries_multi_fit, "theta"))
# nrow(lotteries_theta)
# 
# hist(lotteries_theta$theta.1)
# hist(lotteries_theta$theta.2)
# 
# lotto.sum <- summary(lotteries_multi_fit)
# lotto.sum$summary

## Gathering data
bart_data <- read.csv("bart_pumps.csv", header = TRUE)
bart_subset <- sample_n(bart_data, 1000, replace = FALSE) 

lotteries_theta %>%
  left_join(bart_subset, by = "partid") %>%
  na.omit() %>%
  arrange(partid, trial, block) -> df.data

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
