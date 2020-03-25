### PSY504 Psych Team
### Crystal Lee & Leon Mait
### Created: 03/13/20
### Modified: 03/24/20
### Last edited by: Leon 

### Load packages 
library(tidyr)
library(dplyr)
library(rstan)
#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

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
save(summary.lottery, file = "summary.lottery.RData")


## Calculate posterior cluster probabilities per participant
theta1 <- as.data.frame(summary.lottery$summary)$mean[1]
mu1 <- as.data.frame(summary.lottery$summary)$mean[3]
sd1 <- as.data.frame(summary.lottery$summary)$mean[5]

lotto %>%
  group_by(partid) %>%
  mutate(cluster1 = pnorm(R, mu1, sd1)) -> lotteries_theta

lotto %>%
  group_by(partid) %>%
  mutate(cluster1 = pnorm(R, mu1, sd1),
         cluster2 = 1-cluster1) %>%
  pivot_longer(cols = cluster1:cluster2,
               names_to = "Cluster",
               values_to = "Cluster Probability") -> lotteries_theta_plot

# Plots
stan_hist(lotteries_multi_fit, pars = "mu")
traceplot(lotteries_multi_fit, pars = c("theta"))
ggsave("lotteries_traceplot_theta.png")
traceplot(lotteries_multi_fit, pars = c("mu"))
ggsave("lotteries_traceplot_mu.png")
traceplot(lotteries_multi_fit, pars = c("sigma"))
ggsave("lotteries_traceplot_sigma.png")

ggplot(lotteries_theta_plot, aes(x = R, y = `Cluster Probability`, color = Cluster)) +
  geom_line() +
  xlab("R factor score") +
  ylab("Probability") +
  scale_color_discrete(breaks = c("cluster1", "cluster2"), 
                       labels = c("Cluster 1", "Cluster 2")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("clusters.png")

############################################
##### Linear mixed effects regression ######
############################################

# Gathering data
# lotteries_theta <- as.data.frame(extract(lotteries_multi_fit, "theta"))
# nrow(lotteries_theta)
# 
# hist(lotteries_theta$theta.1)
# hist(lotteries_theta$theta.2)
# 
# lotto.sum <- summary(lotteries_multi_fit)
# lotto.sum$summary

# Gathering data
bart_data <- read.csv("bart_pumps.csv", header = TRUE)

lotteries_theta %>%
  left_join(bart_data, by = "partid") %>%
  na.omit() %>%
  arrange(partid, trial, block) -> df.data

## Model
linear.mixed.effects.stan.prg = "
data {
  int<lower=0> N;   // number of data items
  int S; // unique participants
  vector[N] pumps; // outcome
  vector[N] theta;
  int partid[N]; // the participant id for each row
  }
parameters {
  real intercept;
  real theta_coef;
  real<lower=0> sigma; // residual for the overall regression

  real intercept_adj[S];
  real<lower=0> intercept_rsigma; // random effect variance for intercept
}
model {
  intercept ~ gamma(10, 1);
  pumps ~ gamma(10, 1); 
  theta_coef ~ normal(0, 20);
  sigma ~ cauchy(0, 2.5);
  intercept_rsigma ~ cauchy(0, 2.5);
  for (s in 1:S) {
    intercept_adj[s] ~ normal(0, intercept_rsigma);
  }
  for (n in 1:N) {
    pumps[n] ~ normal((intercept + intercept_adj[partid[n]]) + 
                       theta_coef * theta[n], sigma);
  }
}
"
df.data$partid = as.integer(factor(df.data$partid)) # it's important to use "factor" here to re-level the makes since we took a subset of the data

fit.lmm <- stan(model_code = linear.mixed.effects.stan.prg,
                chains = 3,
                data = list(N = nrow(df.data), S = length(unique(df.data$partid)),
                            theta = df.data$cluster1, pumps = df.data$pumps,
                            partid = df.data$partid), verbose=TRUE)
save(fit.lmm, file = "fitlmm.RData")


summary(fit.lmm, pars=c('intercept', 'theta_coef', 'sigma', 'intercept_rsigma'))

# Plots
stan_plot(fit.lmm, pars = c("intercept", 'intercept_rsigma', "sigma", "theta_coef"), show_density = TRUE, fill_color = "blue")
ggsave("fitlmm.png")

traceplot(fit.lmm, pars = c("theta_coef"))
ggsave("lmm_traceplot_theta.png")

######################################
##### Autoregressive regression ######
######################################

arr.model = "
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
  real pumps_coef;
  real<lower=0> sigma; // residual for the overall regression

  real intercept_adj[S];
  real<lower=0> intercept_rsigma; // random effect variance for intercept
}
model {
  intercept ~ gamma(10, 1);
  pumps ~ gamma(10, 1); 
  explode_coef ~ normal(0, 20);
  theta_coef ~ normal(0, 20);
  sigma ~ cauchy(0, 2.5);
  intercept_rsigma ~ cauchy(0, 2.5);
  for (s in 1:S) {
    intercept_adj[s] ~ normal(0, intercept_rsigma);
  }
  for (n in 2:N) {
    pumps[n] ~ normal((intercept + intercept_adj[partid[n]]) + 
                       explode_coef * explode[n-1] + 
                       theta_coef * theta[n], sigma);
  }
}
"

arr.fit <- stan(model_code = arr.model,
                chains = 3,
                data = list(N = nrow(df.data), S = length(unique(df.data$partid)),
                            theta = df.data$cluster1, pumps = df.data$pumps,
                            partid = df.data$partid, explode = df.data$exploded), verbose=TRUE)

save(arr.fit, file = "arrfit.RData")

summary(arr.fit, pars=c('intercept', 'explode_coef', 'theta_coef', 'sigma', 'intercept_rsigma'))

# Plots
stan_plot(arr.fit, pars = c("intercept", "sigma", "theta_coef", "explode_coef"), show_density = TRUE, fill_color = "blue")
ggsave("arrfit.png")

traceplot(arr.fit, pars = c("theta_coef"))
ggsave("arr_traceplot_theta.png")
traceplot(arr.fit, pars = c("explode_coef"))
ggsave("arr_traceplot_explode.png")
