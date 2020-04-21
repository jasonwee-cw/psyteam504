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
############### HMM ###############
###################################

# Gathering data
bart <- read.csv("bart_pumps.csv", header = TRUE)

# Parameters
K = 2 # hypothesized number of hidden states
N = nrow(bart) #number of instances
V = bart$pumps
T = bart$trial

bart_data <- list(K=K, N=N, V=V, T=T)

# Model
hmm.model = "
data {
  int<lower=1> K;  // num categories
  int<lower=0> N;   // number of data items
  real y[N]
  
  int<lower=1,upper=K> z[N]; // categories
  
  vector<lower=0>[K] alpha;  // transit prior
  vector<lower=0>[V] beta;   // emit prior
}
parameters {
  simplex[K] theta[K];  // transit probs
  simplex[pumps] phi[K];    // emit probs
}
model {
  for (k in 1:K)
    theta[k] ~ dirichlet(alpha);
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);
  for (n in 1:N)
    pumps[n] ~ categorical(phi[z[n]]);
  for (t in 2:N)
    z[n] ~ categorical(theta[z[n - 1]]);
}
"
hmm.fit <- stan(model_code=hmm.model, data = bart_data,
                            verbose=TRUE, chains = 1)






