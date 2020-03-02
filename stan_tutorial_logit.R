library(data.table)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)

#bart_pumps <- read.csv("bart_pumps.csv")
lotteries <- read.csv("lotteries.csv")

df.lotteries <- lotteries
df.lotteries$PX1 <- df.lotteries$PX1/100
df.lotteries$PX2 <- 1 - df.lotteries$PX1
df.lotteries$diffPX <- abs(df.lotteries$PX1 - df.lotteries$PX2)

df.lotteries$PZ1 <- df.lotteries$PZ1/100
df.lotteries$PZ2 <- 1 - df.lotteries$PZ1
df.lotteries$diffPZ <- abs(df.lotteries$PZ1 - df.lotteries$PZ2)

df.lotteries$evDiffX <- abs((df.lotteries$X1 * df.lotteries$PX1) - (df.lotteries$X2 * df.lotteries$PX2))
df.lotteries$evDiffZ <- abs((df.lotteries$Z1 * df.lotteries$PZ1) - (df.lotteries$Z2 * df.lotteries$PZ2))

df.lotteries$evDiff <- df.lotteries$evDiffX - df.lotteries$evDiffZ

df.lotteries_subset <- sample_n(df.lotteries, 1000, replace = TRUE)

### Set up 
y <- df.lotteries_subset$R
x <- df.lotteries_subset$evDiff
N <- as.numeric(nrow(df.lotteries_subset))

### Parse into list for Stan
df.lotteries_data <- list(N=N, x=x, y=y)
# specify model
logit_text <- 
  "data {
int<lower=0> N;   // number of data items
vector[N] x;   // predictor vector
int<lower = 0, upper = 1> y[N];      // outcome vector
}
parameters {
real alpha;           // intercept
real beta;       // coefficients for predictors
}
model {
y ~ bernoulli_logit(alpha + beta * x);  // likelihood
}
"
## Fit the model 
lotteries_fit <- stan(model_code=logit_text, data = df.lotteries_data,
                   verbose=TRUE)

summary(lotteries_fit)
