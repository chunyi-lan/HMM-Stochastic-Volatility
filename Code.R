# MATH6604 Project (Fall 2023)
# York University, Department of Mathematics & Statistics
# Title: Likelihood calculations for Hidden Markov Models in financial mathematics
# Author: Noam Tobin & Chun Yi Lan
# 
# Dataset: 1970-01-01 to 2023-11-25 S&P 500 Daily Closing Price

library(dplyr)
library(astsa)
library(lubridate)
library(stats4)

# Initializing & Descriptive Statistics

SP500 <- read.csv('SPDataset(1970-20231125).csv', header = TRUE)
ClosingPrice <- ts(SP500$Close, start = decimal_date(as.Date("1970-01-01")), frequency = 365)

ts.plot(ClosingPrice, col="#4682B4", main="S&P Daily Closing Price", ylab="Price (Dollars)")

# Daily returns
return <- diff(log(ClosingPrice))

ts.plot(return, col="#4682B4", main="S&P Daily Return", ylab="Return")

summary(ClosingPrice)
summary(return)

# Modelling

# Initial parameter values
initial_params <- c(phi = 0.98, sigma = 0.2, beta = 0.05)

# Set other parameters
m2 <- 2
m10 <- 10
m50 <- 50
m100 <- 100
gmax <- 4
y <- return

# Calculate the log-likelihood of the HMM approx to the SV0 model
likelihood.SV0_reg <- function(y, phi, sigma, beta, m, gmax, reg_term) {
  K <- m + 1
  b <- seq(-gmax, gmax, length = K)
  bs <- (b[-1] + b[-K]) * 0.5
  sey <- beta * exp(bs / 2)
  
  # Check for near-zero standard deviations and replace with a small value
  sey[sey < 1e-6] <- 1e-6 # Avoid having a zero standard deviation as will result in an error when calculating pdf of normal (dnorm) for dividing 0
  
  # Initialize the t.p.m. Gamma
  Gamma <- matrix(0, m, m)
  for (i in 1:m) Gamma[i,] <- diff(pnorm(b, phi * bs[i], sigma))
  Gamma <- Gamma / apply(Gamma, 1, sum) # scale each rows sum to 1
  
  foo <- solve(t(diag(m) - Gamma + 1 + reg_term * diag(m)), rep(1, m))
  llk <- 0
  for (t in 1:length(y)) {
    foo <- foo %*% Gamma * dnorm(y[t], 0, sey)
    sumfoo <- sum(foo)
    if (sumfoo < 1e-6) sumfoo <- 1e-6 # Avoid taking log(0)
    llk <- llk + log(sumfoo)
    foo <- foo / sumfoo
  }
  return(llk)
}

# Define the negative log-likelihood function
neg_log_likelihood <- function(params) {
  phi <- params[1]
  sigma <- params[2]
  beta <- params[3]
  return(-likelihood.SV0_reg(y, phi, sigma, beta, m, gmax, reg_term = 1e-6))
}

# Optimize the negative log-likelihood using the Nelder-Mead method
m <- m2
result_SV0_m2 <- optim(par = initial_params, fn = neg_log_likelihood, method = "Nelder-Mead")
m <- m10
result_SV0_m10 <- optim(par = initial_params, fn = neg_log_likelihood, method = "Nelder-Mead")
m <- m50
result_SV0_m50 <- optim(par = initial_params, fn = neg_log_likelihood, method = "Nelder-Mead")
m <- m100
result_SV0_m100 <- optim(par = initial_params, fn = neg_log_likelihood, method = "Nelder-Mead")

# Extract the optimized parameters
params_SV0_m2 <- result_SV0_m2$par
params_SV0_m10 <- result_SV0_m10$par
params_SV0_m50 <- result_SV0_m50$par
params_SV0_m100 <- result_SV0_m100$par
params_SV0_m2
params_SV0_m10
params_SV0_m50
params_SV0_m100

# AIC
k <- length(params_SV0_m2)  # No. of parameters in the model
n <- length(y)

max_llk_SV0_m2 <- -result_SV0_m2$value  # Maximised log-likelihood (negative because we are minimizing)
aic_SV0_m2 <- 2 * k - 2 * max_llk_SV0_m2
bic_SV0_m2 <- k * log(n) - 2 * max_llk_SV0_m2
aic_SV0_m2
bic_SV0_m2

max_llk_SV0_m10 <- -result_SV0_m10$value
aic_SV0_m10 <- 2 * k - 2 * max_llk_SV0_m10
bic_SV0_m10 <- k * log(n) - 2 * max_llk_SV0_m10
aic_SV0_m10
bic_SV0_m10

max_llk_SV0_m50 <- -result_SV0_m50$value
aic_SV0_m50 <- 2 * k - 2 * max_llk_SV0_m50
bic_SV0_m50 <- k * log(n) - 2 * max_llk_SV0_m50
aic_SV0_m50
bic_SV0_m50

max_llk_SV0_m100 <- -result_SV0_m100$value
aic_SV0_m100 <- 2 * k - 2 * max_llk_SV0_m100
bic_SV0_m100 <- k * log(n) - 2 * max_llk_SV0_m100
aic_SV0_m100
bic_SV0_m100

# Backtesting

simulate_HMM_SV0_model <- function(params, n) {
  beta <- params[1]
  phi <- params[2]
  sigma <- params[3]
  
  # Initialize arrays to store results
  X <- numeric(n)
  Y <- numeric(n)
  
  # Initial value for log-volatility
  X[1] <- rnorm(1)
  
  for (t in 2:n) {
    # Generate innovations
    xi1 <- rnorm(1)
    xi2 <- rnorm(1)
    
    # Update log-volatility
    X[t] <- phi * X[t-1] + sigma * xi2
    
    # Generate observation
    Y[t] <- beta * exp(X[t]/2) * xi1
  }
  
  return(Y)
}

alpha <- 0.01 # 1% VaR

# Simulate the model form the optimized param
simulated_returns <- simulate_HMM_SV0_model(params_SV0_m100, length(y))

# Calculate VaR at Conf Lvl
var_threshold <- quantile(simulated_returns, alpha)

# Identify exceptions
exceptions <- y < var_threshold

# No. of exceptions
n_exceptions <- sum(exceptions)

# Binomial dist'n parameters
binomial_params <- c(n = length(y), alpha = alpha)

# 95th & 99.99th Percentiles for Binomial
green_zone_threshold <- qbinom(0.95, size = binomial_params[1], prob = binomial_params[2])
yellow_zone_lower_threshold <- qbinom(0.95, size = binomial_params[1], prob = binomial_params[2] + 0.005)
yellow_zone_upper_threshold <- qbinom(0.9999, size = binomial_params[1], prob = binomial_params[2])

# Evaluate model performance based on the no. of exceptions
if (n_exceptions < green_zone_threshold) {
  cat("Model is accurate (Green Zone)\n")
} else if (n_exceptions > yellow_zone_upper_threshold) {
  cat("Model is inaccurate (Red Zone)\n")
} else if (n_exceptions > yellow_zone_lower_threshold) {
  cat("Model requires additional information (Yellow Zone)\n")
} else {
  cat("Model performance is within acceptable limits\n")
}