# Hidden Markov Models for Stochastic Volatility Modelling

This repository contains the implementation of Hidden Markov Models (HMMs) for stochastic volatility (SV) modeling, with a specific focus on financial mathematics. The project is based on the research paper by Roland Langrock, Iain L. MacDonald, and Walter Zucchini, and applies HMM approximations to analyze S&P 500 index data. This project is part of the coursework for York University's Mathematics graduate course, MATH 6604 Probability Models.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Requirements](#requirements)
- [Usage](#usage)
- [References](#references)

## Introduction
Stochastic volatility models describe volatility (variance) in financial assets as a stochastic process, addressing the limitations of traditional models that assume constant volatility. This project demonstrates how HMMs can effectively approximate stochastic volatility models (e.g., SV0, SVt) for time series data, enabling parameter estimation, forecasting, and risk assessment in financial markets.

Key highlights include:
- Likelihood calculation for HMM approximations of stochastic volatility models.
- Application to historical S&P 500 index data.
- Backtesting and risk estimation using Value at Risk (VaR).

## Features
- Implementation of HMM approximations for SV models.
- Parameter estimation using optimization techniques.
- Goodness-of-fit evaluation with AIC and BIC.
- Backtesting with Value at Risk (VaR) analysis.

## Requirements
The project is implemented in R. Below are the required packages:
- `dplyr`
- `astsa`
- `lubridate`
- `stats4`

Ensure you have R installed on your system. You can install the required packages using the following command:
```R
install.packages(c("dplyr", "astsa", "lubridate", "stats4"))
```

## Usage
1. **Prepare the Dataset**:  
   - Ensure you have the daily closing prices of the S&P 500 Index from January 1, 1970, to November 25, 2023. Update the dataset path in the R scripts if necessary.

2. **Run the Analysis**:  
   - Open the R scripts provided in the repository.
   - Execute these scripts step-by-step in R or RStudio to perform the analysis:
     - Data preprocessing.
     - Parameter estimation for HMM-based stochastic volatility (SV) models.
     - Goodness-of-fit evaluation using AIC and BIC.
     - Backtesting with Value at Risk (VaR).

3. **View Results**:  
   - The log-likelihood, AIC/BIC values, and backtesting results will be displayed in the console. Generated plots and other outputs will be saved in the results directory.

## References
1. Roland Langrock, Iain L. MacDonald, and Walter Zucchini:  
   "Some nonstandard stochastic volatility models and their estimation using structured hidden Markov models."  
   *Journal of Empirical Finance* (2012) 19: 147-161.  
   [Link to paper](https://www.sciencedirect.com/science/article/pii/S0927539811000661)

2. Walter Zucchini, Iain L. MacDonald, and Roland Langrock:  
   *Hidden Markov Models for Time Series: An Introduction Using R.*  
   2nd ed. Routledge, 2017. ISBN: 9781032179490.

3. S&P 500 Index - Google Finance.  
   Accessed: November 26, 2023.  
   [Link to data](https://www.google.com/finance/quote/.INX:INDEXSP?hl=en)

