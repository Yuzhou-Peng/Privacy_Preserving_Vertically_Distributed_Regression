# Privacy Preserving Vertically Distributed Regression

## Outline

1.  Description
2.  Privacy Preserving Vertically Distributed Regression
3.  Package requirements
4.  Simulation
5.  Results

## Description

This README is prepared for journal peer review of the "hybrid statistical-epidemiological models for prediction in transmission disease with application to SARS-CoV-2" paper.

The proposed Quasi-Score online estimator is proposed for estimating the instantaneous reproduction number of an transmission disease that meets the basic assumptions of time-since-infection model with daily incident cases and covariates data. It allow covariates with measurement error to participate in the model while impose no distributional assumptions, and can update estimators online whenever new data are available.

To demonstrate its usage and performance (as in Section 5 of the paper), we use a simulated data (incident cases of a transmission disease) with time length T=120, initiated incident cases I_0=500 and two simulated covariates data to testify its performance on estimating the association of the instantaneous reproduction number of a transmission disease and potiential covariates. The current code version is written for a special case of the model (as described in Section 4 of the paper) with the link function h being the log function, I_t\|\mathcal{F}\_{t-1} \sim Poisson(R_t \Lambda\_t), the model structure admits AR(1) form and allows two covariates inputs (required to be invertible for its sample covariance matrix with days greater than two.)

Further, the code provides an bootstrap inspection on construct a bootstrap confidence interval for the parameters and a bootstrap confidence band for the instantaneous reproduction number with default tunning parameter block length tunningl=45.

QSOEID is a function coded in QSOEID.R file under [Covid-Quasi-Score](https://github.com/ChorusChow/Covid-Quasi-Score).
