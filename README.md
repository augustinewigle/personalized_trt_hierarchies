# Personalized Treatment Hierarchies in Bayesian Network Meta-Analysis

This repo contains the code and data needed to reproduce the results in the article. 
The individual participant data from STAR*D, EMBARC, and REVAMP studies are confidential and can not be shared. Instead, point estimates of the study-level parameters and their covariances are available in .csv files.
This repository contains the following files:
  - README.md: this file
  - coefficients_embarc.csv, coefficients_revamp.csv, coefficients_stard.csv: The point estimates of the parameters in each study
  - vcov_embarc.csv, vcov_revamp.csv, vcov_stard.csv: The estimated variance-covariance matrices of the parameter estimates from each study
  - fullcovarmatrix.stan: Stan code for the Bayesian NMA model to pool study-level estimates
  - stage2.R: R script to process the .csv files, fit the Bayesian NMA model, and create personalized hierarchies and plots
