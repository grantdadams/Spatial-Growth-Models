# This file creates a hierarchical von Bertalannfy growth function where parameter estimates are random effects of state and sex. The model assumes a left truncated normal distribution.

rm(list=ls())
# Load packages
library(TMB)

#..........................................................................
# Write model
# Remember the language is C++ so comments are indicated by "//" and the end of an expression needs to be explicitly stated by typing ";"
tmb_model = "
// #include <windows.h>
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //------------ DATA PREP ----------//
  DATA_INTEGER(n); // Number of obs
  DATA_INTEGER(J); // Number of states
  DATA_INTEGER(n_pred); // Number of predictors

  DATA_VECTOR(selectivity); // Minimum size for truncation
  DATA_VECTOR(age);  // Age data 
  DATA_VECTOR(y);  // Length data
  DATA_MATRIX(Pred);  // Predictor matrix

  DATA_IVECTOR(region);  // Vector of states

  DATA_MATRIX(W);
  DATA_MATRIX(D);
  
  //------------ PARAMETERS ----------//
  // Beta coefficients for mean of MVN
  PARAMETER_VECTOR(Betas_linf);
  PARAMETER_VECTOR(Betas_k);
  PARAMETER_VECTOR(Betas_t0);

  PARAMETER(log_sigma);
  PARAMETER_VECTOR(log_sigma_params);

  // Mean 0 random effects
  PARAMETER_VECTOR(log_linf_re);
  PARAMETER_VECTOR(log_k_re);
  PARAMETER_VECTOR(t0_re);

  // Transform variance components
  Type sigma=exp(log_sigma);
  Type sigma_sq=exp(2.0*log_sigma);

  vector<Type> sigma_params=exp(log_sigma_params);
  vector<Type> sigma_sq_params=exp(2.0*log_sigma_params);
  vector<Type> tau_vbgf=1/sigma_sq_params;

  // Spatial-correlation factors
  PARAMETER_VECTOR(log_phi_vbgf);
  vector<Type> phi_vbgf = 1/(1+exp(-log_phi_vbgf));

  //------------ MODEL ----------//

  // Integers of index variables and likelihood
  int i; // Initialize indexing variables
  Type nll=0; // Initialize negative log-likelihood

  // Create mean 0 MVN distributions for each parameter
  MVNORM_t<Type> N_0_Linf(matinv(tau_vbgf(0) * (D - phi_vbgf(0) * W)));
  MVNORM_t<Type> N_0_k(matinv(tau_vbgf(1) * (D - phi_vbgf(1) * W)));
  MVNORM_t<Type> N_0_t0(matinv(tau_vbgf(2) * (D - phi_vbgf(2) * W)));

  // Estimate region-specific random effects
  nll -= N_0_Linf(log_linf_re);
  nll -= N_0_k(log_k_re); 
  nll -= N_0_t0(t0_re); 

  // Fit model
  matrix<Type> VBGF_params(n,3); // Temporary matrix for vbgm params

  // Vectors of data and likelihood components
  vector<Type> LengthPred(n);
  vector<Type> CumulativeNorm(n);
  vector<Type> TruncatedNorm(n);

  int ind_row, ind_col; 
  
  for(i = 0; i < n; i++){
  
  ind_row = region(i);

  VBGF_params(i,0) = exp(log_linf_re(ind_row) + Pred(i,) * Betas_linf); // Fill vector of temporary Linf
  VBGF_params(i,1) = exp(log_k_re(ind_row) + Pred(i,) * Betas_linf); // Fill vector of temporary k
  VBGF_params(i,2) = t0_re(ind_row) + Pred(i,) * Betas_linf; // Fill vector of temporary t0

  LengthPred(i) = VBGF_params(i,0) * (1 - exp(-VBGF_params(i,1) * (age(i) - VBGF_params(i,2))));

  // Calculate the CDF for the right truncation minus the CDF for the left truncation
  CumulativeNorm(i) = log(1  - pnorm(selectivity(i),LengthPred(i),sigma));
  
  // Calculate the truncated nromal distribution
  TruncatedNorm(i) = dnorm(y(i),LengthPred(i) , sigma, true) - CumulativeNorm(i);
  
  // Calculate the negative log-likelihood
  nll -= TruncatedNorm(i); 
  }

  // Report transformed parameters
  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(t0);
  ADREPORT(sex_Linf);
  ADREPORT(sex_k);
  ADREPORT(sex_t0);
  ADREPORT(grand_Linf);
  ADREPORT(grand_k);
  ADREPORT(grand_t0);

  return nll;
}"
# Write model to C++
  write(tmb_model, file = "vbgf.cpp")
  
  # Load model template
  compile("vbgf.cpp") # Compile a C++ templated into a shared object file
  dyn.load(dynlib("vbgf"))
  
  #..........................................................................
  # Load data
  source("von Bertalanffy/VGBM Data Prep/VGBM_Data_Prep.R")
  
  # Data list
  Nsex = length(unique(dat$sex))
  J = length(unique(dat$state))
  dataList = list(y = dat$FL_mm, g = dat$state_no - 1, age = dat$fractional_age, n = nrow(dat), J = J, sex = dat$sex - 1, Nsex = Nsex, selectivity = (dat$selectivity), source = dat$source)
  
  # Parameter list
  parametersList <- list(log_sigma = 1, 
                         log_sigma_params = c(1,1,1), 
                         sex_log_sigma_params = c(1,1,1), 
                         Linf_log = matrix(rep(6, J * Nsex), ncol = Nsex, nrow = J), 
                         k_log = matrix(rep(log(.3), J * Nsex), ncol = Nsex, nrow = J),
                         t0_log = matrix(rep(1, J * Nsex), ncol = Nsex, nrow = J),
                         sex_Linf_log = rep(6.2, Nsex),
                         sex_k_log = rep(log(0.3), Nsex),
                         sex_t0_log = rep(1, Nsex),
                         grand_Linf_log = 6,
                         grand_k_log = log(0.3),
                         grand_t0_log = 1,
                         B_k = 0,
                         B_Linf = 0)
  
  # Run the model without the truncated normal
  obj = MakeADFun(
    data = dataList, 
    parameters = parametersList,
    DLL = "vbgfn",
    random = c("sex_Linf_log", "sex_k_log", "sex_t0_log", "Linf_log", "k_log", "t0_log")) # Prepar
  opt <- nlminb(obj$par, obj$fn, obj$gr) # Fit model using nlminb
  rep <- sdreport(obj)
  results <- summary(rep)
  
  # Get the truncated norm inits
  # Parameter list
  parametersListTN <- list(log_sigma = results[which(rownames(results) == "log_sigma"),1], 
                         log_sigma_params = results[which(rownames(results) == "log_sigma_params"),1], 
                         sex_log_sigma_params = results[which(rownames(results) == "sex_log_sigma_params"),1], 
                         Linf_log = matrix(results[which(rownames(results) == "Linf_log"),1], ncol = Nsex, nrow = J), 
                         k_log = matrix(results[which(rownames(results) == "k_log"),1], ncol = Nsex, nrow = J),
                         t0_log = matrix(results[which(rownames(results) == "t0_log"),1], ncol = Nsex, nrow = J),
                         sex_Linf_log = results[which(rownames(results) == "sex_Linf_log"),1],
                         sex_k_log = results[which(rownames(results) == "sex_k_log"),1],
                         sex_t0_log = results[which(rownames(results) == "sex_t0_log"),1],
                         grand_Linf_log = results[which(rownames(results) == "grand_Linf_log"),1],
                         grand_k_log = results[which(rownames(results) == "grand_k_log"),1],
                         grand_t0_log = results[which(rownames(results) == "grand_t0_log"),1],
                         B_k = results[which(rownames(results) == "B_k"),1],
                         B_Linf = results[which(rownames(results) == "B_Linf"),1])
  
  # Run the truncated normal model
  obj = MakeADFun(
    data = dataList, 
    parameters = parametersListTN,
    DLL = "vbgf",
    random = c("sex_Linf_log", "sex_k_log", "sex_t0_log", "Linf_log", "k_log", "t0_log")) # Prepar
  opt <- nlminb(obj$par, obj$fn, obj$gr) # Fit model using nlminb
  rep <- sdreport(obj)
  results <- summary(rep)
  
  
  # Sample the posterior
  library(adnuts)
  init <- function() parametersListTN
  obj = MakeADFun(
    data = dataList, 
    parameters = parametersListTN,
    DLL = "vbgf") # 
  fit1 <- sample_tmb(obj=obj, iter=1000, chains=3, init=init, seeds=1:3,
                     control=list(metric=rep$covar))

  