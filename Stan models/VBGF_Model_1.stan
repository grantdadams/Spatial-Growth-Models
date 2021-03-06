data {
  // ------------------------------------------------------------------------------------ //
  // 1. Assign data to stan objects                                                       //
  // ------------------------------------------------------------------------------------ //
  int<lower = 1>      n_i;            // Sample size
  int<lower = 1>      n_pred;         // Number of predictors
  vector[n_i]         log_length_i;   // Vector of log fork length of fish i
  vector[n_i]         age_i;          // Vector of age of fish i
  matrix[n_i, n_pred] design_mat;     // n_i * n_pred design matrix of linear predictors
}
parameters {
  // ------------------------------------------------------------------------------------ //
  // 2. Specify model parameters                                                          //
  // ------------------------------------------------------------------------------------ //
  // 2.1. Regression coefficients for mean of VBGF parameters
  vector[n_pred] B_log_linf;
  vector[n_pred] B_log_k;
  vector[n_pred] B_t0;
  
  // 2.2. Random error (SD)
  real<lower=0> sigma;
}
model {
  // ------------------------------------------------------------------------------------ //
  // 3. Specify likelihood and priors                                                     //
  // ------------------------------------------------------------------------------------ //
  // 3.1. Temporary model objects to save parameter vectors
  vector[n_i] log_linf_i ;
  vector[n_i] k_i ;
  vector[n_i] t0_i ;
  
  // 3.2. Priors
  // 3.2.1. -- Regression coefficients for mean of VBGF parameters
  B_log_linf  ~ normal(0,10);  
  B_log_k     ~ normal(0,10);
  B_t0        ~ normal(0,10);
  
  // 3.2.2 -- Error components
  sigma       ~ cauchy(0, 5);
  
  // 3.3. Model specification
  // 3.3.1. -- VBGF Parameters
  log_linf_i = design_mat * B_log_linf;
  k_i        = exp(design_mat * B_log_k);
  t0_i       = design_mat * B_t0; 
  
  // 3.3.2. -- Model likelihood
  log_length_i ~ normal( log_linf_i + log1m_exp(-k_i .* (age_i - t0_i)), sigma);
}
generated quantities {
  vector[n_i] log_lik;
  for(i in 1:n_i)
  log_lik[i] = normal_lpdf(log_length_i[i] | (design_mat[i,] * B_log_linf) + log1m_exp(-exp(design_mat[i,] * B_log_k) .* (age_i[i] - (design_mat[i,] * B_t0))), sigma);
}
