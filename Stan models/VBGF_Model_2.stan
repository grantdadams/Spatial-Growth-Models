data {
  // ------------------------------------------------------------------------------------ //
  // 1. Assign data to stan objects                                                       //
  // ------------------------------------------------------------------------------------ //
  int<lower = 1>      n_i;            // Sample size
  int<lower = 1>      n_r;            // Number of regions
  int<lower = 1>      n_pred;         // Number of predictors
  vector[n_i]         log_length_i;   // Vector of log fork length of fish i
  vector[n_i]         age_i;          // Vector of age of fish i
  matrix[n_i, n_pred] design_mat;     // n_i * n_pred design matrix of linear predictors
  int<lower = 0>      r_i[n_i];       // Integer vector of region of fish i
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
  
  // 2.3. Hierarchical VBGF parameter error (SD)
  real<lower = 0> sigma_linf;
  real<lower = 0> sigma_k;
  real<lower = 0> sigma_t0;

  // 2.4. Scaling parameters for non-centered distribution
  vector[n_r] alpha_linf;
  vector[n_r] alpha_k;
  vector[n_r] alpha_t0;
}
transformed parameters{
  // ------------------------------------------------------------------------------------ //
  // 3. Specify derived model parameters to be saved                                      //
  // ------------------------------------------------------------------------------------ //
  // 3.1. VBGF parameter-specific random effects
  vector[n_r] log_linf_re;
  vector[n_r] log_k_re;
  vector[n_r] t0_re;
  
  // 3.2. Get random effects
  log_linf_re = sigma_linf * alpha_linf; // Equivalent log_linf_re ~ normal(0, sigma_linf)
  log_k_re    = sigma_k    * alpha_k;
  t0_re       = sigma_t0   * alpha_t0;
}
model {
  // ------------------------------------------------------------------------------------ //
  // 4. Specify likelihood and priors                                                     //
  // ------------------------------------------------------------------------------------ //
  // 4.1. Temporary model objects to save parameter vectors
  vector[n_i] log_linf_i ;
  vector[n_i] k_i ;
  vector[n_i] t0_i ;
  
  // 4.2. Priors
  // 4.2.1. -- Regression coefficients for mean of VBGF parameters
  B_log_linf  ~ normal(0,10);  
  B_log_k     ~ normal(0,10);
  B_t0        ~ normal(0,10);
  
  // 4.2.2 -- Error components
  sigma       ~ cauchy(0, 5);
  sigma_linf  ~ cauchy(0, 5);
  sigma_k     ~ cauchy(0, 5);
  sigma_t0    ~ cauchy(0, 5);
  
  // 4.2.3 -- Scaling factors for non-centered parameterization
  alpha_linf  ~ normal(0,1);
  alpha_k     ~ normal(0,1);
  alpha_t0    ~ normal(0,1);
  
  // 4.3. Model specification
  // 4.3.1. -- VBGF Parameters
  log_linf_i = design_mat * B_log_linf + log_linf_re[r_i];
  k_i        = exp(design_mat * B_log_k + log_k_re[r_i]);
  t0_i       = design_mat * B_t0 + t0_re[r_i]; 
  
  // 4.3.2. -- Model likelihood
  log_length_i ~ normal( log_linf_i + log1m_exp(-k_i .* (age_i - t0_i)), sigma);
}
generated quantities {
  vector[n_i] log_lik;
  for(i in 1:n_i)
  log_lik[i] = normal_lpdf(log_length_i[i] | (design_mat[i,] * B_log_linf + log_linf_re[r_i[i]]) + log1m_exp(-exp(design_mat[i,] * B_log_k + log_k_re[r_i[i]]) .* (age_i[i] - (design_mat[i,] * B_t0 + t0_re[r_i[i]]))), sigma);
}
