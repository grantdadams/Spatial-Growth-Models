data {
  int<lower = 1> n; // Sample size
  int<lower = 1> J;// Number of regions
  int<lower = 1> n_pred;
  vector[n] log_y;
  vector[n] age;
  matrix[n, n_pred] Pred;
  int<lower = 0> region[n];
}
parameters {
  // Beta coefficients for mean of VBGF parameters
  vector[n_pred] Betas_log_linf;
  vector[n_pred] Betas_log_k;
  vector[n_pred] Betas_t0;
  
  // SD
  real<lower = 0> sigma;
  
  // Scaling parameters for non-centered distribution
  vector[J] alpha_linf;
  vector[J] alpha_k;
  vector[J] alpha_t0;
  
  // Hierarchical Variance
  real<lower = 0> sigma_linf;
  real<lower = 0> sigma_k;
  real<lower = 0> sigma_t0;
}
transformed parameters{
  // Random effects
  vector[J] log_linf_re;
  vector[J] log_k_re;
  vector[J] t0_re;
  
  // Get random effects
  log_linf_re = sigma_linf * alpha_linf; // equivalent to ~ N(0, sigma)
  log_k_re = sigma_k * alpha_k;
  t0_re = sigma_t0 * alpha_t0;
}
model {
  // Parameters
  vector[n] log_linf ;
  vector[n] k ;
  vector[n] t0 ;
  
  // -----------------------------------------------------
  // Priors
  // -----------------------------------------------------
  Betas_log_linf ~ normal(0,10);  
  Betas_log_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  sigma ~ cauchy(0,5);
  sigma_linf ~ cauchy(0,5);
  sigma_k ~ cauchy(0,5);
  sigma_t0 ~ cauchy(0,5);
  
  // Scaling factors
  alpha_linf ~ normal(0,1);   
  alpha_k ~ normal(0,1);
  alpha_t0 ~ normal(0,1);
  
  // Parameters
  log_linf = Pred * Betas_log_linf + log_linf_re[region];
  k = exp(Pred * Betas_log_k + log_k_re[region]);
  t0 = Pred * Betas_t0 + t0_re[region]; 
  
  // Likelihood
  log_y ~ normal( log_linf + log1m_exp(-k .* (age - t0)), sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(log_y[i] | (Pred[i,] * Betas_log_linf + log_linf_re[region[i]]) + log1m_exp(-exp(Pred[i,] * Betas_log_k + log_k_re[region[i]]) .* (age[i] - (Pred[i,] * Betas_t0 + t0_re[region[i]]))), sigma);
}

