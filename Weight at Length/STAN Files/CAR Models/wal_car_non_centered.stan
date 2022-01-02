data {
  int<lower = 1> n;
  int<lower = 1> J;
  int<lower = 1> n_pred;
  vector[n] y;
  vector[n] age;
  real selectivity[n];
  
  matrix[n, n_pred] Pred;
  int<lower = 0> region[n];
  matrix<lower = 0, upper = 1>[J, J] W;
  matrix<lower = 0>[J, J] D;
}
parameters {
  // Scaling parameters for non-centered distribution
  vector[J] alpha_linf;
  vector[J] alpha_k;
  vector[J] alpha_t0;

  // Beta coefficients for mean of MVN
  vector[n_pred] Betas_log_linf;
  vector[n_pred] Betas_log_k;
  vector[n_pred] Betas_t0;
  
  real<lower=0> sigma;
  
  // Hierarchical Variance
  real<lower = 0> tau_linf;
  real<lower = 0> tau_k;
  real<lower = 0> tau_t0;
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_linf;
  real<lower = -1, upper = 1> phi_k;
  real<lower = -1, upper = 1> phi_t0;
}
transformed parameters{
  // Lower triangle matrices of precision
  matrix[J,J] L_linf;
  matrix[J,J] L_k;
  matrix[J,J] L_t0;
  
  // Random effects
  vector[J] log_linf_re;
  vector[J] log_k_re;
  vector[J] t0_re;
  
  // Get cholesky decomp of precision
  L_linf= tau_linf * (D - phi_linf * W);
  L_k = tau_k * (D - phi_k * W);
  L_t0 = tau_t0 * (D - phi_t0 * W);
  
  // Get random effects
  log_linf_re = L_linf \ alpha_linf; // equivalent to inverse(L) * alpha
  log_k_re = L_k \ alpha_k;
  t0_re = L_t0 \ alpha_t0;
}
model {
  vector[n] log_linf ;
  vector[n] k ;
  vector[n] t0 ;
  vector[n] log_y ;
  
  log_y = log(y) ;
  
  alpha_linf ~ normal(0,1); // Scaling factors
  alpha_k ~ normal(0,1); // Scaling factors
  alpha_t0 ~ normal(0,1); // Scaling factors
  
  sigma ~ normal(0, 10);
  
  tau_linf ~ gamma(2, 2);
  tau_k ~ gamma(2, 2);
  tau_t0 ~ gamma(2, 2);
  
  Betas_log_linf ~ normal(0,10);  
  Betas_log_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  // VBGF
  log_linf = Pred * Betas_log_linf + log_linf_re[region];
  k = exp(Pred * Betas_log_k + log_k_re[region]);
  t0 = Pred * Betas_t0 + t0_re[region]; 
    
  log_y ~ normal( log_linf + log1m_exp(-k .* (age - t0)), sigma);
}



