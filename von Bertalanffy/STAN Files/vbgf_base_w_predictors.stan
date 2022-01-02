data {
  int<lower = 1> n;
  int<lower = 1> J;
  int<lower = 1> n_pred;
  vector[n] log_y;
  vector[n] age;
  // real selectivity[n];
  matrix[n, n_pred] Pred;
  int<lower = 0> region[n];
}
parameters {
  // Beta coefficients for mean of MVN
  vector[n_pred] Betas_log_linf;
  vector[n_pred] Betas_log_k;
  vector[n_pred] Betas_t0;
  
  real<lower=0> sigma;
}

model {
  // Temporary model bjects
  vector[n] log_linf ;
  vector[n] k ;
  vector[n] t0 ;
  
  // Priors
  Betas_log_linf ~ normal(0,10);  
  Betas_log_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  sigma ~ cauchy(0,5);
  
  // Likelihood
  log_linf = Pred * Betas_log_linf;
  k = exp(Pred * Betas_log_k);
  t0 = Pred * Betas_t0; 
  
  log_y ~ normal( log_linf + log1m_exp(-k .* (age - t0)), sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(log_y[i] | (Pred[i,] * Betas_log_linf) + log1m_exp(-exp(Pred[i,] * Betas_log_k) .* (age[i] - (Pred[i,] * Betas_t0))), sigma);
}
