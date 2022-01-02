data {
  int<lower = 1> n;
  int<lower = 1> J;
  int<lower = 1> n_pred;
  vector[n] y;
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
transformed parameters{
  vector[n_pred] Betas_linf;
  vector[n_pred] Betas_k;
  
  Betas_linf = exp(Betas_log_linf);
  Betas_k = exp(Betas_log_k);
}
model {
  vector[n] log_linf ;
  vector[n] k ;
  vector[n] t0 ;
  vector[n] log_y ;
  
  Betas_log_linf ~ normal(0,10);  
  Betas_log_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  // VBGF
  log_linf = Pred * Betas_log_linf;
  k = exp(Pred * Betas_log_k);
  t0 = Pred * Betas_t0; 
  
  sigma ~ normal(0,10);
  
  log_y = log(y) ;
  
  log_y ~ normal( log_linf + log1m_exp(-k .* (age - t0)), sigma);
}

