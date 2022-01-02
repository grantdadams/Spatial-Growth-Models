data {
  int<lower = 1> n;
  int<lower = 1> J;
  // int<lower = 1> n_pred;
  vector[n] log_y;
  vector[n] age;
  // real selectivity[n];
  // matrix[n, n_pred] Pred;
  // int<lower = 0> region[n];
}
parameters {
  // VBGF Coefficients
  real log_linf;
  real log_k;
  real t0;
  
  real<lower=0> sigma;
}
transformed parameters{
  real linf;
  real k;
  
  linf = exp(log_linf);
  k = exp(log_k);
}
model {
  log_linf ~ normal(0,10);  
  log_k ~ normal(0,10);
  t0 ~ normal(0,10);
  
  sigma ~ normal(0,10);
  
  log_y ~ normal( log_linf + log1m_exp(-k * (age - t0)), sigma);
}
