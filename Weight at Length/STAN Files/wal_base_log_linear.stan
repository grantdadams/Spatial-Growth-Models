data {
  int<lower = 1> n;
  int<lower = 1> J;
  // int<lower = 1> n_pred;
  real log_weight[n];
  real log_length[n];
  // real selectivity[n];
  // matrix[n, n_pred] Pred;
  // int<lower = 0> region[n];
}
parameters {
  // VBGF Coefficients
  real log_a;
  real log_b;
  
  real<lower=0> sigma;
}
transformed parameters{
  real a;
  real b;
  
  a = exp(log_a);
  b = exp(log_b);
}
model {
  // Priors
  log_a ~ normal(0,10);  
  log_b ~ normal(0,10);

  sigma ~ normal(0,10);
  
  // Likelihood
  for(i in 1:n)
  log_weight[i] ~ normal( log_a + b * log_length[i], sigma);
}
