data {
  int<lower = 1> n;
  int<lower = 1> J;

  real weight[n];
  real length[n];

  int<lower = 1> n_pred;
  matrix[n, n_pred] Pred;
}
parameters {
  // WAL Coefficients
  vector[n_pred] Betas_log_a;
  vector[n_pred] Betas_log_b;
  
  real<lower=0> sigma;
}
model {
  vector[n] a ;
  vector[n] b ;
  
  // Priors
  Betas_log_a ~ normal(0,10);  
  Betas_log_b ~ normal(0,10);

  sigma ~ normal(0,10);
  
  // Likelihood
  a = exp(Pred * Betas_log_a);
  b = exp(Pred * Betas_log_b);
  
  for(i in 1:n)
  weight[i] ~ normal(a[i] * length[i]^b[i], sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(weight[i] | exp(Pred[i,] * Betas_log_a) * length[i]^exp(Pred[i,] * Betas_log_b), sigma);
}
