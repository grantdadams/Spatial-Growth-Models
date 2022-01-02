data {
  int<lower = 1> n;
  int<lower = 1> J;

  vector[n] log_weight;
  vector[n] log_length;

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
  vector[n] log_a ;
  vector[n] b ;
  
  // Priors
  Betas_log_a ~ normal(0,10);  
  Betas_log_b ~ normal(0,10);

  sigma ~ cauchy(0,5);
  
  // Likelihood
  log_a = (Pred * Betas_log_a);
  b = exp(Pred * Betas_log_b);
  

  log_weight ~ normal(log_a + log_length .* b, sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(log_weight[i] | (Pred[i,] * Betas_log_a) + log_length[i] * exp(Pred[i,] * Betas_log_b), sigma);
}
