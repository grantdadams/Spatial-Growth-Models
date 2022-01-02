data {
  int<lower = 1> n;
  int<lower = 1> J;

  real weight[n];
  real length[n];

  int<lower = 1> n_pred;
  matrix[n, n_pred] Pred;
  
  int<lower = 0> region[n];
}
parameters {
  // WAL Coefficients
  vector[n_pred] Betas_log_a;
  vector[n_pred] Betas_log_b;
  
  real<lower=0> sigma;
  
  // Scaling parameters for non-centered distribution
  vector[J] alpha_a;
  vector[J] alpha_b;
  
  // Hierarchical Variance
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
}
transformed parameters{
  // Random effects
  vector[J] log_a_re;
  vector[J] log_b_re;
  
  // Get random effects
  log_a_re = sigma_a * alpha_a; // equivalent to ~ N(0, sigma)
  log_b_re = sigma_b * alpha_b;
}
model {
  vector[n] a ;
  vector[n] b ;
  
  // Priors
  Betas_log_a ~ normal(0,10);  
  Betas_log_b ~ normal(0,10);

  sigma ~ normal(0,10);
  
  // Parameter variance
  sigma_a ~ normal(0,10);
  sigma_b ~ normal(0,10);
  
  // Scaling factors
  alpha_a ~ normal(0,1);
  alpha_b ~ normal(0,1);
  
  // Likelihood
  a = exp(Pred * Betas_log_a + log_a_re[region]);
  b = exp(Pred * Betas_log_b + log_b_re[region]);
  
  for(i in 1:n)
  weight[i] ~ normal(a[i] * length[i]^b[i], sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(weight[i] | exp(Pred[i,] * Betas_log_a + log_a_re[region[i]]) * length[i]^exp(Pred[i,] * Betas_log_b + log_b_re[region[i]]), sigma);
}
