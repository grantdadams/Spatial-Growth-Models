data {
  int<lower = 1> n;
  int<lower = 1> J;

  real weight[n];
  real length[n];

  int<lower = 1> n_pred;
  matrix[n, n_pred] Pred;
  
  int<lower = 0> region[n];
  
  matrix<lower = 0, upper = 1>[J, J] W;
  matrix<lower = 0>[J, J] D;
}
parameters {
  // Scaling parameters for non-centered distribution
  vector[J] alpha_a;
  vector[J] alpha_b;

  // Beta coefficients for mean of MVN
  vector[n_pred] Betas_log_a;
  vector[n_pred] Betas_log_b;
  
  real<lower=0> sigma;
  
  // Hierarchical Variance
  real<lower = 0> tau_a;
  real<lower = 0> tau_b;
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_a;
  real<lower = -1, upper = 1> phi_b;
}
transformed parameters{
  // Lower triangle matrices of precision
  matrix[J,J] L_a;
  matrix[J,J] L_b;
  
  // Random effects
  vector[J] log_a_re;
  vector[J] log_b_re;
  
  // Get cholesky decomp of precision
  L_a = tau_a * (D - phi_a * W);
  L_b = tau_b * (D - phi_b * W);
  
  // Get random effects
  log_a_re = L_a \ alpha_a; // equivalent to inverse(L) * alpha
  log_b_re = L_b \ alpha_b;
}
model {
  // Temporary parameters
  vector[n] a ;
  vector[n] b ;
  
  // Priors
  alpha_a ~ normal(0,1); // Scaling factors
  alpha_b ~ normal(0,1); // Scaling factors
  
  phi_a ~ normal(0, 10);
  phi_b ~ normal(0, 10);
  
  sigma ~ normal(0, 10);
  
  tau_a ~ gamma(2, 2);
  tau_b ~ gamma(2, 2);
  
  Betas_log_a ~ normal(0,10);  
  Betas_log_b ~ normal(0,10);
  
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


