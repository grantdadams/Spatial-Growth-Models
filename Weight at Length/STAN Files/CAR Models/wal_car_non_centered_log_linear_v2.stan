data {
  int<lower = 1> n;
  int<lower = 1> J;

  vector[n] log_weight;
  vector[n] log_length;

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
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
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
  
  // Precision
  real tau_a;
  real tau_b;
    
  tau_a = 1/(sigma_a^2);
  tau_b = 1/(sigma_b^2);
  
  // Get cholesky decomp of precision
  L_a = tau_a * (D - phi_a * W);
  L_b = tau_b * (D - phi_b * W);
  
  // Get random effects
  log_a_re = L_a \ alpha_a; // equivalent to inverse(L) * alpha
  log_b_re = L_b \ alpha_b;
}
model {
  // Temporary parameters
  vector[n] log_a ;
  vector[n] b ;
  
  // Priors
  alpha_a ~ normal(0,1); // Scaling factors
  alpha_b ~ normal(0,1); // Scaling factors
  
  phi_a ~ normal(0, 10);
  phi_b ~ normal(0, 10);
  
  sigma ~ cauchy(0, 5);
  sigma_a ~ cauchy(0, 5);
  sigma_b ~ cauchy(0, 5);
  
  Betas_log_a ~ normal(0,10);  
  Betas_log_b ~ normal(0,10);
  
  // Likelihood
  log_a = (Pred * Betas_log_a + log_a_re[region]);
  b = exp(Pred * Betas_log_b + log_b_re[region]);
  
  log_weight ~ normal(log_a + b .* log_length, sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(log_weight[i] | (Pred[i,] * Betas_log_a + log_a_re[region[i]]) + log_length[i] * exp(Pred[i,] * Betas_log_b + log_b_re[region[i]]), sigma);
}

