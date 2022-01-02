data {
  int<lower = 1> n;
  int<lower = 1> J;

  vector[n] weight;
  vector[n] length;

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
  vector[n] log_a ;
  vector[n] b ;
  vector[n] log_weight;
  vector[n] log_length;
  
  log_weight = log(weight);
  log_length = log(length);
  
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
  log_a = (Pred * Betas_log_a + log_a_re[region]);
  b = exp(Pred * Betas_log_b + log_b_re[region]);
  
  
  log_weight ~ normal(log_a +  log_length .* b, sigma);
}

