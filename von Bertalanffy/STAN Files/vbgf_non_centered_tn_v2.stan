data {
  int<lower = 1> n; // Sample size
  int<lower = 1> J; // Number of regions with data (10)
  int<lower = 1> n_pred; // Number of predictors for mean of hierarchical parameter distributions (3): 1, sex, lat
  vector[n] y; // length
  vector[n] age;
  vector[n] selectivity; // Minimum size from region/sex
  matrix[n, n_pred] Pred; // Predictor matrix: 1, sex, lat
  int<lower = 0> region[n]; // Region ID: 1:10
}
parameters {
  // Scaling parameters for non-centered distribution
  vector[J] alpha_linf;
  vector[J] alpha_k;
  vector[J] alpha_t0;

  // Beta coefficients for mean of MVN
  vector[n_pred] Betas_linf;
  vector[n_pred] Betas_k;
  vector[n_pred] Betas_t0;
  
  real<lower=0> sigma;
  
  // Hyperprior parameter Variance
  real<lower = 0> sigma_vbgf[3];
}
transformed parameters{
  // Random effects
  vector[J] log_linf_re;
  vector[J] log_k_re;
  vector[J] t0_re;
  
  // Get random effects
  log_linf_re = sigma_linf * alpha_linf; // equivalent to ~ N(0, sigma)
  log_k_re = sigma_k * alpha_k;
  t0_re = sigma_t0 * alpha_t0;
}
model {
  vector[n] linf ;
  vector[n] k ;
  vector[n] t0 ;
  vector[n] log_mu ;
  
  // VBGF variance prior
  sigma ~ inv_gamma(2,2);
  sigma_vbgf ~ inv_gamma(2,2);
  
  // Priors on mu VBGF params
  Betas_linf ~ normal(0,10);  
  Betas_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  // VBGF scaling priors
  alpha_linf ~ normal(0,1); // Scaling factors
  alpha_k ~ normal(0,1); // Scaling factors
  alpha_t0 ~ normal(0,1); // Scaling factors
  
  // VBGF Parameters
  linf = exp(Pred * Betas_linf + log_linf_re[region]);
  k = exp(Pred * Betas_k + log_k_re[region]);
  t0 = Pred * Betas_t0 + t0_re[region]; 
  
  // Likelihood vectorized
  log_mu = log(linf .* (1- exp( -k .* (age - t0))));
    
  log(y) ~ normal(log_mu, sigma);
  
  // Add likelihood component of truncation
  target +=  -normal_lccdf(log(selectivity) | log_mu, sigma); // Normal complementary cumulative dist. equivalent to log(1 - cdf(selectivity, mu, sigma))
}



