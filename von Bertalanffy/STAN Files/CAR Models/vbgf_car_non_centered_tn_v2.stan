data {
  int<lower = 1> n; // Sample size
  int<lower = 1> J; // Number of possible regions 11. NOTE: There is no data from region 5
  int<lower = 1> n_pred; // Number of predictors for mean of hierarchical parameter distributions (3): 1, sex, lat
  vector[n] y; // length
  vector[n] age;
  real selectivity[n]; // Minimum size from region/sex
  matrix[n, n_pred] Pred; // Predictor matrix: 1, sex, lat
  int<lower = 0> region[n]; // Region ID: 1:11. NOTE: There is no data from region 5
  matrix<lower = 0, upper = 1>[J, J] W; // 11x11 Neighbor matrix
  matrix<lower = 0>[J, J] D; // 11x11 Neighbor weight matrix
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
  
  // Hyperprior parameter precision
  real<lower = 0> tau_vbgf[3];
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_vbgf[3];
}
transformed parameters{
  // Lower triangle matrices of precision
  matrix[J,J] L_linf;
  matrix[J,J] L_k;
  matrix[J,J] L_t0;
  
  // Random effects
  vector[J] log_linf_re;
  vector[J] log_k_re;
  vector[J] t0_re;
  
  // Get cholesky decomp of precision
  L_linf= cholesky_decompose(tau_vbgf[1] * (D - phi_vbgf[1] * W));
  L_k = cholesky_decompose(tau_vbgf[2] * (D - phi_vbgf[2] * W));
  L_t0 = cholesky_decompose(tau_vbgf[3] * (D - phi_vbgf[3] * W));
  
  // Get random effects
  log_linf_re = inverse(L_linf) * alpha_linf; // Tried L_linf \ alpha_linf, but it does not work with stan_model compiler
  log_k_re = inverse(L_k) * alpha_k;
  t0_re = inverse(L_t0) * alpha_t0;
}
  
  vector[n] linf ;
  vector[n] k ;
  vector[n] t0 ;
  vector[n] log_mu;
  
  // Priors on mu VBGF params
  Betas_linf ~ normal(0,10);
  Betas_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  // VBGF precision prior
  tau_vbgf ~ gamma(2,2);
  sigma ~ inv_gamma(2,2);
  
  // VBGF scaling priors
  alpha_linf ~ normal(0,1); 
  alpha_k ~ normal(0,1); 
  alpha_t0 ~ normal(0,1);
  
  // VBGF Parameters
  linf = exp(log_linf_re[region] + Pred * Betas_linf);
  k = exp(log_k_re[region] + Pred * Betas_k);
  t0 = t0_re[region] + Pred * Betas_t0;
  log_mu = log(linf .* (1- exp( -k .* (age - t0))));
  
  // Likelihood vectorized
  log(y) ~ normal( log_mu, sigma);
  
  target += -normal_lccdf(log(selectivity) | log_mu, sigma);
}


