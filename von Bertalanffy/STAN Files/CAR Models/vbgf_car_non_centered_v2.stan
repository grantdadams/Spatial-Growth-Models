data {
  int<lower = 1> n; // Sample size
  int<lower = 1> J; // Number of possible regions 11
  int<lower = 1> n_pred; // Number of predictors (3): 1, sex, lat
  vector[n] y; // length
  vector[n] age;
  real selectivity[n];
  
  matrix[n, n_pred] Pred; // Predictor matrix: 1, sex, lat
  int<lower = 0> region[n]; // Region ID w/o Georgia: 1:10
  matrix<lower = 0, upper = 1>[J, J] W; // 11x11 Neighbor matrix where row/column 11 is georgia
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
  
  // Hierarchical Variance
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
  log_linf_re = L_linf \ alpha_linf; // equivalent to inverse(L) * alpha
  log_k_re = L_k \ alpha_k;
  t0_re = L_t0 \ alpha_t0;
}
model {
  vector[n] linf ;
  vector[n] k ;
  vector[n] t0 ;
  
  // Priors on mu VBGF params
  Betas_linf ~ normal(0,10);
  Betas_k ~ normal(0,10);
  Betas_t0 ~ normal(0,10);
  
  // VBGF precision prior
  tau_vbgf ~ gamma(2,2);
  
  // VBGF scaling priors
  alpha_linf ~ normal(0,1); 
  alpha_k ~ normal(0,1); 
  alpha_t0 ~ normal(0,1);
  
  // VBGF likelihood
  linf = exp(Pred * Betas_linf + log_linf_re[region]);
  k = exp(Pred * Betas_k + log_k_re[region]);
  t0 = Pred * Betas_t0 + t0_re[region]; 
    
  log(y) ~ normal(log(linf .* (1 - exp( -k .* (age - t0)))), sigma);
}



