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
transformed data{
  real a_min;
  real a_max;
  
  a_min = min(age);
  a_max = max(age);
}
parameters {
  // Scaling parameters for non-centered distribution
  vector[J] alpha_l_min;
  vector[J] alpha_l_max;
  vector[J] alpha_k;

  // Beta coefficients for mean of MVN
  vector[n_pred] Betas_l_min;
  vector[n_pred] Betas_l_max;
  vector[n_pred] Betas_k;
  
  real<lower=0> sigma;
  
  // Hierarchical Variance
  real<lower = 0> tau_vbgf[3];
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_vbgf[3];
}
transformed parameters{
  // Lower triangle matrices of precision
  matrix[J,J] L_l_min;
  matrix[J,J] L_l_max;
  matrix[J,J] L_k;
  
  // Random effects
  vector[J] l_min_re;
  vector[J] l_max_re;
  vector[J] k_re;
  
  // Get cholesky decomp of precision
  L_l_min = cholesky_decompose(tau_vbgf[1] * (D - phi_vbgf[1] * W));
  L_l_max = cholesky_decompose(tau_vbgf[2] * (D - phi_vbgf[2] * W));
  L_k = cholesky_decompose(tau_vbgf[3] * (D - phi_vbgf[3] * W));
  
  // Get random effects
  l_min_re = L_l_min \ alpha_l_min; // equivalent to inverse(L) * alpha
  l_max_re = L_l_max \ alpha_l_max;
  k_re = L_k \ alpha_k;
}
model {
  vector[n] l_min ;
  vector[n] l_max ;
  vector[n] k ;
  vector[n] linf ;
  vector[n] t0 ;
  
  // Priors on mu VBGF params
  Betas_l_min ~ normal(0,10);
  Betas_l_max ~ normal(0,10);
  Betas_k ~ normal(0,10);
  
  // VBGF precision prior
  tau_vbgf ~ gamma(2,2);
  
  // VBGF scaling priors
  alpha_l_min ~ normal(0,1); 
  alpha_l_max ~ normal(0,1); 
  alpha_k ~ normal(0,1);
  
  // VBGF likelihood
  l_min = Pred * Betas_l_min + l_min_re[region];
  l_max = Pred * Betas_l_max + l_max_re[region];
  k = Pred * Betas_k + k_re[region];
  
  linf = l_min + (l_max - l_min) ./ (1-exp(-k*(a_max-a_min)));
  t0 = log(linf ./ (linf  - l_min)) ./ -k + a_min;
    
  log(y) ~ normal(log(l_min + (l_max - l_min) .* (1 - exp( -k .* (age - a_min))) ./ (1 - exp( -k * (a_max - a_min)))  ), sigma);
}



