data {
  int<lower = 1> n;
  int<lower = 1> J;
  int<lower = 1> n_pred;
  vector[n] y;
  vector[n] age;
  vector[n] selectivity;
  
  matrix[n, n_pred] Pred;
  int<lower = 0> region[n];
  matrix<lower = 0, upper = 1>[J, J] W;
  matrix<lower = 0>[J, J] D;
}
transformed data{
  vector[J] zeros;
  zeros = rep_vector(0, J);
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
  real<lower = 0> tau_linf;
  real<lower = 0> tau_k;
  real<lower = 0> tau_t0;
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_linf;
  real<lower = -1, upper = 1> phi_k;
  real<lower = -1, upper = 1> phi_t0;
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
  L_linf= tau_linf * (D - phi_linf * W);
  L_k = tau_k * (D - phi_k * W);
  L_t0 = tau_t0 * (D - phi_t0 * W);
  
  // Get random effects
  log_linf_re = L_linf \ alpha_linf; // equivalent to inverse(L) * alpha
  log_k_re = L_k \ alpha_k;
  t0_re = L_t0 \ alpha_t0;
}
model {
  // Model objects
  vector[n] linf ;
  vector[n] k ;
  vector[n] t0 ;
  
  alpha_linf ~ normal(0,1); // Scaling factors
  alpha_k ~ normal(0,1); // Scaling factors
  alpha_t0 ~ normal(0,1); // Scaling factors
  
  // VBGF
  linf = exp(Pred * Betas_linf + log_linf_re[region]);
  k = exp(Pred * Betas_k + log_k_re[region]);
  t0 = Pred * Betas_t0 + t0_re[region]; 
    
  y ~ normal( linf .* (1- exp( -k .* (age - t0))), sigma)T[selectivity,];
}
generated quantities{
  // Posterior predictive check
  vector[n] new_y;
  vector[n] obs_resid;
  vector[n] new_resid;
  real fit_obs;
  real fit_new;
  real bp;
  
  // Generate fake data
  for(i in 1:n){
    new_y[i] = normal_rng(exp(Pred[i,] * Betas_linf + log_linf_re[region[i]]) .* (1- exp( -exp(Pred[i,] * Betas_k + log_k_re[region[i]]) .* (age[i] - (Pred[i,] * Betas_t0 + t0_re[region[i]])))), sigma)T[selectivity[i],];

  // Residuals
  obs_resid[i] = (y[i] - exp(Pred[i,] * Betas_linf + log_linf_re[region[i]]) .* (1- exp( -exp(Pred[i,] * Betas_k + log_k_re[region[i]]) .* (age[i] - (Pred[i,] * Betas_t0 + t0_re[region[i]])))))^2;
  new_resid[i] =  (new_y[i] - exp(Pred[i,] * Betas_linf + log_linf_re[region[i]]) .* (1- exp( -exp(Pred[i,] * Betas_k + log_k_re[region[i]]) .* (age[i] - (Pred[i,] * Betas_t0 + t0_re[region[i]])))))^2;
    
  }

  fit_obs = sum(obs_resid[]);
  fit_new = sum(new_resid[]);
  
  bp = step(fit_obs - fit_new);
}


