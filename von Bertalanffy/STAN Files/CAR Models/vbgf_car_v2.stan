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
  vector[J] log_linf_re;
  vector[J] log_k_re;
  vector[J] t0_re;

  vector[n_pred] Betas_linf;
  vector[n_pred] Betas_k;
  vector[n_pred] Betas_t0;
  
  real<lower=0> sigma;
  
  // Hierarchical Variance
  real<lower = 0> tau_vbgf[3];
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_vbgf[3];
}
model {
  log_linf_re ~ multi_normal_prec(zeros, tau_vbgf[1] * (D - phi_vbgf[1] * W));
  log_k_re ~ multi_normal_prec(zeros, tau_vbgf[2] * (D - phi_vbgf[2] * W));
  t0_re ~ multi_normal_prec(zeros, tau_vbgf[3] * (D - phi_vbgf[3] * W));
  
  tau_vbgf ~ gamma(2, 2);

  log(y) ~ normal(exp(log_linf_re[region] + Pred * Betas_linf) .* (1- exp( -exp(log_k_re[region] + Pred * Betas_k) .* (age - (t0_re[region] + Pred * Betas_t0)))), sigma);

}

