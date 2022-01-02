data {
  // ------------------------------------------------------------------------------------ //
  // 1. Assign data to stan objects                                                       //
  // ------------------------------------------------------------------------------------ //
  int<lower = 1>      n_i;            // Sample size
  int<lower = 1>      n_r;            // Number of regions
  int<lower = 1>      n_pred;         // Number of predictors
  vector[n_i]         log_length_i;   // Vector of log fork length of fish i
  vector[n_i]         log_weight_i;   // Vector of log weight (g) of fish i
  matrix[n_i, n_pred] design_mat;     // n_i * n_pred design matrix of linear predictors
  int<lower = 0>      r_i[n_i];       // Integer vector of region of fish i
  matrix<lower = 0,   upper = 1>[n_r, n_r] W; // Weights matrix where W_ij = 1 when regions i and j are neighbors
  matrix<lower = 0>   [n_r, n_r] D;       // Matrix where diagonal is number of neighbors of region i; diag(W * 1) (Eq. 10)
}
parameters {
  // ------------------------------------------------------------------------------------ //
  // 2. Specify model parameters                                                          //
  // ------------------------------------------------------------------------------------ //
  // 2.1. Regression coefficients for mean of weight-at-length parameters
  vector[n_pred] B_log_a;
  vector[n_pred] B_log_b;
  
  // 2.2. Random error (SD)
  real<lower=0> sigma;
  
  // 2.3. Hierarchical weight-at-length parameter error (SD)
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;

  // 2.4. Scaling parameters for non-centered distribution
  vector[n_r] alpha_a;
  vector[n_r] alpha_b;
  
  // 2.5. Spatial-correlation factors
  real<lower = -1, upper = 1> phi_a;
  real<lower = -1, upper = 1> phi_b;
}
transformed parameters{
  // ------------------------------------------------------------------------------------ //
  // 3. Specify derived model parameters to be saved                                      //
  // ------------------------------------------------------------------------------------ //
  // 3.1. Inverse of CAR model matrix
  matrix[n_r, n_r] L_a;
  matrix[n_r, n_r] L_b;
 
  // 3.2. Weight-at-length parameter-specific random effects from CAR model
  vector[n_r] log_a_re;
  vector[n_r] log_b_re;

  // 3.3. Precision of hierearchical error
  real tau_a;
  real tau_b;
  
  // 3.4. Get precision
  tau_a = 1/(sigma_a^2);
  tau_b = 1/(sigma_b^2);

  // 3.5. Get inverse of CAR model matrix
  L_a  = tau_a * (D - phi_a * W);
  L_b  = tau_b * (D - phi_b * W);

  // 3.6. Get random effects
  log_a_re = L_a \ alpha_a; // Equivalent to inverse(L) * alpha
  log_b_re = L_b \ alpha_b;

}
model {
  // ------------------------------------------------------------------------------------ //
  // 4. Specify likelihood and priors                                                     //
  // ------------------------------------------------------------------------------------ //
  // 4.1. Temporary model objects to save parameter vectors
  vector[n_i] log_a_i ;
  vector[n_i] b_i ;

  // 4.2. Priors
  // 4.2.1. -- Regression coefficients for mean of weight-at-length parameters
  B_log_a     ~ normal(0,10);  
  B_log_b     ~ normal(0,10);
  
  // 4.2.2 -- Error components
  sigma       ~ cauchy(0, 5);
  sigma_a     ~ cauchy(0, 5);
  sigma_b     ~ cauchy(0, 5);

  // 4.2.3 -- Scaling factors for non-centered parameterization
  alpha_a     ~ normal(0,1);
  alpha_b     ~ normal(0,1);
  
  // 4.2.4. -- Spatial correlation parameters 
  phi_a       ~ normal(0,10);
  phi_b       ~ normal(0,10);

  // 4.3. Model specification
  // 4.3.1. -- Weight-at-length parameters
  log_a_i = design_mat * B_log_a + log_a_re[r_i];
  b_i     = exp(design_mat * B_log_b + log_b_re[r_i]);
  
  // 4.3.2. -- Model likelihood
  log_weight_i ~ normal(log_a_i + log_length_i .* b_i, sigma);
}
generated quantities {
  vector[n_i] log_lik;
  for(i in 1:n_i)
  log_lik[i] = normal_lpdf(log_weight_i[i] | (design_mat[i,] * B_log_a + log_a_re[r_i[i]]) + log_length_i[i] * exp(design_mat[i,] * B_log_b + log_b_re[r_i[i]]), sigma);
}


