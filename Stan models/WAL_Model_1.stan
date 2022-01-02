data {
  // ------------------------------------------------------------------------------------ //
  // 1. Assign data to stan objects                                                       //
  // ------------------------------------------------------------------------------------ //
  int<lower = 1>      n_i;            // Sample size
  int<lower = 1>      n_pred;         // Number of predictors
  vector[n_i]         log_length_i;   // Vector of log fork length of fish i
  vector[n_i]         log_weight_i;   // Vector of log weight (g) of fish i
  matrix[n_i, n_pred] design_mat;     // n_i * n_pred design matrix of linear predictors
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
}
model {
  // ------------------------------------------------------------------------------------ //
  // 3. Specify likelihood and priors                                                     //
  // ------------------------------------------------------------------------------------ //
  // 3.1. Temporary model objects to save parameter vectors
  vector[n_i] log_a_i ;
  vector[n_i] b_i ;
  
  // 3.2. Priors
  // 3.2.1. -- Regression coefficients for mean of weight-at-length parameters
  B_log_a ~ normal(0,10);  
  B_log_b ~ normal(0,10);

  // 3.2.2 -- Error components
  sigma ~ cauchy(0,5);
  
  // 3.3. Model specification
  // 3.3.1. -- Weight-at-length parameters
  log_a_i = (design_mat * B_log_a);
  b_i     = exp(design_mat * B_log_b);
  
  // 3.3.2. -- Model likelihood
  log_weight_i ~ normal(log_a_i + log_length_i .* b_i, sigma);
}
generated quantities {
  vector[n_i] log_lik;
  for(i in 1:n_i)
  log_lik[i] = normal_lpdf(log_weight_i[i] | (design_mat[i,] * B_log_a) + log_length_i[i] * exp(design_mat[i,] * B_log_b), sigma);
}
