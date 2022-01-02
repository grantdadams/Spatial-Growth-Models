data {
  int<lower = 1> n;
  int<lower = 1> J;
  real weight[n];
  real length[n];
  real sex[n];
  int<lower = 0> region[n];
  matrix<lower = 0, upper = 1>[J, J] W;
  matrix<lower = 0>[J, J] D;
  matrix<lower = 0, upper = 1>[J, J] I;
}
transformed data{
  vector[J] zeros;
  zeros = rep_vector(0, J);

}
parameters {
  vector[J] log_a_re;
  vector[J] log_b_re;

  real Beta_a;
  real Beta_b;
  
  real mu_log_a ;
  real mu_log_b ;
  
  real<lower=0> sigma;
  
  real<lower = 0> tau_a;
  real<lower = 0> tau_b;

  real<lower = -1, upper = 1> alpha_a;
  real<lower = -1, upper = 1> alpha_b;
}
transformed parameters{
  vector[J] a;
  vector[J] b;
  
  a = exp(log_a_re + mu_log_a);
  b = exp(log_b_re + mu_log_b);
}
model {
    // Var-cov matrix
  matrix[J,J] C_a;
  matrix[J,J] C_b;
  
  matrix[J,J] Sigma_Inv_a;
  matrix[J,J] Sigma_Inv_b;

  matrix[n,2] WAL_Param;
  
  C_a = I - alpha_a * W';
  C_b = I - alpha_b * W';

  Sigma_Inv_a = tau_a *  C_a' * D * C_a;
  Sigma_Inv_b = tau_b * C_b' * D * C_b;
  
  for(i in 1:n){
    WAL_Param[i,1] = a[region[i]] + Beta_a * sex[i] ;
    WAL_Param[i,2] = a[region[i]] + Beta_b * sex[i] ;
  }
  
  log_a_re ~ multi_normal_prec(zeros, Sigma_Inv_a);
  log_b_re ~ multi_normal_prec(zeros, Sigma_Inv_b);
  
  tau_a ~ gamma(2, 2);
  tau_b ~ gamma(2, 2);
  
  mu_log_a ~ normal(0,10);
  mu_log_b ~ normal(0,10);
  
  Beta_a ~ normal(0,10);
  Beta_b ~ normal(0,10);

  for(i in 1:n)
    weight[i] ~ normal( WAL_Param[i,1] * (length[i] ^ WAL_Param[i,2]), sigma);
}
