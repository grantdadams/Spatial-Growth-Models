data {
  int<lower = 1> n;
  int<lower = 1> J;
  real weight[n];
  real length[n];
  real sex[n];
  int<lower = 0> region[n];
  matrix<lower = 0, upper = 1>[J, J] W;
  matrix<lower = 0>[J, J] D;
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
  log_a_re ~ multi_normal_prec(zeros, tau_a * (D - alpha_a * W));
  log_b_re ~ multi_normal_prec(zeros, tau_b * (D - alpha_b * W));
  
  tau_a ~ gamma(2, 2);
  tau_b ~ gamma(2, 2);
  
  mu_log_a ~ normal(0,10);
  mu_log_b ~ normal(0,10);
  
  Beta_a ~ normal(0,100);
  Beta_b ~ normal(0,100);
  
  for(i in 1:n){
    weight[i] ~ normal((a[region[i]] + Beta_a * sex[i]) * (length[i] ^ (b[region[i]] + Beta_b * sex[i])), sigma);
  }
}

