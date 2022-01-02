functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param phi Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}
data {
  int<lower = 1> n;
  int<lower = 1> J;

  real weight[n];
  real length[n];

  int<lower = 1> n_pred;
  matrix[n, n_pred] Pred;
  
  int<lower = 0> region[n];
  
  matrix<lower = 0, upper = 1>[J, J] W;
  int W_n;
}
transformed data{
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[J] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[J] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
    { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(J - 1)) {
      for (j in (i + 1):J) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
    for (i in 1:J) D_sparse[i] = sum(W[i]);
  {
    vector[J] invsqrtD;  
    for (i in 1:J) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}
parameters {
  vector[J] log_a_re;
  vector[J] log_b_re;
  
  vector[n_pred] Betas_log_a;
  vector[n_pred] Betas_log_b;
  
  real<lower=0> sigma;
  
  real<lower = 0> tau_a;
  real<lower = 0> tau_b;
  
  real<lower = -1, upper = 1> phi_a;
  real<lower = -1, upper = 1> phi_b;
}
model {
  // Temporary parameters
  vector[n] a;
  vector[n] b;
  
  // Priors
  Betas_log_a ~ normal(0,10);  
  Betas_log_b ~ normal(0,10);
  
  tau_a ~ gamma(2, 2);
  tau_b ~ gamma(2, 2);
  
  sigma ~ normal(0, 10);
  
  // Likelihood
  log_a_re ~ sparse_car(tau_a, phi_a, W_sparse, D_sparse, lambda, J, W_n);
  log_b_re ~ sparse_car(tau_b, phi_b, W_sparse, D_sparse, lambda, J, W_n);
  
  a = exp(log_a_re[region] + Pred * Betas_log_a);
  b = exp(log_b_re[region] + Pred * Betas_log_b);
  
  for(i in 1:n)
  weight[i] ~ normal(a[i] * length[i]^b[i], sigma);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n)
  log_lik[i] = normal_lpdf(weight[i] | exp(Pred[i,] * Betas_log_a + log_a_re[region[i]]) * length[i]^exp(Pred[i,] * Betas_log_b + log_b_re[region[i]]), sigma);
}
