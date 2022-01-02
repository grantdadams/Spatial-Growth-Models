functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
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
  int<lower = 1> n_pred;
  real y[n];
  vector[n] age;
  real selectivity[n];
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
  vector[J] log_linf_re;
  vector[J] log_k_re;
  vector[J] log_t0_re;

  vector[n_pred] Betas_linf_loc;
  vector[n_pred] Betas_k_loc;
  vector[n_pred] Betas_t0_loc;
  
  real<lower=0> sigma;
  
  // Hierarchical Variance
  real<lower = 0> tau_vbgf[3];
  
  // Spatial-correlation factors
  real<lower = -1, upper = 1> phi_vbgf[3];
}
model {
  // Temporary parameters
  vector[n] linf;
  vector[n] k;
  vector[n] t0;
  vector[n] log_mu;
  
  log_linf_re ~ sparse_car(tau_vbgf[1], phi_vbgf[1], W_sparse, D_sparse, lambda, J, W_n);
  log_k_re ~ sparse_car(tau_vbgf[2], phi_vbgf[2], W_sparse, D_sparse, lambda, J, W_n);
  log_t0_re ~ sparse_car(tau_vbgf[3], phi_vbgf[3], W_sparse, D_sparse, lambda, J, W_n);

  tau_vbgf ~ gamma(2, 2);
  
  // VBGF
  linf = exp(log_linf_re[region] + Pred * Betas_linf_loc);
  k = exp(log_k_re[region] + Pred * Betas_k_loc);
  t0 = log_t0_re[region] + Pred * Betas_t0_loc;
  log_mu = log(linf .* (1- exp( -k .* (age - t0))));
  
  // Likelihood vectorized
  log(y) ~ normal( log_mu, sigma);
  
  target += -normal_lccdf(log(selectivity) | log_mu, sigma);
}

