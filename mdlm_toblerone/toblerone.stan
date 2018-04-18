functions {
  #include mdlm_utils.stan
  
  // ILR transform
  vector[,] ILR_INV(matrix[] X, matrix contrast_matrix_t) {
    int N_timepoints = size(X); 
    int q = rows(contrast_matrix_t);
    int r = rows(X[1]);
    matrix[r, q+1] y;
    vector[q+1] y_a[N_timepoints, r];
    
    for (i in 1:N_timepoints){
      y = exp(X[i]*contrast_matrix_t);
      for (j in 1:r){
        y_a[i,j] = (y[j]/sum(y[j]))';
      }
    }
    return(y_a);
  }
  
  // ILR transform that accepts vectors
  vector[] ILR_INV_V(vector[] x, matrix contrast_matrix, int N_samples, int K) {
    vector[K+1] y[N_samples];

    for (i in 1:N_samples){
      y[i] = exp(contrast_matrix*x[i]);
      y[i] = y[i]/sum(y[i]);
    }
    return(y);
  }
}
data {
  int<lower=1> N_timepoints_total; // number of samples in each time-series (missing + non_missing)
  int<lower=1> N_timepoints_observed; // number of samples acctually observed (non_missing)
  int<lower=1> N_timepoints_sample; // number of time-points to sample from (non_missing with selected missing)
  
  int<lower=1> r; // number of time series i
  int<lower=2> N_species; // number of species (e.g., OTUs in dataset)
  int<lower=1> p; // Dimention of each series' state-space
  
  int<lower=0> TT_obs_ind[N_timepoints_total, r]; // indicator of time points of observations
  int<lower=0> TT_rep_ind[N_timepoints_total, r]; // indicator of timepoints that are replicates
  int<lower=0> TT_sample_ind[N_timepoints_total]; // indicator of timepoints to sample from

  
  int<lower=0> Y[N_timepoints_observed, r, N_species]; // counts
  int<lower=0, upper=1> Y_obs[N_timepoints_observed, r];
  matrix[N_species,N_species-1] contrast_matrix;
  
  // DLM matricies
  int d_F;
  int d_G;
  matrix[p,N_species-1] FF[d_F];
  matrix[p,p] GG[d_G];

  // initialization parameters and Priors
  vector[p] m0[r];
  cov_matrix[p] C0;
  real W_scale_mean_prior[p];
  real V_scale_mean_prior[N_species-1];
  real<lower=0> W_scale_var_prior[p];
  real<lower=0> V_scale_var_prior[N_species-1];
  real<lower=0.00001> W_lkj_prior;
  real<lower=0.00001> V_lkj_prior;
}
transformed data {
  int<lower=1> q = N_species-1; // number of ILR components
  matrix[q, q+1] contrast_matrix_t = contrast_matrix';
  int<lower=0> TT_obs_combined_ind[N_timepoints_total]; // For Indexing Y
  
  for(i in 1:N_timepoints_total)
    TT_obs_combined_ind[i] = sum(TT_obs_ind[i,])>0;
}
parameters {
  vector[q] eta[N_timepoints_observed, r];
  vector<lower=0>[p] W_scale[1];
  vector<lower=0>[q] V_scale[1];
  cholesky_factor_corr[p] W_corr[1];
  cholesky_factor_corr[q] V_corr[1];
}
transformed parameters {
  simplex[N_species] pi[N_timepoints_observed, r];
  cholesky_factor_cov[p] L_W[1];
  cholesky_factor_cov[q] L_V[1];
  
  L_W[1] = diag_pre_multiply(W_scale[1], W_corr[1]);
  L_V[1] = diag_pre_multiply(V_scale[1], V_corr[1]);
  for (i in 1:r)
    pi[,i] = ILR_INV_V(eta[,i], contrast_matrix, N_timepoints_observed, q);
}
model {
  W_scale[1] ~ lognormal(W_scale_mean_prior, W_scale_var_prior);
  V_scale[1] ~ lognormal(V_scale_mean_prior, V_scale_var_prior);
  W_corr[1] ~ lkj_corr_cholesky(W_lkj_prior);
  V_corr[1] ~ lkj_corr_cholesky(V_lkj_prior);
  for (i in 1:r) // 4 DLMs are "independent"
    target += multivariate_dlm_cholesky_lpdf(eta[,i] | FF, GG, L_V, L_W, m0[i], C0, 
                                                       TT_obs_combined_ind,
                                                       // TT_obs_ind[,i], 
                                                       TT_rep_ind[,i]);
  for (i in 1:N_timepoints_observed){
    for (j in 1:r){
      if (Y_obs[i,j]==1){
        Y[i, j] ~ multinomial(pi[i, j]);        
      }
    }
  }
}
generated quantities {
  vector[p] theta[N_timepoints_sample, r];
  cov_matrix[q] V;
  cov_matrix[p] W;
  for (i in 1:r) // 4 DLMs are "independent"
    theta[,i] = multivariate_dlm_cholesky_simsmo_rng(eta[,i], FF, GG, L_V, L_W, m0[i], C0, 
                                                     TT_obs_combined_ind,
                                                     // TT_obs_ind[,i],
                                                     TT_rep_ind[,i],
                                                     TT_sample_ind,
                                                     N_timepoints_total, 
                                                     N_timepoints_sample);
  V = multiply_lower_tri_self_transpose(L_V[1]);
  W = multiply_lower_tri_self_transpose(L_W[1]);
}
