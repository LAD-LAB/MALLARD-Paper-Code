matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}


/* 
* 
* @param X array of vectors[q] of observed data (one array entry for each timepoint)
* @param F array of matricies [] of 
*
*/
real multivariate_dlm_cholesky_lpdf(vector[] X, matrix[] F, matrix[] G, 
                              matrix[] L_V, matrix[] L_W, 
                              vector m0, matrix C0, int[] TT_obs, int[] TT_rep){
  int TT = size(TT_obs); // number of timepoints
  int q = num_elements(X[1]); // number of species
  int p = rows(G[1]); // number of dimentions of univariate state space
  
  // sizes for arrays (to determine if time-varying)
  int d_F = size(F);
  int d_G = size(G);
  int d_L_V = size(L_V);
  int d_L_W = size(L_W);
  
  // system matricies for each iteration
  // matrix[p,q] F_t;
  // matrix[p,p] G_t;
  matrix[q,q] V_t;
  matrix[p,p] W_t;
  
  // internals to filter
  vector[p] m;
  matrix[p,p] C;
  vector[p] a; 
  matrix[p,p] R;
  matrix[p,q] A;
  vector[q] f;
  matrix[q,q] Q;
  vector[q] e;
  
  real log_prob_accum = 0;
  int pos=0;
  
  // initialize
  m = m0;
  C = C0;
  
  for (t in 1:TT){
    // if (d_F >= t)
    //   F_t = F[t];
    // if (d_G >= t)
    //   G_t = G[t];
    if (d_L_V >= t)
      V_t = multiply_lower_tri_self_transpose(L_V[t]);
    if (d_L_W >= t)
      W_t = multiply_lower_tri_self_transpose(L_W[t]);
    
    // prior (theta_t | D_{t-1})
    a = m;
    if (TT_rep[t] == 0){ // not a replicate (duplicate of prior sample)
      R = to_symmetric_matrix(C + W_t);
    } else {  // is a replicate (duplicate of prior sample)
      R = to_symmetric_matrix(C);
    }

    if (TT_obs[t]==1){ // Observation present
      pos = pos + 1;
      //forecast (Y_t | D_{t-1})
      f = a;
      Q = to_symmetric_matrix(R + V_t);
      
      // update (theta_t | D_t)
      A = mdivide_right_spd(R, Q);
      e = X[pos] - f;
      m = a + A*e;
      C = to_symmetric_matrix(R - quad_form_sym(Q, A'));

      // update log probability
      log_prob_accum = log_prob_accum + multi_normal_lpdf(X[pos] | f, Q);
    } else {
      // no real update due to missing observation + no change to log probability
      m = a;
      C = R;
    }
  }
  return(log_prob_accum);
}


/* 
* 
* @param X array of vectors[q] of observed data (one array entry for each timepoint)
* @param F array of matricies [] of 
*
*/
vector[] multivariate_dlm_cholesky_simsmo_rng(vector[] X, matrix[] F, matrix[] G,
                                          matrix[] L_V, matrix[] L_W, 
                                          vector m0, matrix C0, int[] TT_obs, int[] TT_rep, 
                                          int[] TT_sample_ind, 
                                          int N_timepoints_total, int N_timepoints_sample){
  int TT = N_timepoints_total; // number of timepoints
  int q = num_elements(X[1]); // number of species
  int p = rows(G[1]); // number of dimentions of univariate state space

  // sizes for arrays (to determine if time-varying)
  int d_F = size(F);
  int d_G = size(G);
  int d_L_V = size(L_V);
  int d_L_W = size(L_W);
  
  // system matricies for each iteration
  matrix[q,q] V_t;
  matrix[p,p] W_t;
  
  // internals to filter
  vector[p] m[TT+1];
  matrix[p,p] C[TT+1];
  vector[p] a[TT]; 
  matrix[p,p] R[TT];
  matrix[p,q] A;
  vector[q] f;
  matrix[q,q] Q;
  vector[q] e;
  matrix[p,p] B;
  
  // return variable
  vector[p] theta[N_timepoints_sample];
  
  int pos=0;
  int pos_s=N_timepoints_sample; // marker for position in theta for sampling
  
  // initialize
  m[1] = m0;
  C[1] = C0;
  
  for (t in 1:TT){

    if (d_L_V >= t)
      V_t = multiply_lower_tri_self_transpose(L_V[t]);
    if (d_L_W >= t)
      W_t = multiply_lower_tri_self_transpose(L_W[t]);
    
    
    // FORWARDS FILTERING
    // prior (theta_t | D_{t-1})
    a[t] = m[t];

    if (TT_rep[t] == 0){
      R[t] = to_symmetric_matrix(C[t] + W_t);
    } else {
      R[t] = to_symmetric_matrix(C[t]);
    }
    
    
    if (TT_obs[t]==1){ // Observation present
      pos = pos + 1;
      //forecast (Y_t | D_{t-1})
      f = a[t];
      Q = to_symmetric_matrix(R[t] + V_t);
      
      // update (theta_t | D_t)
      A = mdivide_right_spd(R[t], Q);
      e = X[pos] - f;
      m[t+1] = a[t] + A*e;
      C[t+1] = to_symmetric_matrix(R[t] - quad_form_sym(Q, A'));
      
    } else {
      // no real update due to missing observation + no change to log probability
      m[t+1] = a[t];
      C[t+1] = R[t];
    }
  }
  // SMOOTHING AND SAMPLING (Not using Backwards Sampling as this can crash 
  // with replicate samples)
  if (TT_sample_ind[TT] == 1){
    theta[pos_s] = multi_normal_cholesky_rng(m[TT+1], cholesky_decompose(C[TT+1]));
    pos_s = pos_s-1;
  }
  
  for (t in 1:(TT-1)){
    B = mdivide_right_spd(C[TT-t+1], R[TT-t+1]);
    m[TT-t+1] = m[TT-t+1] + B*(m[TT-t+2]-a[TT-t+1]);
    C[TT-t+1] = C[TT-t+1] + quad_form(C[TT-t+2]-R[TT-t+1], B');
    if (TT_sample_ind[TT-t]==1){ // only sample if going to use that sample
      theta[pos_s] = multi_normal_cholesky_rng(m[TT-t+1], cholesky_decompose(C[TT-t+1]));
      pos_s = pos_s-1; 
    }
  }
  return(theta);
}
