data {
  
  // Poll data
  int N;
  int id_cand;
  int tot_eff[N];
  real<lower=0,upper=1> vshare_raw[N];
  real sum_share[N];
  int rounding_ind[N];
  int r_1[N];
  int N_1;
  int r_2[N];
  int N_2;
  int r_3[N];
  int N_3;
  int r_4[N];
  int N_4;
  int r_5[N];
  int N_5;
  int id_poll[N];
  int P;
  int id_date[N];
  int id_month[N];
  int id_house[N];
  int F;
  matrix[N,3] X;
  
  // EZ and CT adjustement
  int isn_z[N];
  int is_t[N];
  
  // Splines
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int num_basis;
  int D;
  matrix[num_basis,D] S;
  
}
parameters {
  
  // Splines
  real alpha0;
  real alpha_raw[num_basis];
  real<lower=0> tau_alpha;
  
  // Covariates
  real mu[P];
  real lambda[F];
  real beta[3];
  real nu[2];
  real<lower=0> tau_mu;
  real<lower=0> tau_lambda;
  
  // EZ and CT adjustment
  real gamma_z_tilde;
  real gamma_t_tilde;
  
  // Rounding error
  real<lower=-0.005,upper=0.005> epsilon1[N_1];
  real<lower=-0.01,upper=0.01> epsilon2[N_2];
  real<lower=0.001,upper=0.015> epsilon3[N_3];
  real<lower=0.001,upper=0.01> epsilon4[N_4];
  real<lower=0.001,upper=0.005> epsilon5[N_5];
  
}
transformed parameters {
  
  real alpha[num_basis];
  real gamma_z;
  real gamma_t;
  real<lower=0,upper=1> theta[N];
  real<lower=-.015,upper=.015> epsilon[N];
  real<lower=0,upper=1> vshare[N];
  
  // Spline coefficients, specified as a random walk
  // to avoid overfit
  alpha[1] = alpha_raw[1];
  for (i in 2:num_basis)
    alpha[i] = alpha[i-1] + alpha_raw[i] * tau_alpha;
  
  // EZ and CT adjustment
  if (id_cand < 12) {
    gamma_z = gamma_z_tilde;
    gamma_t = gamma_t_tilde;
  } else {
    gamma_z = 0;
    gamma_t = gamma_t_tilde;
  }
  
  
  for (i in 1:N) {
    
    // Theta
    if (id_month[i] == 1) {
      theta[i] = inv_logit(alpha0 * id_date[i] + to_row_vector(alpha) * S[,id_date[i]] + // Spline
                                  tau_mu * mu[id_poll[i]] + // Poll effect
                                  tau_lambda * lambda[id_house[i]] + // House effect
                                  X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3] +
                                  gamma_z * isn_z[i]); // EZ adjustment
    } else if (id_month[i] < 4 || id_month[i] > 5) {
      theta[i] = inv_logit(alpha0 * id_date[i] + to_row_vector(alpha) * S[,id_date[i]] + // Spline
                                  tau_lambda * lambda[id_house[i]] + // House effect
                                  X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3]);
    } else {
      theta[i] = inv_logit(alpha0 * id_date[i] + to_row_vector(alpha) * S[,id_date[i]] + // Spline
                                  tau_mu * mu[id_poll[i]] + // Poll effect
                                  tau_lambda * lambda[id_house[i]] + // House effect
                                  X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3] +
                                  gamma_t * is_t[i]); // CT adjustment
    }
    
    
    // Rounding error
    if (rounding_ind[i] == 1) {
      epsilon[i] = epsilon1[r_1[i]];
    } else if (rounding_ind[i] == 2) {
      epsilon[i] = epsilon2[r_2[i]];
    } else if (rounding_ind[i] == 3) {
      epsilon[i] = epsilon3[r_3[i]];
    } else if (rounding_ind[i] == 4) {
      epsilon[i] = epsilon4[r_4[i]];
    } else {
      epsilon[i] = epsilon5[r_5[i]];
    }
    vshare[i] = (vshare_raw[i] + epsilon[i]) / sum_share[i];
    
    
  }
  
}
model {
  
  
  // Priors
  
    // Splines
    alpha0 ~ normal(0, 1);
    alpha_raw ~ normal(0, 2);
    tau_alpha ~ student_t(3, 0, 1);
    
    // Poll effect
    mu ~ normal(0, 1);
    tau_mu ~ student_t(3, 0, 1);
    
    // House effect
    lambda ~ normal(0, 1);
    tau_lambda ~ student_t(3, 0, 1);
    
    // Population definition and poll type effects
    beta ~ normal(0, 1);
    nu ~ normal(0, 1);
    
    // EZ and CT adjustment
    gamma_z_tilde ~ normal(0, 1);
    gamma_t_tilde ~ normal(0, 1);
  
  
  // Likelihood
  for (i in 1:N)
    target += beta_proportion_lpdf(vshare[i] | theta[i], tot_eff[i]);
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  real<lower=0,upper=1> prob[D];
  for (d in 1:D)
      prob[d] = inv_logit(alpha0 * d + to_row_vector(alpha) * S[,d]);
      
}
