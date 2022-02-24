data {
  
  // Poll data
  int N;
  int tot_eff[N];
  real<lower=0,upper=1> vshare_raw[N];
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
  int id_cand[N];
  int C;
  int id_poll[N];
  int P;
  int id_date[N];
  int id_month[N];
  int id_house[N];
  int F;
  matrix[N,3] X;
  
  // EZ and CT adjustement
  int isn_z[N];
  int isn_t[N];
  
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
  real alpha0[C];
  matrix[num_basis,C] alpha_raw;
  real<lower=0> tau_alpha[C];
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_alpha0;
  
  // Covariates
  matrix[P,C] mu;
  real tau_mu_tilde[C];
  matrix[F,C] lambda;
  real tau_lambda_tilde[C];
  matrix[4,C] tau_beta_tilde;
  matrix[3,C] tau_nu_tilde;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_beta[4];
  real<lower=0> sigma_nu[3];
  
  // EZ and CT adjustment
  real gamma_z_tilde[C-2];
  real gamma_t_tilde[C-1];
  real tau_gamma_z;
  real tau_gamma_t;
  
  // Rounding error
  real<lower=-0.005,upper=0.005> epsilon1[N_1];
  real<lower=-0.01,upper=0.01> epsilon2[N_2];
  real<lower=0.001,upper=0.015> epsilon3[N_3];
  real<lower=0.001,upper=0.01> epsilon4[N_4];
  real<lower=0.001,upper=0.005> epsilon5[N_5];
  
}
transformed parameters {
  
  matrix[num_basis,C] alpha;
  real<lower=0> tau_mu[C];
  real<lower=0> tau_lambda[C];
  matrix[3,C] beta;
  matrix[2,C] nu;
  real gamma_z[C];
  real gamma_t[C];
  real<lower=0,upper=1> theta[N];
  real epsilon[N];
  real<lower=0,upper=1> vshare[N];
  
  // Spline coefficients, specified as a random walk
  // to avoid overfit
  for (c in 1:C) {
    alpha[1,c] = alpha_raw[1,c];
    for (i in 2:num_basis)
      alpha[i,c] = alpha[i-1,c] + alpha_raw[i,c] * tau_alpha[c];
  }
  
  // Uncentered parametrization
  for (c in 1:C) {
    tau_mu[c] = exp(sigma_mu * tau_mu_tilde[c]); 
    tau_lambda[c] = exp(sigma_lambda * tau_lambda_tilde[c]);
    for (x in 1:3)
      beta[x,c] = sigma_beta[x] * tau_beta_tilde[x,c];
    for (x in 1:2)
      nu[x,c] = sigma_nu[x] * tau_nu_tilde[x,c];
    if (c < 11) {
      gamma_z[c] = tau_gamma_z * gamma_z_tilde[c];
      gamma_t[c] = tau_gamma_t * gamma_t_tilde[c];
    } else if (c == 11) {
      gamma_z[c] = 0;
      gamma_t[c] = tau_gamma_t * gamma_t_tilde[c];
    } else {
      gamma_z[c] = 0;
      gamma_t[c] = 0;
    }
  }
  
  
  for (i in 1:N) {
    
    // Theta
    if (id_month[i] == 1) {
      theta[i] = inv_logit(alpha0[id_cand[i]] * id_date[i] + to_row_vector(alpha[,id_cand[i]]) * S[,id_date[i]] + // Spline
                                  tau_mu[id_cand[i]] * mu[id_poll[i],id_cand[i]] + // Poll effect
                                  tau_lambda[id_cand[i]] * lambda[id_house[i],id_cand[i]] + // House effect
                                  X[i,1] * (beta[1,id_cand[i]] + nu[1,id_cand[i]] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2,id_cand[i]] + nu[2,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3,id_cand[i]] +
                                  gamma_z[id_cand[i]] * isn_z[i]); // EZ omission adjustment
    } else if (id_month[i] < 4 || id_month[i] > 5) {
      theta[i] = inv_logit(alpha0[id_cand[i]] * id_date[i] + to_row_vector(alpha[,id_cand[i]]) * S[,id_date[i]] + // Spline
                                  tau_lambda[id_cand[i]] * lambda[id_house[i],id_cand[i]] + // House effect
                                  X[i,1] * (beta[1,id_cand[i]] + nu[1,id_cand[i]] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2,id_cand[i]] + nu[2,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3,id_cand[i]]);
    } else {
      theta[i] = inv_logit(alpha0[id_cand[i]] * id_date[i] + to_row_vector(alpha[,id_cand[i]]) * S[,id_date[i]] + // Spline
                                  tau_mu[id_cand[i]] * mu[id_poll[i],id_cand[i]] + // Poll effect
                                  tau_lambda[id_cand[i]] * lambda[id_house[i],id_cand[i]] + // House effect
                                  X[i,1] * (beta[1,id_cand[i]] + nu[1,id_cand[i]] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2,id_cand[i]] + nu[2,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3,id_cand[i]] +
                                  gamma_t[id_cand[i]] * isn_t[i]); // CT omission adjustment
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
    vshare[i] = vshare_raw[i] + epsilon[i];
    
    
  }
  
}
model {
  
  // Priors
  alpha0 ~ normal(0, sigma_alpha0);
  tau_alpha ~ normal(0, sigma_alpha);
  sigma_alpha0 ~ student_t(3, 0, 1);
  sigma_alpha ~ student_t(3, 0, 1);
  for (x in 1:3)
    tau_beta_tilde[x,] ~ normal(0, 1);
  for (x in 1:2)
    tau_nu_tilde[x,] ~ normal(0, 1);
  tau_mu_tilde ~ normal(0, 1);
  tau_lambda_tilde ~ normal(0, 1);
  tau_gamma_z ~ normal(0, 1);
  tau_gamma_t ~ normal(0, 1);
  gamma_z_tilde ~ normal(0, 1);
  gamma_t_tilde ~ normal(0, 1);
  sigma_mu ~ student_t(3, 0, 1);
  sigma_lambda ~ student_t(3, 0, 1);
  sigma_beta ~ student_t(3, 0, 1);
  sigma_nu ~ student_t(3, 0, 1);
  for (c in 1:C) {
    alpha_raw[,c] ~ normal(0, 2);
    mu[,c] ~ normal(0, 1);
    lambda[,c] ~ normal(0, 1);
  }
  
  // Likelihood
  for (i in 1:N)
    target += beta_proportion_lpdf(vshare[i] | theta[i], tot_eff[i]);
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  matrix<lower=0,upper=1>[D,C] prob;
  for (c in 1:C)
    for (d in 1:D)
        prob[d,c] = inv_logit(alpha0[c] * d + to_row_vector(alpha[,c]) * S[,d]);
      
}
