data {
  
  // Poll data
  int N;
  int tot_eff[N];
  int vote_eff[N];
  int id_cand[N];
  int C;
  int id_poll[N];
  int id_date[N];
  int id_month[N];
  int M;
  int id_firm[N];
  int F;
  matrix[N,4] X;
  
  // EZ adjustement
  int isn_z[N];
  
  // Splines
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int num_basis;
  int D;
  matrix[num_basis,D] B;
  
}
parameters {
  
  // Splines
  real alpha0[C];
  matrix[num_basis,C] alpha_raw;
  real<lower=0> tau_alpha[C];
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_alpha0;
  
  // Covariates
  matrix[12,C] mu;
  real tau_mu_tilde[C];
  matrix[F,C] lambda;
  real tau_lambda_tilde[C];
  matrix[4,C] tau_beta_tilde;
  matrix[3,C] tau_nu_tilde;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_beta[3];
  real<lower=0> sigma_nu[3];
  
  // EZ adjustment
  real gamma_tilde[C-1];
  real tau_gamma;
  
}
transformed parameters {
  
  matrix[num_basis,C] alpha;
  real<lower=0> tau_mu[C];
  real<lower=0> tau_lambda[C];
  matrix[4,C] beta;
  matrix[3,C] nu;
  real gamma[C];
  
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
    for (x in 1:4)
      beta[x,c] = sigma_beta[x] * tau_beta_tilde[x,c];
    for (x in 1:3)
      nu[x,c] = sigma_nu[x] * tau_nu_tilde[x,c];
    if (c < 12) {
      gamma[c] = tau_gamma * gamma_tilde[c];
    } else {
      gamma[c] = 0;
    }
  }
  
}
model {
  
  // Priors
  alpha0 ~ normal(0, sigma_alpha0);
  tau_alpha ~ normal(0, sigma_alpha);
  sigma_alpha0 ~ student_t(3, 0, 1);
  sigma_alpha ~ student_t(3, 0, 1);
  for (x in 1:4)
    tau_beta_tilde[x,] ~ normal(0, 1);
  for (x in 1:3)
    tau_nu_tilde[x,] ~ normal(0, 1);
  tau_mu_tilde ~ normal(0, 1);
  tau_lambda_tilde ~ normal(0, 1);
  tau_gamma ~ normal(0, 1);
  gamma_tilde ~ normal(0, 1);
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
  for (i in 1:N) {
    if (id_month[i] == 2) {
      target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  alpha0[id_cand[i]] * id_date[i] + to_row_vector(alpha[,id_cand[i]]) * B[,id_date[i]] + // Spline
                                  tau_mu[id_cand[i]] * mu[id_poll[i],id_cand[i]] + // Poll effect
                                  tau_lambda[id_cand[i]] * lambda[id_firm[i],id_cand[i]] + // Firm effect
                                  X[i,1] * (beta[1,id_cand[i]] + nu[1,id_cand[i]] * (id_date[i] - 1)) + // Sample size and population definition effects
                                  X[i,2] * (beta[2,id_cand[i]] + nu[2,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,3] * (beta[3,id_cand[i]] + nu[3,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,4] * beta[4,id_cand[i]] +
                                  gamma[id_cand[i]] * isn_z[i]); // EZ omission adjustment
    } else {
      target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  alpha0[id_cand[i]] * id_date[i] + to_row_vector(alpha[,id_cand[i]]) * B[,id_date[i]] + // Spline
                                  tau_lambda[id_cand[i]] * lambda[id_firm[i],id_cand[i]] + // Firm effect
                                  X[i,1] * (beta[1,id_cand[i]] + nu[1,id_cand[i]] * (id_date[i] - 1)) + // Sample size and population definition effects
                                  X[i,2] * (beta[2,id_cand[i]] + nu[2,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,3] * (beta[3,id_cand[i]] + nu[3,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,4] * beta[4,id_cand[i]]);
    }
  }
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  matrix<lower=0,upper=1>[D,C] prob;
  for (d in 1:D)
    for (c in 1:C)
      prob[d,c] = inv_logit(alpha0[c] * d + to_row_vector(alpha[,c]) * B[,d] + beta[3,c] + nu[3,c] * (d - 1));
      
}
