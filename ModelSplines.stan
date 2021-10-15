data {
  
  // Poll data
  int N;
  int tot_eff[N];
  int vote_eff[N];
  int id_cand[N];
  int C;
  int id_date[N];
  int id_firm[N];
  int F;
  matrix[N,3] X;
  
  // EZ adjustement
  int isn_z[N];
  real zemcorr_mean[C-1];
  real zemcorr_sd[C-1];
  
  // Prior from previous polls
  real prior_mu[C];
  
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
  matrix[num_basis,C] a_raw;
  real a0[C];
  real<lower=0> tau[C];
  
  // Covariates
  matrix[F,C] lambda;
  matrix[3,C] beta;
  real<lower=0> sigma_lambda[C];
  real<lower=0> sigma_beta[3];
  real zemm_corr[C];
  real eta[C];
  
  // Uncertainty
  real epsilon[N];
  real<lower=0> sigma_epsilon;
  
}
transformed parameters {
  
  // Spline coefficients, specified as a random walk
  // to avoid overfit
  matrix[num_basis,C] a;
  for (c in 1:C) {
    a[1,c] = prior_mu[c] + a_raw[1,c] * tau[c];
    for (i in 2:num_basis)
      a[i,c] = a[i-1,c] + a_raw[i,c] * tau[c];
  }
  
}
model {
  
  // Priors
  a0 ~ normal(0, 1);
  tau ~ normal(0, 1);
  for (x in 1:3)
    beta[x,] ~ normal(0, sigma_beta[x]);
  sigma_lambda ~ cauchy(0, 2);
  sigma_beta ~ cauchy(0, 2);
  for (c in 1:C) {
    a_raw[,c] ~ normal(0, 1);
    lambda[,c] ~ normal(0, sigma_lambda[c]);
      if (c < 12) {
        zemm_corr[c] ~ normal(zemcorr_mean, zemcorr_sd);
      } else {
        zemm_corr[12] ~ normal(0, .01);
      }
  }
  eta ~ cauchy(0, 2);
  epsilon ~ normal(0, sigma_epsilon);
  sigma_epsilon ~ cauchy(0, 2);
  
  // Likelihood
  for (i in 1:N)
    vote_eff[i] ~ binomial(tot_eff[i],
                           inv_logit(a0[id_cand[i]]*id_date[i] + to_row_vector(a[,id_cand[i]])*B[,id_date[i]] + // Spline
                                       lambda[id_firm[i],id_cand[i]] + // Firm effect
                                       X[i,] * beta[,id_cand[i]] + // Sample size and population definition effects
                                       isn_z[i] * zemm_corr[id_cand[i]] + // EZ omission adjustment
                                       eta[id_cand[i]] * epsilon[i])); // Additional uncertainty
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  matrix[D,C] prob;
  for (d in 1:D)
    for (c in 1:C)
      prob[d,c] = inv_logit(a0[c]*d + to_row_vector(a[,c])*B[,d] + beta[3,c]);
      
}