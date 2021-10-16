data {
  
  // Poll data
  int N;
  int tot_eff[N];
  int vote_eff[N];
  int id_cand[N];
  int C;
  int id_poll[N];
  int P;
  int id_date[N];
  int id_firm[N];
  int F;
  matrix[N,3] X;
  
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
  real a0[C];
  matrix[num_basis,C] a_raw;
  real<lower=0> tau_a[C];
  real<lower=0> sigma_a;
  
  // Covariates
  matrix[P,C] mu;
  real<lower=0> tau_mu[C];
  matrix[F,C] lambda;
  real<lower=0> tau_lambda[C];
  matrix[3,C] beta;
  real<lower=0> sigma_mu[C];
  real<lower=0> sigma_lambda[C];
  real<lower=0> sigma_beta[3];
  
  // EZ adjustment
  real gamma[C];
  real<lower=0> sigma_gamma;
  
}
transformed parameters {
  
  // Spline coefficients, specified as a random walk
  // to avoid overfit
  matrix[num_basis,C] a;
  for (c in 1:C) {
    a[1,c] = a_raw[1,c];
    for (i in 2:num_basis)
      a[i,c] = a[i-1,c] + a_raw[i,c] * tau_a[c];
  }
  
}
model {
  
  // Priors
  a0 ~ normal(0, 2);
  tau_a ~ normal(0, sigma_a);
  sigma_a ~ student_t(3, 0, 1);
  for (x in 1:3)
    beta[x,] ~ normal(0, sigma_beta[x]);
  tau_mu ~ normal(0, sigma_mu);
  tau_lambda ~ normal(0, sigma_lambda);
  sigma_mu ~ student_t(3, 0, 1);
  sigma_lambda ~ student_t(3, 0, 1);
  sigma_beta ~ student_t(3, 0, 1);
  sigma_gamma ~ student_t(3, 0, 1);
  for (c in 1:C) {
    a_raw[,c] ~ normal(0, 2);
    mu[,c] ~ normal(0, 1);
    lambda[,c] ~ normal(0, 1);
    if (c < 12) {
      gamma[c] ~ normal(0, sigma_gamma);
    } else {
      gamma[12] ~ normal(0, .01);
    }
  }
  
  // Likelihood
  for (i in 1:N)
    target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  a0[id_cand[i]] * id_date[i] + to_row_vector(a[,id_cand[i]]) * B[,id_date[i]] + // Spline
                                  tau_mu[id_cand[i]] * mu[id_poll[i],id_cand[i]] + // Poll effect
                                  tau_lambda[id_cand[i]] * lambda[id_firm[i],id_cand[i]] + // Firm effect
                                  X[i,] * beta[,id_cand[i]] + // Sample size and population definition effects
                                  gamma[id_cand[i]] * isn_z[i]); // EZ omission adjustment
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  matrix<lower=0,upper=1>[D,C] prob;
  for (d in 1:D)
    for (c in 1:C)
      prob[d,c] = inv_logit(a0[c] * d + to_row_vector(a[,c]) * B[,d] + beta[3,c]);
      
}