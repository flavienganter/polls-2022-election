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
  int id_month[N];
  int M;
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
  real<lower=0> sigma_a0;
  
  // Covariates
  matrix[P,C] mu;
  real<lower=0> tau_mu[C];
  matrix[F,C] lambda;
  real<lower=0> tau_lambda[C];
  matrix[3,C] beta;
  matrix[3,C] nu;
  real<lower=0> sigma_mu[C];
  real<lower=0> sigma_lambda[C];
  real<lower=0> sigma_beta[3];
  real<lower=0> sigma_nu[3];
  
  // EZ adjustment
  matrix[M,C] gamma;
  real<lower=0> tau_gamma[C];
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_gamma0;
  
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
  a0 ~ normal(0, sigma_a0);
  tau_a ~ normal(0, sigma_a);
  sigma_a0 ~ student_t(3, 0, 1);
  sigma_a ~ student_t(3, 0, 1);
  for (x in 1:3) {
    beta[x,] ~ normal(0, sigma_beta[x]);
    nu[x,] ~ normal(0, sigma_nu[x]);
  }
  tau_mu ~ normal(0, sigma_mu);
  tau_lambda ~ normal(0, sigma_lambda);
  tau_gamma ~ normal(0, sigma_gamma);
  sigma_mu ~ student_t(3, 0, 1);
  sigma_lambda ~ student_t(3, 0, 1);
  sigma_beta ~ student_t(3, 0, 1);
  sigma_nu ~ student_t(3, 0, 1);
  sigma_gamma0 ~ student_t(3, 0, 1);
  sigma_gamma ~ student_t(3, 0, 1);
  for (c in 1:C) {
    a_raw[,c] ~ normal(0, 2);
    mu[,c] ~ normal(0, 1);
    lambda[,c] ~ normal(0, 1);
    if (c < 12) {
      gamma[1,c] ~ normal(0, sigma_gamma0);
      for (m in 2:M)
        gamma[m,c] ~ normal(gamma[1,c], tau_gamma);
    } else {
      gamma[,12] ~ normal(0, .01);
    }
  }
  
  // Likelihood
  for (i in 1:N)
    target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  a0[id_cand[i]] * id_date[i] + to_row_vector(a[,id_cand[i]]) * B[,id_date[i]] + // Spline
                                  tau_mu[id_cand[i]] * mu[id_poll[i],id_cand[i]] + // Poll effect
                                  tau_lambda[id_cand[i]] * lambda[id_firm[i],id_cand[i]] + // Firm effect
                                  X[i,1] * (beta[1,id_cand[i]] + nu[1,id_cand[i]] * (id_date[i] - 1)) + // Sample size and population definition effects
                                  X[i,2] * (beta[2,id_cand[i]] + nu[2,id_cand[i]] * (id_date[i] - 1)) +
                                  X[i,3] * (beta[3,id_cand[i]] + nu[3,id_cand[i]] * (id_date[i] - 1)) +
                                  gamma[id_month[i],id_cand[i]] * isn_z[i]); // EZ omission adjustment
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  matrix<lower=0,upper=1>[D,C] prob;
  for (d in 1:D)
    for (c in 1:C)
      prob[d,c] = inv_logit(a0[c] * d + to_row_vector(a[,c]) * B[,d] + beta[3,c] + nu[3,c] * (d - 1));
      
}