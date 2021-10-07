data {
  int N;
  int tot_eff[N];
  int vote_eff[N];
  int id_cand[N];
  int C;
  int id_date[N];
  int id_firm[N];
  int F;
  matrix[N,3] X;
  real prior_mu[C];
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int num_basis;
  int D;
  matrix[num_basis,D] B;
}
parameters {
  matrix[num_basis,C] a_raw;
  real a0[C];
  matrix[F,C] lambda;
  matrix[3,C] beta;
  real epsilon[N];
  real<lower=0> sigma_lambda[C];
  real<lower=0> sigma_beta[3];
  real<lower=0> sigma_epsilon;
  real eta[C];
  real<lower=0> tau[C];
}
transformed parameters {
  matrix[num_basis,C] a;
  for (c in 1:C) {
    a[1,c] = prior_mu[c] + a_raw[1,c] * tau[c];
    for (i in 2:num_basis)
      a[i,c] = a[i-1,c] + a_raw[i,c] * tau[c];
  }
}
model {
  // priors
  eta ~ cauchy(0, 2);
  epsilon ~ normal(0, sigma_epsilon);
  sigma_lambda ~ cauchy(0, 2);
  sigma_beta ~ cauchy(0, 2);
  sigma_epsilon ~ cauchy(0, 2);
  tau ~ normal(0, 1);
  a0 ~ normal(0, 1);
  for (c in 1:C) {
    lambda[,c] ~ normal(0, sigma_lambda[c]);
    a_raw[,c] ~ normal(0, 1);
  }
  for (x in 1:3)
    beta[x,] ~ normal(0, sigma_beta[x]);
  
  // likelihood
  for (i in 1:N)
    vote_eff[i] ~ binomial(tot_eff[i],
                           inv_logit(a0[id_cand[i]]*id_date[i] + to_row_vector(a[,id_cand[i]])*B[,id_date[i]] +
                                       lambda[id_firm[i],id_cand[i]] +
                                       X[i,] * beta[,id_cand[i]] +
                                       eta[id_cand[i]] * epsilon[i]));
}
generated quantities {
  matrix[D,C] prob;
  for (d in 1:D)
    for (c in 1:C)
      prob[d,c] = inv_logit(a0[c]*d + to_row_vector(a[,c])*B[,d]);
}