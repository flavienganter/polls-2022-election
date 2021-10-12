data {
  int N;
  int tot_eff[N];
  int vote_eff[N];
  int id_cand[N];
  int C;
  int id_poll[N];
  int P;
  int isn_z[N];
  real prior_mu[C];
}
parameters {
  real mu[C];
  real lambda[P];
  real kappa[C];
}
model {
  // priors
  for (c in 1:C)
    mu[c] ~ normal(prior_mu[c], 2);
  lambda ~ normal(0, 3);
  kappa ~ normal(0, 2);
  
  // likelihood
  for (i in 1:N)
    vote_eff[i] ~ binomial(tot_eff[i],
                           inv_logit(mu[id_cand[i]] +
                                       kappa[id_cand[i]] * isn_z[i] +
                                       lambda[id_poll[i]]));
}