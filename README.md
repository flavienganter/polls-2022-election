# A Poll Aggregator for the 2022 French Presidential Election

Estimates of the voting intentions for the first round of the 2022 French Presidential Election (medians of the posterior distributions, and 95%, 90%, 80%, and 50% high density credible intervals) on April 8, 2022, two days before the first round:

![](https://github.com/flavienganter/polls-2022-election/blob/main/PollsFrance2022_latest.png?raw=true)

Evolutions of voting intentions since September 2021 (medians of the posterior distributions, and 95% and 50% high density credible intervals) and official results:

![](https://github.com/flavienganter/polls-2022-election/blob/main/PollsFrance2022_evolution_final.png?raw=true)

## Model

I use data from all voting intention polls fielded since September 1, 2021, based on the survey reports available on the [Commission des sondages website](https://www.commission-des-sondages.fr/notices/). I build on [Heidemanns, Gelman and Morris (2020)](https://hdsr.mitpress.mit.edu/pub/nw1dzd02/release/1) to build a poll aggregator that does just that—aggregating polls—with no prediction intention whatsoever. The model is estimated with Stan.

For each scenario $i$ (part of poll $p\ =\ p_{[i]}$) and each candidate $c$, $s_{ci}^{\*}$ is the (adjusted) share of respondents who indicated support for candidate $c$. $s_{ci}^{\*}$ is typically not available, as polling firm round their estimates, so that one can only observe $s_{ci}$. To account for the uncertainty induced by the rounding, I model $s_{ci}^{\*}$ as a latent parameter defined by
$$ s_{ci}^{\*}\ =\ s_{ci}\ +\ \varepsilon_{ci} $$
with $\varepsilon_{ci}\ \sim\ \mathcal{U}\[b_{ci}^l;b_{ci}^u\]$, where $b_{ci}^l$ and $b_{ci}^u$ define the interval around $s_{ci}$ in which $s_{ci}^{\*}$ can be.

Noting $N_i$ the total number of respondents disclosing their voting intentions in the scenario $i$, I model the latent variable $s_{ci}^{\*}$ with a Beta distribution:
$$
s_{ci}^{\*}\ \sim\ \text{Beta}(\theta_{ci}N_i,(1-\theta_{ic})N_i)
$$
where $\theta_{ci}$ is defined as
$$ \theta_{ci}\ \equiv\ \text{logit}^{-1}(\psi_c(date_i)\ +\ \mu_{cp\[i\]}\ +\ \lambda_{ch\[i\]}\ +\ M_i\lambda_{cm\[date_i\]}\ +\ X_i\beta_c(date_i)) $$
and $date_i$ is the date, centered so that $date_i\ =\ 1$ on September 1, 2021.

Two elements motivate the choice of approximating the posterior distribution by a series of Beta regression models, and not by a single multinomial model, which would arguably be the most logical choice:
1. The Beta distribution allows me to model voting intentions directly instead of indirectly, via the number of respondents who indicated support for a given candidate. This is particularly convenient as the only estimates available, $s_{ci}$, is a rounded and adjusted proportion.
2. Having a series of model, rather than one multinomial model, accommodates very easily the fact that not all polls include all current candidates at every point in time.

### Splines

I model the evolutions of voting intentions over time with a spline of degree 3 with $K\ =\ 12$ knots:
$$ \psi_c(date_i)\ =\ \alpha_{c0}\cdot date_i\ +\ \sum_{k=1}^K\alpha_{ck}B_{3k}(date_i) $$
where $(B_{3k}(\cdot))\_k$ is a sequence of $B$-splines. I define a poll's date as the median day of the fielding period, or as the day immediately following the median when that median is not properly defined. To enforce smoothness and prevent the model from overfitting, I impose a random-walk prior on $(\alpha_{ck})\_{ck}$: $\alpha_{c0}\ \sim\ \mathcal{N}(0,\sigma_{\alpha 0})$ and $\alpha_{ck}\ \sim\ \mathcal{N}(\alpha_{ck-1},\tau_{c\alpha})$, $\forall k>0$.

### Poll and House Effects

In order to partially pool information among the various scenarios of the same poll, I include a candidate-specific poll effect $(\mu_{cp\[i\]})$, and I also adjust for house effects $(\lambda_{ch\[i\]})$.

### Zemmour and Taubira Adjustments

$M_i$ is a vector of (standardized) dummies that flag whether Éric Zemmour and Christiane Taubira were among the tested candidates. I include this vector because not all polls included scenarios that tested Éric Zemmour as candidate in September. In the same vein, polls before December 15 did not include scenarios with Christiane Taubira, and not all of them did between December 15 and January 15.

The current specification leverages polls that include scenarios both with and without Zemmour (or Taubira) to estimate the impact of not including him (or her) among the potential candidates, and thereby allows me to adjust the results from polls that do not include Zemmour (or Taubira) scenarios.

### Other Covariates

The vector $X_i$ includes two (standardized) dummies that adjust for the subsample of respondents that the polling firm calculated their estimates on (all respondents, only respondents who are absolutely sure that they will vote in April 2022, or an intermediary subsample), and one additional (standardized) dummy that adjusts for whether the poll is a rolling poll. To allow the effect of these covariates to vary as the election date gets closer, these coefficients (except for the rolling poll dummy, for now) incorporate a time trend:
$$ \beta_c^{(u)}(date_i)\ =\ \beta_{c0}^{(u)}\ +\ \nu_c(date_i\ -\ 1). $$
