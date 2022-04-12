# A Poll Aggregator for the 2022 French Presidential Election

Estimates of the voting intentions for the first round of the 2022 French Presidential Election (medians of the posterior distributions, and 95%, 90%, 80%, and 50% high density credible intervals) on April 8, 2022, two days before the first round:

![](https://github.com/flavienganter/polls-2022-election/blob/main/PollsFrance2022_latest.png?raw=true)

Evolutions of voting intentions since September 2021 (medians of the posterior distributions, and 95% and 50% high density credible intervals) and official results:

![](https://github.com/flavienganter/polls-2022-election/blob/main/PollsFrance2022_evolution_final.png?raw=true)

## Model

I use data from all voting intention polls fielded since September 1, 2021, based on the survey reports available on the [Commission des sondages website](https://www.commission-des-sondages.fr/notices/). I build on [Heidemanns, Gelman and Morris (2020)](https://hdsr.mitpress.mit.edu/pub/nw1dzd02/release/1) to build a poll aggregator that does just that—aggregating polls—with no prediction intention whatsoever. The model is estimated with Stan.

For each scenario _i_ (part of poll _p_ = _p[i]_) and each candidate _c_, _s<sup>*</sup><sub>ci</sub>_ is the (adjusted) share of respondents who indicated support for candidate _c_. _s<sup>*</sup><sub>ci</sub>_ is typically not available, as polling firm round their estimates, so that one can only observe _s<sub>ci</sub>_. To account for the uncertainty induced by the rounding, I model _s<sup>*</sup><sub>ci</sub>_ as a latent parameter defined by

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/latent_s.png?raw=true)

with

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/epsilon.png?raw=true)

where _b<sub>ci</sub><sup>l</sup>_ and _b<sub>ci</sub><sup>u</sup>_ define the interval around _s<sub>ci</sub>_ in which _s<sup>*</sup><sub>ci</sub>_ can be.

Noting _N<sub>i</sub>_ the total number of respondents disclosing their voting intentions in the scenario _i_, I model the latent variable _s<sup>*</sup><sub>ci</sub>_ with a Beta distribution:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/beta_m.png?raw=true)

where 𝜽<sub>_ci_</sub> is defined as:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/theta.png?raw=true)

and _date<sub>i</sub>_ is the date, centered so that _date<sub>i</sub>_ = 1 on September 1, 2021.

Two elements motivate the choice of approximating the posterior distribution by a series of Beta regression models, and not by a single multinomial model, which would arguably be the most logical choice:
1. The Beta distribution allows me to model voting intentions directly instead of indirectly, via the number of respondents who indicated support for a given candidate. This is particularly convenient as the only estimates available, _s<sub>ci</sub>_, is a rounded and adjusted proportion.
2. Having a series of model, rather than one multinomial model, accommodates very easily the fact that not all polls include all current candidates at every point in time.

### Splines

I model the evolutions of voting intentions over time with a spline of degree 3 with _K_ = 4 knots (as of today):

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/spline.png?raw=true)

where (_B<sub>3k</sub>_(.))<sub>_k_</sub> is a sequence of _B_-splines. I define a poll's date as the median day of the fielding period, or as the day immediately following the median when that median is not properly defined. To enforce smoothness and prevent the model from overfitting, I impose a random-walk prior on (𝛼<sub>_ck_</sub>)<sub>_ck_</sub>:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_alpha1.png?raw=true)

### Poll and House Effects

In order to partially pool information among the various scenarios of the same poll, I include a candidate-specific poll effect 𝜇<sub>_cp[i]_</sub>, and I also adjust for house effects (𝜆<sub>_ch[i]_</sub>).

### Zemmour and Taubira Adjustments

_M<sub>i</sub>_ is a vector of (standardized) dummies that flag whether Éric Zemmour and Christiane Taubira were among the tested candidates. I include this vector because not all polls included scenarios that tested Éric Zemmour as candidate in September. In the same vein, polls before December 15 did not include scenarios with Christiane Taubira, and not all of them did between December 15 and January 15.

The current specification leverages polls that include scenarios both with and without Zemmour (or Taubira) to estimate the impact of not including him (or her) among the potential candidates, and thereby allows me to adjust the results from polls that do not include Zemmour (or Taubira) scenarios.

### Other Covariates

The vector _X<sub>i</sub>_ includes two (standardized) dummies that adjust for the subsample of respondents that the polling firm calculated their estimates on (all respondents, only respondents who are absolutely sure that they will vote in April 2022, or an intermediary subsample), and one additional (standardized) dummy that adjusts for whether the poll is a rolling poll. To allow the effect of these covariates to vary as the election date gets closer, these coefficients (except for the rolling poll dummy, for now) incorporate a time trend:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/beta.png?raw=true)
