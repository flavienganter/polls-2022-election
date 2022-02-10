# A Poll Aggregator for the 2022 French Presidential Election

Current estimates of the voting intentions for the first round of the 2022 French Presidential Election (medians of the posterior distributions, and 95%, 90%, 80%, and 50% high density credible intervals):

![](https://github.com/flavienganter/polls-2022-election/blob/main/PollsFrance2022_latest.png?raw=true)

Evolutions of voting intentions since September 2021 (medians of the posterior distributions, and 95% and 50% high density credible intervals):

![](https://github.com/flavienganter/polls-2022-election/blob/main/PollsFrance2022_evolution.png?raw=true)

## Model

I use data from all voting intention polls fielded since September 1, 2021, based on the survey reports available on the [Commission des sondages website](https://www.commission-des-sondages.fr/notices/). I build on [Heidemanns, Gelman and Morris (2020)](https://hdsr.mitpress.mit.edu/pub/nw1dzd02/release/1) to build a poll aggregator that does just that—aggregating polls—with no prediction intention whatsoever. The model is estimated with Stan.

For each scenario _i_ (part of poll _p_ = _p[i]_) and each candidate _c_, _y<sub>ci</sub>_ is the number of respondents who indicated support for candidate _c_, and _n<sub>i</sub>_ is the total number of respondents supporting any of the candidates included in scenario _i_. I fit a binomial model:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/binomial.png?raw=true)

where 𝜽<sub>_ci_</sub> is modeled as

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/theta.png?raw=true)

and _date<sub>i</sub>_ is the date, centered so that _date<sub>i</sub>_ = 1 on September 1, 2021.

Note that I do not model voting intentions as a single multinomial model, which would arguably be the most logical choice. Instead, I approximate the posterior distribution by estimating a series of binomial models which accommodates very easily the fact that not all polls include all current candidates at every point in time.

### Splines

I model the evolutions of voting intentions over time with a spline of degree 3 with _K_ = 4 knots (as of today):

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/spline.png?raw=true)

where (_B<sub>3k</sub>_(.))<sub>_k_</sub> is a sequence of _B_-splines. I define a poll's date as the median day of the fielding period, or as the day immediately following the median when that median is not properly defined. To enforce smoothness and prevent the model from overfitting, I impose a random-walk prior on (𝛼<sub>_ck_</sub>)<sub>_ck_</sub>:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_alpha1.png?raw=true)

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_alpha2.png?raw=true)

### Poll and Polling Firm Effects

In order to partially pool information among the various scenarios of the same poll, I include a candidate-specific poll effect 𝜇<sub>_cp[i]_</sub>, and I also adjust for polling firm effects (𝜆<sub>_co[i]_</sub>). Both effects are given a hierarchical structure and estimated with weakly informative priors:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_mu.png?raw=true)

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_lambda.png?raw=true)

### Zemmour and Taubira Adjustments

_z<sub>i</sub>_ is a (standardized) dummy that flags whether Éric Zemmour was among the tested candidates. I include it because not all polls included scenarios that tested Éric Zemmour as candidate in September. Keeping only scenarios that included Zemmour—and thus completely discarding polls that did not include a Zemmour scenario—would make the September estimates very noisy, and it would be a waste of information. Yet it remains that scenarios without Zemmour are not very informative per se, as he will very likely be candidate, and his inclusion among the candidates significantly changes voting intentions for other candidates. The current specification leverages September polls that include scenarios both with and without Zemmour to estimate the impact of not including him on the voting intentions for other candidates, and thereby allows me to adjust the results from polls that do not include Zemmour scenarios.

In the same vein, polls before December 15 did not include scenarios with Christiane Taubira, and not all of them did between December 15 and January 15. I adjust estimates similarly as for Éric Zemmour.

### Other Covariates

The vector _X<sub>i</sub>_ includes log(_n<sub>i</sub>_) to adjust for a potential sample size effect, two (standardized) dummies that adjust for the subsample of respondents that the polling firm calculated their estimates on (all respondents, only respondents who are absolutely sure that they will vote in April 2022, or an intermediary subsample), and one additional (standardized) dummy that adjusts for whether the poll is a rolling poll. To allow the effect of these covariates to vary as the election date gets closer, these coefficients (except for the rolling poll dummy) incorporate a time trend:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/beta.png?raw=true)
