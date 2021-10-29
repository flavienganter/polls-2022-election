# A Poll Aggregator for the 2022 French Presidential Election

(just for fun—work in progress)

## Data

I use data from all voting intention polls fielded since September 1, 2021, based on the survey reports available on the [Commission des sondages website](https://www.commission-des-sondages.fr/notices/). Up to now, I have kept data from all scenarios testing one of the three main Les Républicains candidates—Xavier Bertrand, Valérie Pécresse, and Michel Barnier–and average them; I will adjust in early December, when we know who will eventually be the candidate for Les Républicains. I have also kept scenarios both including and excluding Éric Zemmour (see Model for more details).

## Model

I build on [Heidemanns, Gelman and Morris (2020)](https://hdsr.mitpress.mit.edu/pub/nw1dzd02/release/1) to build a poll aggregator that does just that—aggregating polls—with no prediction intention whatsoever. The model is estimated with Stan.

For each scenario _i_ (part of poll _p_ = _p[i]_) and each candidate _c_, _y<sub>ci</sub>_ is the number of respondents who indicated support for candidate _c_, and _n<sub>i</sub>_ is the total number of respondents supporting any of the candidates tested in scenario _i_. I fit a binomial model:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/binomial.png?raw=true)

where 𝜽<sub>_ci_</sub> is modeled as

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/theta.png?raw=true)

and _date<sub>i</sub>_ is the date, centered so that _date<sub>i</sub>_ = 1 on September 1, 2021.

### Splines

I model the evolutions of voting intentions over time with a spline of degree 3 with _K_ = 4 knots (as of today):

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/spline.png?raw=true)

where (_B<sub>3k</sub>_(.))<sub>_k_</sub> is a sequence of _B_-splines. I define a poll's date as the median day of the fielding period, or as the day immediately following the median when that median is not properly defined. To enforce smoothness and prevent the model from overfitting, I impose a random-walk prior on (𝛼<sub>_ck_</sub>)<sub>_ck_</sub>:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_alpha1.png?raw=true)

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_alpha2.png?raw=true)

### Poll and Organization Effects

In order to partially pool information among the various scenarios of the same poll, I include a candidate-specific poll effect 𝜇<sub>_cp[i]_</sub>, and I also adjust for polling organization effects (𝜆<sub>_co[i]_</sub>). Both are given a hierarchical structure and estimated with weakly informative priors:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_mu.png?raw=true)

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_lambda.png?raw=true)

### Zemmour Adjustment

_z<sub>i</sub>_ is a (standardized) dummy that flags whether Éric Zemmour was among the tested candidates. I include it because not all polls included scenarios that tested Éric Zemmour as candidate in early September. Keeping only scenarios that included Zemmour—and thus completely discarding polls that did not include a Zemmour scenario—would make the September estimates very noisy, and it would be a waste of information. Yet it remains that scenarios without Zemmour are not very informative per se, as he will very likely be candidate, and his inclusion among the candidates significantly changes voting intentions for other candidates (especially Marine Le Pen and the candidate from Les Républicains). The current specification leverages polls that include scenarios both with and without Zemmour to estimate the impact of not including him on the voting intentions for other candidates, and thereby allows me to adjust the results from polls that do not include Zemmour scenarios. I allow the coefficient to vary over time by estimating a distinct coefficient every month, with a random-walk prior:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_gamma1.png?raw=true)

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/prior_gamma2.png?raw=true)

### Other Covariates

The vector _X<sub>i</sub>_ includes log(_n<sub>i</sub>_) to adjust for a potential sample size effect, and two (standardized) dummies that adjust for the subsample of respondents that the polling organization calculated their estimates on (all respondents, only respondents who are absolutely sure that they will vote in April 2022, or an intermediary subsample). To allow the effect of these covariates to vary as the election date gets closer, these coefficients incorporate a time trend:

![](https://github.com/flavienganter/polls-2022-election/blob/main/img/beta.png?raw=true)

## Quantity of Interest

My quantity of interest is the share of the voting intentions that each candidate gets at each date since September 1, 2021, in scenarios where Éric Zemmour is a candidate, and considering not only respondents who are absolutely sure that they will vote on April 10, 2022, but excluding those who are relatively certain that they will _not_ vote. Ideally, I would weight each respondent by how sure they are that they will vote for the first round of the 2022 presidential election, but I do not have the data to implement that strategy.

## Potential Improvements

- Assign each poll with the entire fielding period, instead of just with the median day of the period the poll was fielded
- Propagate the uncertainty related to the truncation of the voting intention share of some "small" candidates (would not make a big difference—mostly aesthetic)
