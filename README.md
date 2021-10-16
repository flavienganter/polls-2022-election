# A Poll Aggregator for the 2022 French Presidential Election

(just for fun—work in progress)

## Data

I use data from all voting intention polls fielded since September 1, 2021, based on the survey reports available on the [Commission des sondages website](https://www.commission-des-sondages.fr/notices/). Up to now, I have kept data from all scenarios testing one of the three main Les Républicains candidates—Xavier Bertrand, Valérie Pécresse, and Michel Barnier–and I will adjust in early December, when we know who will eventually be the candidate for Les Républicains. I have also kept scenarios both including and excluding Éric Zemmour (see Model for more details).

## Model

I build on [Heidemanns, Gelman and Morris (2020)](https://hdsr.mitpress.mit.edu/pub/nw1dzd02/release/1) to build a poll aggregator that does just that—aggregating polls—with no prediction intention whatsoever. The model is estimated with Stan.

For each scenario _i_ (part of poll _p_ = _p[i]_) and each candidate _c_, _y<sub>ci</sub>_ is the number of respondents who indicated support for candidate _c_, and _n<sub>i</sub>_ is the total number of respondents supporting any of the candidates tested in scenario _i_. I fit a binomial model:
![](https://github.com/flavienganter/polls-2022-election/blob/main/img/binomial.png?raw=true)
where 𝜽<sub>_ci_</sub> is modeled as
