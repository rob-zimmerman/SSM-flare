# SSM-flare

R code to accompany:

> Zimmerman, R., van Dyk, D. A., Kashyap, V. L., & Siemiginowska, A. (2024). *Separating states in astronomical sources using hidden Markov models: with a case study of flaring and quiescence on EV Lac*. Monthly Notices of the Royal Astronomical Society. 534(3):2142–2167. [doi:10.1093/mnras/stae2082](https://doi.org/10.1093/mnras/stae2082)

This repository contains analysis code for separating quiescent versus flaring intervals in multiband (soft/hard) Poisson light curves using a two-stage estimation approach:
1. **Stage 1:** fit a continuous-state hidden Markov model (HMM) for a latent continuous state process $X_t$ that drives expected Poisson count rates in multiple passbands.
2. **Stage 2:** dichotomize the fitted continuous states with a finite mixture model to classify time intervals into quiescence versus activity (and estimate the long-run flaring fraction, flare intervals, and enable per-state spectral comparisons).

The paper provides full modelling details, inference strategies, and astrophysical interpretations.

---

## Method summary

- Observations are binned into discrete time intervals; the observed data are count light curves in multiple passbands (soft/hard), modelled in the Poisson regime.
- The latent process is a continuous-space Markov chain (several autoregressive models are considered in the paper), yielding a continuous-state HMM.
- Because continuous-state HMM likelihoods can be computationally challenging, the paper develops a discretized approximation that can be made arbitrarily accurate.
- The fitted/decoded latent states are classified into quiescent versus flaring intervals using a finite mixture model on the marginal fitted states.

---

## In this repo

At the top level you’ll find:

- `01885/` — analysis for EV Lac ObsID **01885** (one epoch).
- `10679/` — analysis for EV Lac ObsID **10679** (another epoch).

These ObsIDs correspond to the two EV Lac datasets analyzed in the paper. Each ObsID folder contains the scripts/data needed to reproduce that epoch’s results. A typical workflow is:

1. Load and prepare light curves $\boldsymbol{Y}_t$ (soft/hard bands, time binning, etc.)
2. Stage 1 fit: fit the continuous-space HMM and obtain filtered/smoothed state estimates $\hat{X}_t$
3. Stage 2 fit: fit a finite mixture model to the $\hat{X}_t$ and classify time bins into quiescent versus flaring
4. Summaries/plots: flare fraction, inferred flare intervals and diagnostic plots 
