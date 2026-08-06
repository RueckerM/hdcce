# hdcce

Estimation and inference for high-dimensional panel data models with interactive
fixed effects.

The package implements the lasso-type high-dimensional common correlated effects
(HD-CCE) estimator and its desparsified counterpart, which yields asymptotically
valid confidence intervals for individual coefficients, from

> Rücker, M., Vogt, M., Linton, O. and Walsh, C. (2025). Estimation and inference
> in high-dimensional panel data models with interactive fixed effects.
> *Quantitative Economics* **16**(4), 1457–1509.
> [doi:10.3982/QE2308](https://doi.org/10.3982/QE2308)

together with the HD-CCE estimator and significance test for a parametric 
additive specification from

> Rücker, M., Vogt, M. and Linton, O. (2026). High-dimensional panel data models
> with interactive fixed effects: beyond the linear case.
> [arXiv:2608.02055](https://arxiv.org/abs/2608.02055)

## Installation

```r
install.packages("hdcce")

# development version
# remotes::install_github("RueckerM/hdcce")
```

## The model

For units `i = 1, ..., N` and periods `t = 1, ..., T`,

```
Y_it = X_it' beta + gamma_i' F_t + eps_it
```
as in Rücker et al. (2025) with `p` regressors, `K` unobserved factors `F_t`, and unit-specific loadings
`gamma_i`. Both `p` and the number of non-zero coefficients may be large relative
to `NT`. The factors are removed by projecting on the space spanned by the
cross-sectional averages of the regressors, after which the slope vector is
estimated by lasso.

## Data layout

**The panel must be sorted by unit.** Rows `(i-1)*T + 1, ..., i*T` of `x` and `y`
belong to unit `i`:

```
row 1        unit 1, period 1
...
row T        unit 1, period T
row T+1      unit 2, period 1
...
```

This is not checked, and a panel sorted the other way will return numbers rather
than an error.

## Estimation

```r
library(hdcce)
data("data_estimation")

fit <- hdcce_estimator(data_estimation, obs_N = 20, obs_T = 20)

fit$K_hat                    # 3
head(fit$coefs, 6)           # 0.925 0.911 0.877 1.054 0.000 0.000
```

The number of factors is chosen data-driven, as the count of normalised
eigenvalues of the covariance matrix of the averaged regressors exceeding
`TRUNC`; pass `NFACTORS` to fix it instead, and `scree_plot = TRUE` to inspect
the eigenvalue decay.

## Inference

`hdcce_inference` returns a desparsified estimate, its standard error and a
confidence interval for each coefficient in `COEF_INDEX_VEC`, one entry per
index in `$results`:

```r
data("data_inference")

dat <- list(x = data_inference$x, y = data_inference$y[, 1])
inf <- hdcce_inference(dat, obs_N = 20, obs_T = 20, COEF_INDEX_VEC = 1)

inf$results[["1"]]$confidence_band
#>            conf_band_min conf_band_max
#> alpha=0.01       -0.1660        0.0872
#> alpha=0.05       -0.1357        0.0569
#> alpha=0.1        -0.1202        0.0414
```

Rows of `confidence_band` follow the order of `alpha`, which defaults to
`c(0.01, 0.05, 0.10)`.

## Estimation and Inference in an additive model

Supplying a dictionary switches to the parametric additive model from Rücker et 
al. (2026)

```
Y_it = sum_j m_j(X_itj) + gamma_i' F_t + eps_it,     m_j(x) = phi_j(x)' beta_j
```

and provides estimation for the components `m_j` and tests `H0: m_j = 0` with a 
max-type statistic whose critical values come from a Gaussian coupling. 
`Phi` holds the transformed regressors and `group` 
records which covariate each column belongs to.

In this branch `$se` and `$confidence_band` are `NULL`, and `$profile` reports
the statistic at each location `w` together with `n_eff`, the number of nodewise
residuals falling in that bump. If `n_eff` is small the Gaussian approximation is
unreliable — increase the bandwidth `h` or reduce `C`.

## Simulated data

`generate_data` reproduces the design of the simulation study and was used to
build the two data sets shipped with the package:

```r
sim <- generate_data(obs_N = 20, obs_T = 20, p = 61,
                     mu = rep(1, 3 + 61 - 1), RHO = 0.25)

sim$data_estimation    # $x, $y
sim$data_inference     # $x, and $y with one column per signal c** = 0, 0.1, 0.2
```

`mu` has length `K + p - 1` with `K = 3`, since it is the mean of the factor
loadings rather than of the coefficients.

## A note on identification

The factor space is recovered from the cross-sectional averages of the
regressors, which requires the mean loading matrix to have rank `K`. If it does
not, `K_hat` collapses towards one, the estimated error standard deviation
inflates, and the procedures lose power. Comparing `fit$K_hat` against a scree
plot is a quick check.

## Citation

```r
citation("hdcce")
```

## License

GPL (>= 2)
