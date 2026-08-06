# hdcce

Estimation and inference for high-dimensional panel data models with interactive
fixed effects.

> Rücker, M., Vogt, M., Linton, O. and Walsh, C. (2025). Estimation and inference
> in high-dimensional panel data models with interactive fixed effects.
> *Quantitative Economics* **16**(4), 1457–1509.
> [doi:10.3982/QE2308](https://doi.org/10.3982/QE2308)
>
> Rücker, M., Vogt, M. and Linton, O. (2026). High-dimensional panel data models
> with interactive fixed effects: beyond the linear case.
> [arXiv:2608.02055](https://arxiv.org/abs/2608.02055)

## Installation

```r
install.packages("hdcce")
# remotes::install_github("RueckerM/hdcce")   # development version
```

---

## 1. The linear model

For units $i = 1,\dots,n$ and periods $t = 1,\dots,T$,

$$Y_{it} = X_{it}^\top \beta + \gamma_i^\top F_t + \varepsilon_{it},$$

with $X_{it} \in \mathbb{R}^p$, unobserved factors $F_t \in \mathbb{R}^K$ and
unit-specific loadings $\gamma_i \in \mathbb{R}^K$. Stacking over $t$, with
$Y_i \in \mathbb{R}^T$, $\boldsymbol X_i \in \mathbb{R}^{T\times p}$ and
$\boldsymbol F \in \mathbb{R}^{T \times K}$,

$$Y_i = \boldsymbol X_i \beta + \boldsymbol F \gamma_i + \varepsilon_i .$$

The dimension $p$ may exceed $nT$; $\beta$ is assumed sparse. Neither
$\boldsymbol F$ nor $\gamma_i$ is estimated.

## 2. Estimation — `hdcce_estimator`

**Step 1: the factor space.** Let
$\bar{\boldsymbol X} = n^{-1}\sum_{i=1}^n \boldsymbol X_i \in \mathbb{R}^{T\times p}$
and $\widehat{\boldsymbol\Sigma} = \bar{\boldsymbol X}^\top\bar{\boldsymbol X}/T$,
with eigenvalues $\hat\lambda_1 \ge \dots \ge \hat\lambda_p \ge 0$ and
eigenvectors $\widehat U_{(1)},\dots,\widehat U_{(p)}$. The number of factors is
estimated by

$$\widehat K = \\#\\{\, j : \hat\lambda_j / \hat\lambda_1 > \tau \,\\},$$

with the truncation $\tau$ set by `TRUNC` (default `0.01`), or fixed directly
through `NFACTORS`. Writing
$\widehat{\boldsymbol U} = (\widehat U_{(1)} \cdots \widehat U_{(\widehat K)})$
and $\widehat{\boldsymbol W} = \bar{\boldsymbol X}\widehat{\boldsymbol U} \in
\mathbb{R}^{T\times\widehat K}$, the projection is

$$\widehat{\boldsymbol\Pi} = \boldsymbol I_T -
  \widehat{\boldsymbol W}\big(\widehat{\boldsymbol W}^\top\widehat{\boldsymbol W}\big)^{-}
  \widehat{\boldsymbol W}^\top ,$$

computed from an orthonormal basis of $\widehat{\boldsymbol W}$, so that
$\widehat{\boldsymbol\Pi}$ is exactly idempotent and
$\widehat{\boldsymbol\Pi}\boldsymbol F \approx 0$.

**Step 2: penalised least squares on the projected data.** With
$\widetilde Y_i = \widehat{\boldsymbol\Pi}Y_i$ and
$\widetilde{\boldsymbol X}_i = \widehat{\boldsymbol\Pi}\boldsymbol X_i$,

$$\widehat\beta_\lambda \in \arg\min_{b\in\mathbb{R}^p}
  \left\\{ \frac{1}{nT}\sum_{i=1}^n \big\\|\widetilde Y_i -
  \widetilde{\boldsymbol X}_i b\big\\|^2 + \lambda\\|b\\|_1 \right\\},$$

with $\lambda$ chosen by cross-validation over folds that partition the
**cross-section** (`NFOLDS`, or `foldid`), or supplied through `lambda`. Setting
`variant = "LS"` replaces the penalty by ordinary least squares, which requires
$p < nT$.

> `glmnet` minimises $(2nT)^{-1}\sum_i\|\cdot\|^2 + \lambda_g\|b\|_1$, so
> $\lambda = 2\lambda_g$. The `lambda` argument and the returned `Lambda` are
> both on the scale of the display above.

```r
fit <- hdcce_estimator(data, obs_N = n, obs_T = T)
fit$coefs        # beta_hat
fit$K_hat        # K_hat
fit$Lambda       # selected penalty
fit$eigenvalues  # normalised eigenvalues used for K_hat
```

## 3. Inference — `hdcce_inference`

For a target coordinate $j$ the lasso estimate is debiased along the direction of
the residual of a **nodewise regression** of the $j$-th projected regressor on
the remaining ones:

$$\widehat\theta_\kappa \in \arg\min_{\vartheta}
  \left\\{\frac{1}{nT}\sum_{i=1}^n\big\\|\widetilde X_{i(j)} -
  \widetilde{\boldsymbol X}_{i(-j)}\vartheta\big\\|^2 +
  \kappa\\|\vartheta\\|_1\right\\},
  \qquad
  \hat u_i = \widetilde X_{i(j)} - \widetilde{\boldsymbol X}_{i(-j)}\widehat\theta_\kappa .$$

With lasso residuals
$\hat r_i = \widetilde Y_i - \widetilde{\boldsymbol X}_i\widehat\beta_\lambda$,
the **desparsified estimator** is

$$\widehat\beta^{\text{de}}_j = \widehat\beta_{\lambda,j} +
  \frac{\sum_{i=1}^n \hat u_i^\top \hat r_i}
       {\sum_{i=1}^n \hat u_i^\top \widetilde X_{i(j)}}
  \qquad\text{and}\qquad
  \frac{\widehat\beta^{\text{de}}_j - \beta_j}{\widehat{se}_j}
  \xrightarrow{d} N(0,1).$$

Writing $D = \sum_i \hat u_i^\top \widetilde X_{i(j)}$ and
$\hat\sigma_i^2 = (T-\widehat K)^{-1}\\|\hat r_i\\|^2$, the `HAC` argument selects

| `HAC` | $\widehat{se}_j^{\,2}$ | assumption on $\varepsilon_{it}$ |
|---|---|---|
| `1` | $D^{-2}\,\hat\sigma^2\sum_i\\|\hat u_i\\|^2$, $\hat\sigma^2 = (n T - n\widehat K)^{-1}\sum_i\\|\hat r_i\\|^2$ | homoscedastic, no serial correlation |
| `2` (default) | $D^{-2}\sum_i \hat\sigma_i^2\\|\hat u_i\\|^2$ | heteroscedastic, no serial correlation |
| `3` | $D^{-2}\sum_i \big(\sum_t \hat u_{it}\hat r_{it}\big)^2$ | heteroscedastic and serially correlated |

Confidence intervals are
$\widehat\beta^{\text{de}}_j \pm \widehat{se}_j\,\Phi^{-1}(1-\alpha/2)$, one row
per entry of `alpha`, and the reported $p$-value is
$2\{1-\Phi(|\widehat\beta^{\text{de}}_j|/\widehat{se}_j)\}$ for
$H_0\colon \beta_j = 0$ — so $p \le \alpha$ exactly when $0$ falls outside the
$(1-\alpha)$ interval.

```r
inf <- hdcce_inference(data, obs_N = n, obs_T = T, COEF_INDEX_VEC = c(1, 5))
inf$results[["1"]]$coef_despar
inf$results[["1"]]$se
inf$results[["1"]]$confidence_band
```

## 4. The additive model and the significance test

Supplying `dictionaries` replaces the linear index by

$$Y_{it} = \sum_{j=1}^p m_j(X_{it,j}) + \gamma_i^\top F_t + \varepsilon_{it},
  \qquad m_j(x) = \phi_j(x)^\top\beta_j,$$

where $\phi_j = (\phi_{j1},\dots,\phi_{jL_j})^\top$ is a known vector of
transformations. Collecting them gives $\boldsymbol\Phi_i \in
\mathbb{R}^{T\times d}$ with $d = \sum_j L_j$: `Phi` is the stacked
$nT \times d$ matrix and `group` records, for each of its $d$ columns, which
covariate it transforms.

The test of $H_0\colon m_j = 0$ against $H_1\colon m_j \neq 0$ proceeds as
follows.

**(i)** The projection $\check{\boldsymbol\Pi}$ is built as in §2 but from
$\bar{\boldsymbol X}_{(-j)}$, the cross-sectional average with column $j$
removed. It is therefore recomputed for each tested index.

**(ii)** All quantities are empirically centred across units and projected:
$\check Y_i = \check{\boldsymbol\Pi}(Y_i - \bar Y)$,
$\check{\boldsymbol\Phi}_i = \check{\boldsymbol\Pi}(\boldsymbol\Phi_i - \bar{\boldsymbol\Phi})$,
$\check X_{i(j)} = \check{\boldsymbol\Pi}(X_{i(j)} - \bar X_{(j)})$.

**(iii)** The tested covariate is assumed to obey a nodewise equation
$X_{it,j} = \sum_{j'\neq j}\phi^{\text{nw}}_{j'}(X_{it,j'})^\top\theta_{j'} +
F_t^\top\nu_i + u_{it}$ with sparse $\theta$; a lasso on the projected data
yields residuals $\check u_i$. Supply `Phi_nw` and `group_nw` for a nodewise
dictionary different from `Phi`.

**(iv)** With $\check\beta_\lambda$ the lasso on
$(\check Y, \check{\boldsymbol\Phi})$ and
$R_i = \check Y_i - \check{\boldsymbol\Phi}_{i(-j)}\check\beta_{\lambda,-j}$ the
residual under the null, let $\tau$ be a kernel supported on $[-1,1]$ and
$\tau_{w,h}(u) = \tau\\{(u-w)/h\\}/\sqrt h$. For
$\mathcal W = \\{w \in [-C,C] : w = -C + (2\ell-1)h,\ \ell\in\mathbb N\\}$,

$$\Psi_{w,h} = \frac{\sum_{i=1}^n R_i^\top \tau_{w,h}(\check u_i)}
   {\big\\{\sum_{i=1}^n\\|\check{\boldsymbol\Pi}\tau_{w,h}(\check u_i)\\|^2\big\\}^{1/2}},
  \qquad \Psi = \max_{w\in\mathcal W}\big|\Psi_{w,h}\big| .$$

**(v)** Critical values come from a Gaussian coupling. With

$$d_{i,w} = \frac{\check{\boldsymbol\Pi}\tau_{w,h}(\check u_i)}
  {\big\\{(nT)^{-1}\sum_i\\|\check{\boldsymbol\Pi}\tau_{w,h}(\check u_i)\\|^2\big\\}^{1/2}},
  \qquad
  \hat\sigma^2 = \frac{\sum_i\\|\check Y_i - \check{\boldsymbol\Phi}_i\check\beta_\lambda\\|^2}
                      {n(T-\widehat K)},$$

and $G_{it}$ i.i.d. $N(0,\hat\sigma^2)$ independent of the data,

$$\Psi^{\text{Gauss}} = \max_{w\in\mathcal W}
  \Big|(nT)^{-1/2}\sum_{i,t} d_{it,w}\,G_{it}\Big| ,$$

drawn `B` times with $d_{i,w}$, $\check u_i$ and $\hat\sigma$ held fixed. The
test rejects when $\Psi \ge \hat c_{1-\alpha}$, the empirical $(1-\alpha)$
quantile of $\Psi^{\text{Gauss}}$, and the reported $p$-value is
$(1 + \\#\\{b : \Psi^{\text{Gauss}}_b \ge \Psi\\})/(B+1)$.

```r
tst <- hdcce_inference(data, obs_N = n, obs_T = T, COEF_INDEX_VEC = 1,
                       dictionaries = list(Phi = Phi, group = group))
tst$results[["1"]]$statistic        # Psi
tst$results[["1"]]$critical_values  # c_hat, one per alpha
tst$results[["1"]]$p_value
tst$results[["1"]]$profile          # Psi_{w,h} over the grid, with bump counts
```

In this branch `$se` and `$confidence_band` are `NULL`; in the linear branch
`$profile` is `NULL`. `profile$n_eff` gives the number of nodewise residuals
falling in each bump $[w-h, w+h]$; when the smallest is below roughly 30 the
Gaussian approximation is unreliable and the bandwidth should be increased.

---

## Data layout

**The panel must be sorted by unit**: rows $(i-1)T+1,\dots,iT$ of `x`, `y` and
`Phi` belong to unit $i$. This is not checked, and a panel sorted by period
returns numbers rather than an error.

## Quick start

```r
library(hdcce)
data("data_estimation"); data("data_inference")

## estimation
fit <- hdcce_estimator(data_estimation, obs_N = 20, obs_T = 20)
fit$K_hat                                    # 3
head(fit$coefs, 4)                           # 0.925 0.911 0.877 1.054

## inference on beta_1 (true value 0 in column 1)
dat <- list(x = data_inference$x, y = data_inference$y[, 1])
inf <- hdcce_inference(dat, obs_N = 20, obs_T = 20, COEF_INDEX_VEC = 1)
inf$results[["1"]]$confidence_band
#>            conf_band_min conf_band_max
#> alpha=0.01       -0.1660        0.0872
#> alpha=0.05       -0.1357        0.0569
#> alpha=0.1        -0.1202        0.0414

## test H0: m_1 = 0 with phi_j(x) = (x, x^2)
X   <- data_inference$x; p <- ncol(X)
Phi <- do.call(cbind, lapply(1:p, function(j) cbind(X[, j], X[, j]^2)))
tst <- hdcce_inference(list(x = X, y = data_inference$y[, 1]),
                       obs_N = 20, obs_T = 20, COEF_INDEX_VEC = 1,
                       dictionaries = list(Phi = Phi, group = rep(1:p, each = 2)))
tst$results[["1"]]$p_value
```

## Identification

The factor space is recovered from $\bar{\boldsymbol X}$, which requires the mean
loading matrix to have rank $K$ — the usual CCE rank condition. When it fails,
$\bar{\boldsymbol X}$ has a systematic part of lower rank, $\widehat K$ is driven
below $K$, the residual retains a factor component and $\hat\sigma$ inflates;
size is largely preserved by the self-normalisation but power is lost. Comparing
`fit$K_hat` against `fit$eigenvalues` is a quick diagnostic.

## Simulated data

`generate_data` reproduces the design of the simulation study and produced the
two shipped data sets. Note that `mu` has length $K + p - 1$ with $K = 3$: it is
the mean of the factor loadings, not of the coefficients.

```r
sim <- generate_data(obs_N = 20, obs_T = 20, p = 61, mu = rep(1, 63), RHO = 0.25)
sim$data_estimation   # $x (nT x p), $y (nT)
sim$data_inference    # $y has one column per signal c** = 0, 0.1, 0.2
```

## Citation

```r
citation("hdcce")
```

## License

GPL (>= 2)