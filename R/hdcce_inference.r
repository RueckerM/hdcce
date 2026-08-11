#' Inference for HD panels with IFE
#'
#' \code{hdcce_inference} performs high-dimensional inference with either
#' (i) the desparsified HD-CCE estimator of Rücker et al. (2025), giving
#' confidence intervals for individual coefficients, or (ii) the max-type
#' significance test of Rücker et al. (2026) for the null that a covariate has
#' no effect in the parametric additive model, when a dictionary is supplied.
#'
#' @param data List containing the balanced panel data with \code{data$y} the
#'   dependent variable and \code{data$x} the regressors. Both are sorted such
#'   that the first \code{obs_T} observations are those for unit 1, followed by
#'   the \code{obs_T} observations for unit 2, etc.
#' @param obs_N,obs_T The number of cross-section units (\code{obs_N}) and the
#'   time series length (\code{obs_T}).
#' @param COEF_INDEX_VEC Indices of the regressors to be tested.
#' @param dictionaries Either \code{NULL} (linear model, the default) or a list
#'   with components \code{Phi}, an \code{obs_N * obs_T} x \code{d} matrix of
#'   dictionary transformations of \code{data$x}, and \code{group}, a
#'   \code{d}-vector with entries in \code{1,...,p} indicating which covariate
#'   each column transforms. Optionally \code{Phi_nw} and \code{group_nw} give a
#'   separate dictionary for the nodewise regression; they default to
#'   \code{Phi} and \code{group}.
#' @param TRUNC The truncation parameter tau used to estimate the number of
#'   factors. Default is 0.01.
#' @param NFACTORS Allows to set the number of factors used in estimation.
#'   Default is NULL, so that the data driven choice with TRUNC is used.
#' @param NFOLDS Number of folds for cross-validation in \code{cv.glmnet}.
#' @param foldid Integer vector (\code{obs_N * obs_T}-dimensional) with a fold
#'   label for each observation.
#' @param alpha Vector of significance levels.
#' @param HAC Integer 1, 2 or 3, selecting the homoscedastic, heteroscedasticity
#'   robust, or HAC variance estimator. Only used in the linear case.
#' @param standardize Logical; whether \code{glmnet} is called with standardized
#'   projected data.
#' @param s_rule Either \code{"lambda.1se"} (default) or \code{"lambda.min"},
#'   the rule used to select the penalties from cross-validation.
#' @param kernel Kernel with support \code{[-1,1]} used for the test, either
#'   \code{"epanechnikov"}, \code{"biweight"}, or a function.
#' @param C Half-width of the region of locations scanned by the test, in the
#'   units of the nodewise residuals: the statistic is maximised over
#'   \eqn{w \in [-C, C]}. The default is \code{NULL}, in which case \code{C} is set to the
#'   \code{QUANT} quantile of the absolute nodewise residuals, so that the
#'   interval \eqn{[-C, C]} contains that share of them.
#' @param QUANT Quantile used for the data-driven choice of \code{C}; ignored
#'   when \code{C} is supplied. The default is 0.9. Larger values reach further
#'   into the tails at the cost of thinner outer bumps.
#' @param h Bandwidth of the kernel weights, in the same units as \code{C}. The
#'   default is \code{NULL}, in which case \code{h} is chosen by Silverman's
#'   rule of thumb applied to the nodewise residuals,
#'   \eqn{0.9 \hat\sigma_u (nT)^{-1/5}}.
#' @param B Number of Monte Carlo draws for the Gaussian coupling.
#' @return A list with components \code{results} (a list with one entry per
#'   element of \code{COEF_INDEX_VEC}, named by the index), \code{alpha},
#'   \code{obs_N}, \code{obs_T}, \code{type} and \code{call}. Each entry of
#'   \code{results} is a list with
#' \describe{
#'   \item{coef_despar}{Desparsified estimate; \code{NULL} in the dictionary case.}
#'   \item{se}{Estimated standard error; \code{NULL} in the dictionary case.}
#'   \item{statistic}{z-statistic for H0: beta_j = 0, or the max-type statistic in the dictionary case..}
#'   \item{critical_values}{Critical value for each level in \code{alpha}.}
#'   \item{p_value}{Two-sided normal p-value, or Monte Carlo p-value in the dictionary case..}
#'   \item{reject}{Logical, one entry per level in \code{alpha}.}
#'   \item{confidence_band}{Confidence interval per level; \code{NULL} in the
#'     dictionary case.}
#'   \item{profile}{Data frame with the grid of locations, the corresponding
#'     Psi_w and the number of residuals in each bump; \code{NULL} in the
#'     linear case.}
#'   \item{K_hat, sigma_hat}{Estimated number of factors and error standard
#'     deviation.}
#' }
#' @details
#' The grid of locations is
#' \eqn{\mathcal{W} = \{ w \in [-C, C] : w = -C + (2\ell - 1) h,\ \ell = 1, 2,
#' \ldots \}}, so that the supports of the kernel weights are disjoint and the
#' number of locations is \eqn{L = \lfloor C / h + 1/2 \rfloor}
#'
#' @references Rücker, M., Vogt, M., Linton, O. and Walsh, C. (2025).
#'   Estimation and inference in high-dimensional panel data models with
#'   interactive fixed effects. \emph{Quantitative Economics}, 16(4),
#'   1457--1509. \doi{10.3982/QE2308}
#'
#'   Rücker, M., Vogt, M. and Linton, O. (2026). High-dimensional panel data
#'   models with interactive fixed effects: beyond the linear case.
#'   \url{https://arxiv.org/abs/2608.02055}
#' @examples
#' data("data_inference")
#' obs_N <- 20
#' obs_T <- 20
#' dat <- list(x = data_inference$x, y = data_inference$y[, 1])
#' fit <- hdcce_inference(dat, obs_N, obs_T, COEF_INDEX_VEC = 1)
#' fit$results[["1"]]$confidence_band
#' @export
hdcce_inference <- function(data, obs_N, obs_T, COEF_INDEX_VEC,
                            dictionaries = NULL,
                            TRUNC = 0.01, NFACTORS = NULL,
                            NFOLDS = 10, foldid = NULL,
                            alpha = c(0.01, 0.05, 0.1), HAC = 2,
                            standardize = TRUE,
                            s_rule = "lambda.1se",
                            kernel = "epanechnikov", C = NULL, QUANT = 0.9,
                            h = NULL,
                            B = 5000){

  s_rule <- match.arg(s_rule)
  cl     <- match.call()

  #==========================================================================#
  # Initial checks
  #==========================================================================#
  if(is.numeric(foldid)){
    if(!all(foldid == floor(foldid)))
      stop("Provided vector for CV must contain integers only.")
    if(length(foldid) != (obs_N * obs_T))
      stop("Provided vector for CV has wrong dimension.")
  }
  if(TRUNC > 1 | TRUNC <= 0)
    stop("Supplied truncation must be in (0,1].")
  if(NFOLDS >= obs_N)
    stop("Supplied number NFOLDS must be less than obs_N.")
  if(NFOLDS != floor(NFOLDS))
    stop("Supplied number NFOLDS must be integer valued.")
  if(length(HAC) != 1 || !(HAC %in% c(1, 2, 3)))
    stop("HAC must be equal to 1, 2 or 3.")
  if(any(alpha <= 0) || any(alpha >= 1))
    stop("All significance levels in alpha must lie in (0,1).")

  X_data <- as.matrix(data$x)
  Y_data <- as.vector(data$y)
  p      <- ncol(X_data)
  n_obs  <- obs_N * obs_T

  if(nrow(X_data) != n_obs || length(Y_data) != n_obs)
    stop("Invalid data dimensions.")
  if(!is.null(NFACTORS)){
    if(NFACTORS >= p)      stop("Supplied number NFACTORS must be less than p.")
    if(NFACTORS >= obs_T)  stop("Supplied number NFACTORS must be less than obs_T.")
    if(NFACTORS != floor(NFACTORS))
      stop("Supplied number NFACTORS must be integer valued.")
  }
  if(any(COEF_INDEX_VEC > p | COEF_INDEX_VEC < 1 |
         floor(COEF_INDEX_VEC) != COEF_INDEX_VEC))
    stop("COEF_INDEX_VEC must contain integers between 1 and p.")

  #--- dictionary set-up ----------------------------------------------------#
  use_dict <- !is.null(dictionaries)
  if(use_dict){
    Phi <- as.matrix(dictionaries$Phi)
    grp <- dictionaries$group
    if(nrow(Phi) != n_obs)
      stop("dictionaries$Phi has wrong number of rows.")
    if(length(grp) != ncol(Phi))
      stop("dictionaries$group must have length ncol(dictionaries$Phi).")
    if(!all(grp %in% seq_len(p)))
      stop("dictionaries$group must contain entries in 1,...,p.")
    Phi_nw <- if(is.null(dictionaries$Phi_nw)) Phi else as.matrix(dictionaries$Phi_nw)
    grp_nw <- if(is.null(dictionaries$group_nw)) grp else dictionaries$group_nw
    if(length(grp_nw) != ncol(Phi_nw))
      stop("dictionaries$group_nw must have length ncol(dictionaries$Phi_nw).")
    for(j in COEF_INDEX_VEC){
      if(!any(grp == j))
        stop("No dictionary column belongs to covariate ", j, ".")
    }
  } else {
    Phi <- NULL
    grp <- seq_len(p)                  # makes the linear case a special case
  }

  #==========================================================================#
  # Fold assignment (computed once; never overwrites the user's foldid)
  #==========================================================================#
  if(is.null(foldid)){
    if(obs_N %% NFOLDS > 0 & obs_N > NFOLDS){
      fold_vec <- c(rep(1:(obs_N %% NFOLDS),
                        each = ((floor(obs_N/NFOLDS) + 1) * obs_T)),
                    rep((obs_N %% NFOLDS + 1):NFOLDS,
                        each = floor(obs_N/NFOLDS) * obs_T))
      message("Fold assignment based on NFOLDS modulo.")
    } else if(obs_N > NFOLDS){
      fold_vec <- rep(1:NFOLDS, each = (obs_N/NFOLDS * obs_T))
      message("Fold assignment based on NFOLDS.")
    } else {
      fold_vec <- rep(1:obs_N, each = obs_T)
      message("Leave one out cv used.")
    }
  } else {
    fold_vec <- foldid
    message("User specified folds for CV selected.")
  }

  #==========================================================================#
  # Step 1 a): number of factors, from the full cross-sectional averages
  #==========================================================================#
  X_bar <- matrix(NA_real_, nrow = obs_T, ncol = p)
  for(t in 1:obs_T){
    X_bar[t, ] <- colMeans(X_data[seq(t, n_obs, by = obs_T), , drop = FALSE])
  }

  eig_full <- eigen(crossprod(X_bar)/obs_T, symmetric = TRUE)
  eig_norm <- eig_full$values / eig_full$values[1]

  if(!is.null(NFACTORS)){
    K_hat <- NFACTORS
  } else {
    K_hat <- sum(TRUNC < eig_norm)
    message(paste("Number of factors estimated given by 'K_hat' =", K_hat))
  }
  if(K_hat >= obs_T) stop("Estimated number of factors exceeds obs_T - 1.")

  kfun <- .hdcce_kernel(kernel)
  if(is.null(h)) h <- n_obs^(-1/5)

  Results <- vector(mode = "list", length = length(COEF_INDEX_VEC))
  names(Results) <- as.character(COEF_INDEX_VEC)

  #==========================================================================#
  # Loop over the tested indices
  #==========================================================================#
  for(k in seq_along(COEF_INDEX_VEC)){

    COEF_INDEX <- COEF_INDEX_VEC[k]

    #------------------------------------------------------------------------#
    # Step 1 b)-c): projection matrix from X_bar without the j-th column
    #------------------------------------------------------------------------#
    X_bar_mj <- X_bar[, -COEF_INDEX, drop = FALSE]
    eig_mj   <- eigen(crossprod(X_bar_mj)/obs_T, symmetric = TRUE)

    W_tilde  <- X_bar_mj %*% eig_mj$vectors[, seq_len(K_hat), drop = FALSE]
    Q_tilde  <- qr.Q(qr(W_tilde))                      # orthonormal basis
    Pi_tilde <- diag(1, obs_T, obs_T) - tcrossprod(Q_tilde)

    #------------------------------------------------------------------------#
    # Step 1 d): centre (dictionary case only) and project
    #------------------------------------------------------------------------#
    if(use_dict){
      Design_raw <- .hdcce_center_units(Phi,    obs_N, obs_T)
      Nodew_raw  <- .hdcce_center_units(Phi_nw, obs_N, obs_T)
      Y_raw      <- .hdcce_center_units(Y_data, obs_N, obs_T)
      Xj_raw     <- .hdcce_center_units(X_data[, COEF_INDEX, drop = FALSE],
                                        obs_N, obs_T)
    } else {
      Design_raw <- X_data
      Nodew_raw  <- X_data
      Y_raw      <- Y_data
      Xj_raw     <- X_data[, COEF_INDEX, drop = FALSE]
    }

    X_tilde  <- .hdcce_project_units(Design_raw, Pi_tilde, obs_N, obs_T)
    Y_tilde  <- as.vector(.hdcce_project_units(Y_raw, Pi_tilde, obs_N, obs_T))
    Xj_tilde <- as.vector(.hdcce_project_units(Xj_raw, Pi_tilde, obs_N, obs_T))
    Nw_tilde <- if(use_dict){
      .hdcce_project_units(Nodew_raw, Pi_tilde, obs_N, obs_T)
    } else X_tilde

    #------------------------------------------------------------------------#
    # Step 2 a): main lasso
    #------------------------------------------------------------------------#
    fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde, foldid = fold_vec,
                                   family = "gaussian", alpha = 1,
                                   intercept = FALSE, standardize = standardize)
    coefs_Lasso <- as.vector(stats::coef(fit_Lasso, s = s_rule))[-1]
    lambda_cv   <- fit_Lasso[[s_rule]]

    resid_Lasso <- as.vector(Y_tilde -
                               stats::predict(fit_Lasso, newx = X_tilde,
                                              type = "response", s = s_rule))

    #------------------------------------------------------------------------#
    # Step 2 b): nodewise lasso
    #   linear     : response = projected x_j, design = projected x_{-j}
    #   dictionary : response = projected raw x_j, design = projected Phi_nw
    #------------------------------------------------------------------------#
    if(use_dict){
      nw_x    <- Nw_tilde[, (grp_nw != COEF_INDEX), drop = FALSE]
    } else {
      nw_x    <- X_tilde[, (grp != COEF_INDEX), drop = FALSE]
    }

    fit_node <- glmnet::cv.glmnet(x = nw_x, y = Xj_tilde, foldid = fold_vec,
                                  family = "gaussian", alpha = 1,
                                  intercept = FALSE, standardize = standardize)
    kappa_cv <- fit_node[[s_rule]]

    resid_node <- as.vector(Xj_tilde -
                              stats::predict(fit_node, newx = nw_x,
                                             type = "response", s = s_rule))

    #========================================================================#
    # Step 3 a): linear case - desparsified estimate and confidence band
    #========================================================================#
    if(!use_dict){

      denom <- as.numeric(crossprod(Xj_tilde, resid_node))

      if(HAC == 1){
        sigma_eps <- (1/n_obs) * (obs_T/(obs_T - K_hat)) * sum(resid_Lasso^2)
        var_scaled <- sigma_eps * sum(resid_node^2) / denom^2
      }
      if(HAC == 2){
        tmp <- numeric(obs_N)
        for(i in 1:obs_N){
          idx <- ((i - 1) * obs_T + 1):(i * obs_T)
          sigma_eps <- (1/obs_T) * (obs_T/(obs_T - K_hat)) * sum(resid_Lasso[idx]^2)
          tmp[i] <- sigma_eps * sum(resid_node[idx]^2)
        }
        var_scaled <- sum(tmp) / denom^2
      }
      if(HAC == 3){
        W_mat <- matrix(1, nrow = obs_T, ncol = obs_T)
        tmp <- numeric(obs_N)
        for(i in 1:obs_N){
          idx <- ((i - 1) * obs_T + 1):(i * obs_T)
          rp  <- resid_node[idx] * resid_Lasso[idx]
          tmp[i] <- as.numeric(crossprod(rp, W_mat %*% rp))
        }
        var_scaled <- sum(tmp) / denom^2
      }

      despar_beta <- coefs_Lasso[COEF_INDEX] +
        as.numeric(crossprod(resid_node, resid_Lasso)) / denom
      se <- sqrt(var_scaled)

      z_stat <- despar_beta / se
      crit   <- stats::qnorm(1 - alpha/2)
      pval   <- 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)

      confidence_band <- cbind(conf_band_min = despar_beta + se * stats::qnorm(alpha/2),
                               conf_band_max = despar_beta + se * stats::qnorm(1 - alpha/2))
      rownames(confidence_band) <- paste0("alpha=", alpha)

      sigma_hat <- sqrt((1/n_obs) * (obs_T/(obs_T - K_hat)) * sum(resid_Lasso^2))

      Results[[k]] <- list(index           = COEF_INDEX,
                           coef_despar     = despar_beta,
                           se              = se,
                           statistic       = z_stat,
                           critical_values = stats::setNames(crit, paste0("alpha=", alpha)),
                           p_value         = pval,
                           reject          = stats::setNames(abs(z_stat) >= crit,
                                                             paste0("alpha=", alpha)),
                           confidence_band = confidence_band,
                           profile         = NULL,
                           K_hat           = K_hat,
                           sigma_hat       = sigma_hat,
                           lambda          = lambda_cv,
                           kappa           = kappa_cv)

      #========================================================================#
      # Step 3 b): dictionary case - max-type significance test
      #========================================================================#
    } else {

      main_keep <- (grp != COEF_INDEX)

      # Residuals under H0 (j-th block removed) and full residuals for sigma^2
      R_vec <- as.vector(Y_tilde -
                           X_tilde[, main_keep, drop = FALSE] %*% coefs_Lasso[main_keep])
      sigma2_hat <- sum(resid_Lasso^2) / (obs_N * (obs_T - K_hat))

      #  Grid of locations: tau_{w,h}(u) = tau((u-w)/h)/sqrt(h)
      #  and W = { w in [-C,C] : w = -C + (2l-1)h },
      #   with C and h constants in the units of the nodewise residuals.
      if (is.null(C) || is.null(h)) {
        sd_u <- sqrt(mean(resid_node^2))          # only used to set defaults
        if (is.null(C)) C <- stats::quantile(abs(resid_node), QUANT, names = FALSE)
        if (is.null(h)) h <- sd_u * n_obs^(-1/5)
      }
      L <- floor(C / h + 0.5)
      if (L < 1) stop("Empty grid: h must not exceed C.")
      w_grid <- -C + (2 * seq_len(L) - 1) * h

      Psi_w <- numeric(L)
      n_eff <- integer(L)

      D     <- matrix(0, nrow = n_obs, ncol = L)

      for (l in seq_len(L)) {
        arg      <- (resid_node - w_grid[l]) / h
        tau_l    <- kfun(arg) / sqrt(h)
        n_eff[l] <- sum(abs(arg) <= 1)
        Pi_tau   <- as.vector(Pi_tilde %*% matrix(tau_l, nrow = obs_T, ncol = obs_N))
        ss       <- sum(Pi_tau^2)
        if (ss > .Machine$double.eps^0.5) {
          Psi_w[l] <- sum(R_vec * tau_l) / sqrt(ss)
          D[, l]   <- Pi_tau / sqrt(ss / n_obs)
        } else {
          Psi_w[l] <- NA_real_
        }
      }

      active <- which(!is.na(Psi_w))
      if(length(active) == 0L)
        stop("All grid points are degenerate; adjust C or h.")
      if(min(n_eff[active]) < 30)
        warning("Fewer than 30 residuals in the sparsest bump; ",
                "the Gaussian approximation may be unreliable. ",
                "Consider a larger h or a smaller C.")

      Psi <- max(abs(Psi_w[active]))

      # Gaussian coupling: randomness over the errors only, D held fixed
      Z         <- matrix(stats::rnorm(n_obs * B), nrow = n_obs, ncol = B)
      lin       <- crossprod(D[, active, drop = FALSE], Z) / sqrt(n_obs)
      Psi_gauss <- sqrt(sigma2_hat) * apply(abs(lin), 2, max)

      crit <- stats::quantile(Psi_gauss, probs = 1 - alpha, type = 1, names = FALSE)
      pval <- (1 + sum(Psi_gauss >= Psi)) / (B + 1)

      Results[[k]] <- list(index           = COEF_INDEX,
                           coef_despar     = NULL,
                           se              = NULL,
                           statistic       = Psi,
                           critical_values = stats::setNames(crit, paste0("alpha=", alpha)),
                           p_value         = pval,
                           reject          = stats::setNames(Psi >= crit,
                                                             paste0("alpha=", alpha)),
                           confidence_band = NULL,
                           profile         = data.frame(w = w_grid, Psi_w = Psi_w,
                                                        n_eff = n_eff),
                           K_hat           = K_hat,
                           sigma_hat       = sqrt(sigma2_hat),
                           bandwidth       = h,
                           C               = C,
                           lambda          = lambda_cv,
                           kappa           = kappa_cv)
    }
  }

  out <- list(results = Results,
              alpha   = alpha,
              obs_N   = obs_N,
              obs_T   = obs_T,
              type    = if(use_dict) "dictionary" else "linear",
              call    = cl)
  out
}
