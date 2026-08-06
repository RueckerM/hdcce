#' Estimate HD panels with IFE
#'
#' \code{hdcce_estimator} fits (i) the linear HD-CCE estimator of Rücker et al.
#' (2025), or (ii) the HD-CCE estimator for a parametric additive (dictionary)
#' design as in Rücker et al. (2026), when \code{dictionaries} is supplied.
#'
#' @param data List containing the balanced panel data with \code{data$y} the
#'   dependent variable and \code{data$x} the p regressors. Both are sorted such
#'   that the first \code{obs_T} observations are those for unit 1, followed by
#'   the \code{obs_T} observations for unit 2, etc.
#' @param obs_N,obs_T The number of cross-section units (\code{obs_N}) and the
#'   time series length (\code{obs_T}).
#' @param dictionaries Either \code{NULL} (linear model, the default), an
#'   \code{obs_N * obs_T} x \code{d} matrix of dictionary transformations of
#'   \code{data$x}, or a list with components \code{Phi} (that matrix) and,
#'   optionally, \code{group}, a \code{d}-vector with entries in
#'   \code{1,...,p} indicating which covariate each column transforms. When a
#'   dictionary is supplied, the design and the response are empirically
#'   centred across units before projection.
#' @param TRUNC The truncation parameter tau used to estimate the number of
#'   factors. Default is 0.01.
#' @param NFACTORS Allows to set the number of factors used in estimation.
#'   Default is NULL, so that the data driven choice with TRUNC is used.
#' @param variant Either \code{"Lasso"} (default) or \code{"LS"} for the least
#'   squares variant. \code{"LS"} requires the design to have full column rank.
#' @param lambda User specified lambda grid, on the scale of the penalised
#'   criterion in the paper (see Details).
#' @param NFOLDS The number of folds (partitioned along the cross-section) used
#'   for cross-validation. Default is 10. Fold size can vary by one if
#'   \code{obs_N} is not divisible by \code{NFOLDS}.
#' @param foldid Integer vector (\code{obs_N * obs_T}-dimensional) containing a
#'   fold label for each observation.
#' @param scree_plot Logical; whether a scree plot of the eigendecomposition of
#'   \eqn{\Sigma} should be shown. The default is FALSE.
#' @param standardize Logical; whether \code{glmnet} is called with the
#'   standardized projected data. The default is TRUE.
#' @details The criterion in the paper is
#'   \eqn{(nT)^{-1} \sum_i \| \tilde{Y}_i - \tilde{X}_i b \|^2 + \lambda \|b\|_1},
#'   whereas \code{glmnet} minimises
#'   \eqn{(2nT)^{-1} \sum_i \| \tilde{Y}_i - \tilde{X}_i b \|^2 + \lambda_g \|b\|_1}.
#'   The two penalties are related by \eqn{\lambda = 2 \lambda_g}. Both the
#'   \code{lambda} argument and the returned \code{Lambda} are on the scale of
#'   the paper.
#' @return A list with
#' \describe{
#'   \item{coefs}{The coefficient estimates. A numeric vector, except when a
#'     \code{lambda} grid is supplied, in which case a \code{d} x
#'     \code{length(lambda)} matrix with one column per penalty.}
#'   \item{K_hat}{The estimated number of factors, or \code{NFACTORS} if set.}
#'   \item{Lambda}{The penalty selected by cross-validation, on the scale of the
#'     paper; the supplied grid if \code{lambda} was given; \code{NULL} for
#'     \code{variant = "LS"}.}
#'   \item{eigenvalues}{The normalised eigenvalues used to determine
#'     \code{K_hat}.}
#'   \item{type}{Either \code{"linear"} or \code{"dictionary"}.}
#'   \item{call}{The matched call.}
#' }
#' @references Rücker, M., Vogt, M., Linton, O. and Walsh, C. (2025).
#'   Estimation and inference in high-dimensional panel data models with
#'   interactive fixed effects. \emph{Quantitative Economics}, 16(4),
#'   1457--1509. \doi{10.3982/QE2308}
#'
#'   Rücker, M., Vogt, M. and Linton, O. (2026). High-dimensional panel data
#'   models with interactive fixed effects: beyond the linear case.
#'   \url{https://arxiv.org/abs/2608.02055}
#' @examples
#' data("data_estimation")
#' obs_N <- 20
#' obs_T <- 20
#' fit <- hdcce_estimator(data = data_estimation, obs_N = obs_N, obs_T = obs_T)
#' fit$K_hat
#' head(fit$coefs, 10)
#' @export
hdcce_estimator <- function(data, obs_N, obs_T, dictionaries = NULL,
                            TRUNC = 0.01, NFACTORS = NULL,
                            variant = c("Lasso", "LS"),
                            lambda = NULL, NFOLDS = 10, foldid = NULL,
                            scree_plot = FALSE, standardize = TRUE){

  variant <- match.arg(variant)
  cl      <- match.call()

  #==========================================================================#
  # Initial checks
  #==========================================================================#
  if(is.numeric(foldid)){
    if(!all(foldid == floor(foldid)))
      stop("Provided vector for CV must contain integers only.")
    if(length(foldid) != (obs_N * obs_T))
      stop("Provided vector for CV has wrong dimension.")
  }
  if(length(TRUNC) != 1 || !is.numeric(TRUNC) || TRUNC > 1 || TRUNC <= 0)
    stop("TRUNC must be a single number in (0,1].")
  if(NFOLDS != floor(NFOLDS))
    stop("Supplied number NFOLDS must be integer valued.")
  if(NFOLDS >= obs_N)
    stop("Supplied number NFOLDS must be less than obs_N.")
  if(!is.null(lambda) && any(lambda <= 0))
    stop("Supplied lambda grid must be positive.")

  X_data <- as.matrix(data$x)
  Y_data <- as.vector(data$y)
  p      <- ncol(X_data)
  n_obs  <- obs_N * obs_T

  if(nrow(X_data) != n_obs || length(Y_data) != n_obs)
    stop("Supplied dimensions differ.")
  if(anyNA(X_data) || anyNA(Y_data))
    stop("Data must not contain missing values.")
  if(!is.null(NFACTORS)){
    if(NFACTORS != floor(NFACTORS))
      stop("Supplied number NFACTORS must be integer valued.")
    if(NFACTORS < 1 || NFACTORS >= min(p, obs_T))
      stop("NFACTORS must be an integer in [1, min(p, obs_T)).")
  }

  #--- dictionary set-up ----------------------------------------------------#
  use_dict <- !is.null(dictionaries)
  if(use_dict){
    if(is.list(dictionaries)){
      Phi <- as.matrix(dictionaries$Phi)
      grp <- dictionaries$group
    } else {
      Phi <- as.matrix(dictionaries)
      grp <- NULL
    }
    if(is.null(Phi))
      stop("dictionaries must contain a component 'Phi'.")
    if(nrow(Phi) != n_obs)
      stop("dictionaries$Phi must have obs_N * obs_T rows.")
    if(anyNA(Phi))
      stop("dictionaries$Phi must not contain missing values.")
    if(!is.null(grp)){
      if(length(grp) != ncol(Phi))
        stop("dictionaries$group must have length ncol(dictionaries$Phi).")
      if(!all(grp %in% seq_len(p)))
        stop("dictionaries$group must contain entries in 1,...,p.")
    }
  }

  #==========================================================================#
  # Step 1 a): eigendecomposition of the empirical covariance matrix
  #==========================================================================#
  X_bar <- matrix(NA_real_, nrow = obs_T, ncol = p)
  for(t in 1:obs_T){
    X_bar[t, ] <- colMeans(X_data[seq(t, n_obs, by = obs_T), , drop = FALSE])
  }

  Cov_X_bar_eigen <- eigen(crossprod(X_bar)/obs_T, symmetric = TRUE)
  eigen_values    <- Cov_X_bar_eigen$values / Cov_X_bar_eigen$values[1]

  #==========================================================================#
  # Step 1 b): number of factors
  #==========================================================================#
  if(is.numeric(NFACTORS)){
    K_hat <- NFACTORS
    message(paste("User-supplied number of factors given by 'NFACTORS' = ",
                  NFACTORS, " is used in estimation."))
    if(isTRUE(scree_plot)){
      graphics::plot(eigen_values, ylim = c(0, 1),
                     ylab = "Normalized Eigenvalues",
                     main = paste("Number of factors set to ", K_hat))
      graphics::points(eigen_values[1:K_hat], col = "red")
    }
  } else {
    K_hat <- sum(TRUNC < eigen_values)
    message(paste("Number of factors estimated given by 'K_hat' =", K_hat))
    if(isTRUE(scree_plot)){
      graphics::plot(eigen_values, ylim = c(0, 1),
                     ylab = "Normalized Eigenvalues",
                     main = paste("Estimated number of factors = ", K_hat))
      graphics::abline(h = TRUNC, col = "red")
      graphics::legend("topright", legend = c("Truncation"), lty = 1, col = "red")
    }
  }
  if(K_hat < 1)
    stop("No eigenvalue exceeds TRUNC; the estimated number of factors is zero.")
  if(K_hat >= obs_T)
    stop("Estimated number of factors exceeds obs_T - 1.")

  #==========================================================================#
  # Step 1 c): projection matrix, via an orthonormal basis of W_hat
  #==========================================================================#
  W_hat  <- X_bar %*% Cov_X_bar_eigen$vectors[, seq_len(K_hat), drop = FALSE]
  Q_hat  <- qr.Q(qr(W_hat))
  Pi_hat <- diag(1, obs_T, obs_T) - tcrossprod(Q_hat)

  #==========================================================================#
  # Step 1 d): centre (dictionary case only) and project
  #==========================================================================#
  if(use_dict){
    Design_raw <- .hdcce_center_units(Phi, obs_N, obs_T)
    Y_raw      <- .hdcce_center_units(Y_data, obs_N, obs_T)
  } else {
    Design_raw <- X_data
    Y_raw      <- Y_data
  }

  X_hat <- .hdcce_project_units(Design_raw, Pi_hat, obs_N, obs_T)
  Y_hat <- as.vector(.hdcce_project_units(Y_raw, Pi_hat, obs_N, obs_T))

  coef_names <- if(use_dict && !is.null(grp)){
    paste0("x", grp, ".", stats::ave(grp, grp, FUN = seq_along))
  } else if(use_dict){
    paste0("phi", seq_len(ncol(X_hat)))
  } else {
    cn <- colnames(X_data)
    if(!is.null(cn) && all(nzchar(cn)) && !anyDuplicated(cn)) cn else NULL
  }

  #==========================================================================#
  # Step 2: estimation on the transformed data
  #==========================================================================#
  Lambda <- NULL

  if(variant == "LS"){
    if(ncol(X_hat) >= nrow(X_hat))
      stop("The least squares variant requires fewer columns than observations; ",
           "use variant = 'Lasso'.")
    coefs <- stats::lm.fit(x = X_hat, y = Y_hat)$coefficients
    if(anyNA(coefs))
      warning("The projected design is rank deficient; ",
              "some least squares coefficients are not identified.")
    message("LS variant selected.")

  } else {

    if(!is.null(lambda)){
      # glmnet penalty is half the penalty of the paper
      fit_Lasso <- glmnet::glmnet(x = X_hat, y = Y_hat, family = "gaussian",
                                  alpha = 1, lambda = lambda/2,
                                  intercept = FALSE, standardize = standardize)
      coefs  <- as.matrix(fit_Lasso$beta)
      Lambda <- 2 * fit_Lasso$lambda
      message("User specified lambda grid selected.")

    } else {

      if(is.null(foldid)){
        message(paste("User-supplied number of folds given by 'NFOLDS' = ",
                      NFOLDS, " is used to create fold vector."))
        if(obs_N %% NFOLDS > 0 & obs_N > NFOLDS){
          fold_vec <- c(rep(1:(obs_N %% NFOLDS),
                            each = ((floor(obs_N/NFOLDS) + 1) * obs_T)),
                        rep((obs_N %% NFOLDS + 1):NFOLDS,
                            each = floor(obs_N/NFOLDS) * obs_T))
          message("Fold assignment based on NFOLDS with modulo.")
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

      fit_Lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = fold_vec,
                                     family = "gaussian", alpha = 1,
                                     intercept = FALSE, standardize = standardize)
      coefs  <- as.vector(stats::coef(fit_Lasso, s = "lambda.min"))[-1]
      Lambda <- 2 * fit_Lasso$lambda.min
    }
  }

  if(!is.null(coef_names) && is.null(dim(coefs)) && length(coefs) == length(coef_names))
    names(coefs) <- coef_names
  if(!is.null(coef_names) && !is.null(dim(coefs)) && nrow(coefs) == length(coef_names))
    rownames(coefs) <- coef_names

  out <- list(coefs       = coefs,
              K_hat       = K_hat,
              Lambda      = Lambda,
              eigenvalues = eigen_values,
              type        = if(use_dict) "dictionary" else "linear",
              call        = cl)
  class(out) <- "hdcce_fit"
  out
}
