#' Inference for HD panels with IFE
#'
#' \code{hdcce_inference} High-dimensional inference with the desparsified CCE
#' estimator proposed in Linton, 0., Ruecker, M., Vogt, M., Walsh, C. (2024)
#' "Estimation and Inference in High-Dimensional Panel Data Models with Interactive
#' Fixed Effects".
#'
#' @param data List containing the balanced panel data with data$y containing
#'   the dependent variables and data$x the regressors. Both are sorted such
#'   that first T observations are those for unit 1, followed by the T
#'   observations for unit 2, etc.
#' @param obs_N,obs_T The number of cross-section units (obs_N) and (obs_T) in
#'   the data.
#' @param TRUNC The truncation parameter tau used to estimate the number of
#'   factors. Default is 0.01.
#' @param NFACTORS Allows to set the number of factors use in estimation.
#'   Default is NULL so that the data driven choice with TRUNC = 0.01 is used.
#' @param NFOLDS Specify the number of folds to use when using cross-validation
#'   to select lasso penalty based on cv.glmnet. Default is the default in
#'   cv.glmnet, which is 10.
#'@param foldid Integer vector (obs_N*obs_T - dimensional) containing a label
#'   for each sample to determine its fold for CV.
#'@param COEF_INDEX Index of regressor to be tested.
#'@param alpha Vector containing significance levels.
#' @return The function returns the following:
#' \itemize{
#'  \describe{
#'   \item{$coef_despar}{Desparsified lasso estimate.}
#'   \item{$Avar}{Asymptotic variance of desparsified lasso.}
#'   \item{$conf_bands}{Symmetric confidence bands for given alpha.}
#'   \item{$var_est_Lasso}{Estimated variance of idiosyncratic regression error.}
#'   }
#' }
#' @examples
#' # Load the data set
#' data("data_example")
#'
#' # Set the dimensions of the data
#' obs_N <- 50
#' obs_T <- 50
#'
#' # Do the estimation
#' estimate_model <- hdcce_inference(data, obs_N, obs_T, TRUNC = 0.01, NFACTORS = NULL,
#'  NFOLDS = NULL, foldid = NULL, COEF_INDEX, alpha = c(0.01, 0.05, 0.1))
#'
#' @references  Linton, 0., Ruecker, M., Vogt, M., Walsh, C. (2024) "Estimation
#' and Inference in High-Dimensional Panel Data Models with Interactive
#' Fixed Effects"




#' @export
hdcce_inference <- function(data, obs_N, obs_T, TRUNC = 0.01, NFACTORS = NULL,
                            NFOLDS = NULL, foldid = NULL, COEF_INDEX,
                            alpha = c(0.01, 0.05, 0.1)){



# Initial Checks
#--------------------------------------------------------------------------

  # Check to see whether a correct variant has been supplied
  if ((variant != "LS") * (variant != "Lasso") == 1){
    stop('Variant must be set to "LS" or "Lasso"')
  }
  # Interception for foldid
  if(class(foldid) == "numeric"){
    if(all(foldid == floor(foldid)) == FALSE){
      stop('Provided vector for CV must contain integers only')
    }
    if(length(foldid) != (obs_N * obs_T)){
      stop('Provided vector for CV has wrong dimension.')
    }
  }
  # Interception for trunc
  if(is.null(TRUNC) == FALSE){
    if(TRUNC > 1 | TRUNC <= 0){
      stop("Supplied truncation invalid. Must be in (0,1]")
    }
  }

  # Interception for NFOLDS
  if(is.null(NFOLDS) == FALSE){
    if(NFOLDS >= obs_N){
      stop("Supplied number NFOLDS must be less then obs_N")
    }
    if(NFOLDS != floor(NFOLDS)){
      stop("Supplied number NFOLDS must be integer valued")
    }
  }


  # Pull out the data
  X_data <- data$x
  Y_data <- data$y
  p <- dim(X_data)[2]

  # Interception for the data
  if(length(X_data[,1]) != obs_N * obs_T){
    stop("Supplied dimensions differ.")
  }
  if(length(Y_data) != obs_N * obs_T){
    stop("Supplied dimensions differ.")
  }

  # Interception for NFACTORS
  if(is.null(NFACTORS) == FALSE){
    if(NFACTORS >= p){
      stop("Supplied number NFACTORS must be less then p")
    }
    if(NFACTORS != floor(NFACTORS)){
      stop("Supplied number NFACTORS must be integer valued")
    }
  }

  #============================================================================#
  # Step 1 a):  Calculate the projection matrix needed for inference
  #============================================================================#

  # Cross-sectional averages of the regressors
  X_bar <- matrix(NA, ncol = p, nrow = obs_T)
  for(t in 1:obs_T){
    indices <- seq(t, obs_N * obs_T, by = obs_T)
    X_bar[t,] <- colMeans(X_data[indices,])
  }

  # Empirical covariance matrix and eigenstructure without the "j-th" regressor
  Cov_X_bar <- 1/obs_T * t(X_bar[,-COEF_INDEX]) %*% X_bar[,-COEF_INDEX]
  Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric = TRUE)


  #============================================================================#
  # Step 1 b):  Estimation of number of factors
  #============================================================================#

  # Normalize the eigenvalues by the largest one
  eigen_values <- Cov_X_bar_eigen$values/Cov_X_bar_eigen$values[1]

  if(is.null(NFACTORS) == TRUE){
    K_hat <- NFACTORS
  }
  ## Number of normalized eigenvalues larger than TRUNC
  else{
    K_hat <- sum((TRUNC < eigen_values))
    message(paste("Number of factors estimated given by 'K_hat' =", K_hat))
  }

  #============================================================================#
  # Step 1 (c): Computation of projection matrix
  #============================================================================#

  W_tilde <- X_bar[,-COEF_INDEX] %*% Cov_X_bar_eigen$vectors[,1:K_hat]
  Pi_tilde <- diag(1, obs_T, obs_T) -  W_tilde %*% solve(t(W_tilde) %*% W_tilde)  %*% t(W_tilde)


  #============================================================================#
  # Step 2 a): Transform the data
  #============================================================================#

  Y_tilde <- rep(NA, obs_T * obs_N)
  X_tilde <- matrix(NA, nrow = obs_N * obs_T, ncol = p)

  for(i in 1:obs_N){
    indices <- ((i-1) * obs_T + 1) : (i * obs_T)
    Y_tilde[indices] <- Pi_tilde %*% t(t(Y_data[indices]))
    X_tilde[indices,] <- Pi_tilde %*% t(t(X_data[indices,]))
  }


  #============================================================================#
  # Step 2 b): Run the MAIN Lasso regression on the transformed data and get variance
  #            estimate
  #============================================================================#

    # Check for user-specified number of folds
    if(is.null(NFOLDS) == FALSE & is.null(foldid) == TRUE){
      message(paste("User-supplied number of folds given by 'NFOLDS' = ",
                    NFOLDS," is used to create fold vector."))
      if(obs_N %% NFOLDS > 0){ # If not divisible fill first folds with one leftover
        foldid <- c(rep(1:(obs_N %% NFOLDS), each = ((floor(obs_N/NFOLDS)+1) * obs_T)),
                    rep((obs_N %% NFOLDS + 1):NFOLDS, each = floor(obs_N/NFOLDS)*obs_T))
      } else{
        foldid <- rep(1:NFOLDS, each =  (obs_N/NFOLDS * obs_T))
      }

      fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde, nfolds = NFOLDS,
                                        foldid = foldid, family = "gaussian",
                                        alpha = 1, intercept = FALSE)

      coefs_Lasso <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
      Lambda_CV  <- fit_Lasso.cv$lambda.min


    # Check if there is user-specified set of folds
  } else if(is.null(NFOLDS) == FALSE & is.null(foldid) == FALSE){
    fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde, nfolds = NFOLDS,
                                      foldid = foldid, family = "gaussian",
                                      alpha = 1, intercept = FALSE)
    coefs_Lasso <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
    Lambda_CV  <- fit_Lasso.cv$lambda.min
    message("User specified folds for CV selected.")

  } else{ # Neither specified NFOLDS nor foldid
    if(obs_N %% 10 > 0){ # If not divisible fill first folds with one leftover
      foldid <- c(rep(1:(obs_N %% 10), each = ((floor(obs_N/10)+1) * obs_T)),
                  rep((obs_N %% 10 + 1):10, each = floor(obs_N/10)*obs_T))
    } else{
      foldid <- rep(1:10, each =  (obs_N/10 * obs_T))
    }

    fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde, family = "gaussian",
                                      foldid = foldid, alpha = 1,
                                      intercept = FALSE)
    coefs_Lasso <- stats::coef(fit_Lasso.cv, s = "lambda.min")[-1]
    Lambda_CV  <- fit_Lasso.cv$lambda.min
    message("10 Fold CV as explained in paper.")
  }



  # Variance estimate using cross-validated penalty
  yhat_Lasso <- stats::predict(fit_Lasso, newx = X_tilde,
                        type = "response", s = "lambda.min")
  resid_Lasso <- Y_tilde - yhat_Lasso
  var_est_Lasso <- obs_T/(obs_T - K_hat) * mean(resid_Lasso^2)

  #============================================================================#
  # Step 2(b): Run the nodewise Lasso estimation on transformed data
  #============================================================================#

  # Check for user-specified number of folds
  if(is.null(NFOLDS) == FALSE & is.null(foldid) == TRUE){
    message(paste("User-supplied number of folds given by 'NFOLDS' = ",
                  NFOLDS," is used to create fold vector."))
    if(obs_N %% NFOLDS > 0){ # If not divisible fill first folds with one leftover
      foldid <- c(rep(1:(obs_N %% NFOLDS), each = ((floor(obs_N/NFOLDS)+1) * obs_T)),
                  rep((obs_N %% NFOLDS + 1):NFOLDS, each = floor(obs_N/NFOLDS)*obs_T))
    } else{
      foldid <- rep(1:NFOLDS, each =  (obs_N/NFOLDS * obs_T))
    }
    fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX], nfolds = NFOLDS,
                                   foldid = foldid, family = "gaussian",
                                   alpha = 1, intercept = FALSE)

    coefs_node_Lasso <- stats::coef(fit_node_Lasso, s = "lambda.min")[-1]
    kappa_cv  <- fit_node_Lasso$lambda.min


  # Check if there is user-specified set of folds
  } else if(is.null(NFOLDS) == FALSE & is.null(foldid) == FALSE){
    fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX], nfolds = NFOLDS,
                                   foldid = foldid, family = "gaussian",
                                   alpha = 1, intercept = FALSE)
    coefs_node_Lasso <- stats::coef(fit_node_Lasso, s = "lambda.min")[-1]
    kappa_cv  <- fit_node_Lasso$lambda.min
    message("User specified folds for CV selected.")

  } else{ # Neither specified NFOLDS nor foldid
    if(obs_N %% 10 > 0){ # If not divisible fill first folds with one leftover
      foldid <- c(rep(1:(obs_N %% 10), each = ((floor(obs_N/10)+1) * obs_T)),
                  rep((obs_N %% 10 + 1):10, each = floor(obs_N/10)*obs_T))
    } else{
      foldid <- rep(1:10, each =  (obs_N/10 * obs_T))
    }

    fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX], family = "gaussian",
                                   foldid = foldid, alpha = 1,
                                   intercept = FALSE)
    coefs_node_Lasso <- stats::coef(fit_node_Lasso, s = "lambda.min")[-1]
    kappa_cv  <- fit_node_Lasso$lambda.min
    message("10 Fold CV as explained in paper.")
  }


  # Kappa grid that cv.glmnet uses for cross-validation
  kappa_grid <- fit_node_Lasso$lambda

  # Index for the CV chosen kappa
  kappa_cv_idx <- fit_node_Lasso$index[1]

  #============================================================================#
  # Step 3: Calculate the de-sparsified estimator
  #============================================================================#

  #============================================================================#
  # Step 3(a): Choose the nodewise penalty parameter
  #============================================================================#
  # Kappa grid length
  kappa_grid_len <- length(kappa_grid)

  # Calculate the "scaled asymptotic variance"
  var_scaled <- rep(0, len = kappa_grid_len)
  for(i in 1:kappa_grid_len){
    # (i) Nodewise lasso residuals
    yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                               type = "response", s = kappa_grid[i])
    resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso

    # (ii) "Scaled asymptotic variance
    var_scaled[i] <-  t(resid_node_Lasso) %*% resid_node_Lasso / (t(X_tilde[,COEF_INDEX]) %*% resid_node_Lasso)^2
  }


  # 25% increase of the scaled variance for truncation
  V_trunc <- 1.25 * var_scaled[kappa_cv_idx]

  # Indicator whether the scaled variance is less than the threshold
  trunc_ind <- (var_scaled < V_trunc)
  kappa_idx <- sum(trunc_ind) # CARE IS VAR SCALED REALLY STRICTLY INCREASING?
  KAPPA <- fit_node_Lasso.cv$lambda[kappa_idx]


  #============================================================================#
  # Step 3(b): Construction of the de-sparsified estimator
  #============================================================================#
  yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                             type = "response", s = kappa_grid[kappa_idx])
  resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso

  despar_beta <- coefs_Lasso[COEF_INDEX] + t(resid_node_Lasso) %*% resid_Lasso / t(resid_node_Lasso) %*% X_tilde[, COEF_INDEX]

  # Collect results
  Avar <- sqrt(var_est_Lasso * var_scaled[kappa_idx])

  conf_band_min <- rep(despar_beta, length(alpha)) + Avar * qnorm(alpha/2)
  conf_band_max <- rep(despar_beta, length(alpha)) + Avar * qnorm(1-alpha/2)

  confidence_band <- cbind(conf_band_min, conf_band_max)


  return(list(coef_despar = despar_beta, Avar = Avar,
              confidence_band = confidence_band, var_est = var_est_Lasso))
}


