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
#' @param obs_N,obs_T The number of cross-section units (obs_N) and time series
#'   length (obs_T) in the data.
#' @param TRUNC The truncation parameter tau used to estimate the number of
#'   factors. Default is 0.01.
#' @param NFACTORS Allows to set the number of factors use in estimation.
#'   Default is NULL so that the data driven choice with TRUNC = 0.01 is used.
#' @param NFOLDS Specifies the number of folds to use when using cross-validation
#'   to select the lasso penalty based on cv.glmnet. Default is 10.
#'   Fold size can vary by one if obs_N is not divisible by NFOLDS.
#'@param foldid Integer vector (obs_N*obs_T - dimensional) containing a label
#'   for each sample to determine its fold for CV.
#'@param COEF_INDEX_VEC Indeces of regressors to be tested.
#'@param alpha Vector containing significance levels.
#' @return The function returns a list for where for each COEF_INDEX_VEC it contains
#' the following:
#' \itemize{
#'  \describe{
#'   \item{$coef_despar}{Desparsified lasso estimate.}
#'   \item{$Avar}{Asymptotic variance of desparsified lasso.}
#'   \item{$confidence_band}{Symmetric confidence interval (lower and upper interval
#'   value) for given alpha.}
#'   }
#' }
#' @examples
#' # Load the data set
#' data("data_inference")
#'
#' # Set the dimensions of the data
#' obs_N <- 50
#' obs_T <- 50
#'
#' # signal = 1,2,3 corresponds to c** = 0, 0.1, 0.2
#' signal <- 1
#'
#' # Use correct column of Y
#' data_inference <- list(x = data_inference$x, y = data_inference$y[,signal])
#'
#' # Regressor to be tested
#' COEF_INDEX_VEC <- 1
#'
#' # Inference example
#' inference_model <- hdcce_inference(data_inference, obs_N, obs_T, TRUNC = 0.01,
#'  NFACTORS = NULL, NFOLDS = 10, foldid = NULL, COEF_INDEX_VEC,
#'  alpha = c(0.01, 0.05, 0.1))
#'
#' print(inference_model$confidence_band)
#'
#' @references  Linton, 0., Ruecker, M., Vogt, M., Walsh, C. (2024) "Estimation
#' and Inference in High-Dimensional Panel Data Models with Interactive
#' Fixed Effects"




#' @export
hdcce_inference <- function(data, obs_N, obs_T, TRUNC = 0.01, NFACTORS = NULL,
                            NFOLDS = 10, foldid = NULL, COEF_INDEX_VEC,
                            alpha = c(0.01, 0.05, 0.1)){



# Initial Checks
#--------------------------------------------------------------------------
# Interception for foldid
  if(class(foldid) == "numeric"){
    if(!all(foldid == floor(foldid))){
      stop('Provided vector for CV must contain integers only')
    }
    if(length(foldid) != (obs_N * obs_T)){
      stop('Provided vector for CV has wrong dimension.')
    }
  }
  # Interception for TRUNC
    if(TRUNC > 1 | TRUNC <= 0){
      stop("Supplied truncation must be in (0,1].")
    }

  # Interception for NFOLDS
    if(NFOLDS >= obs_N){
      stop("Supplied number NFOLDS must be less then obs_N.")
    }
    if(NFOLDS != floor(NFOLDS)){
      stop("Supplied number NFOLDS must be integer valued.")
    }


  # Access the data
  X_data <- data$x
  Y_data <- data$y
  p <- dim(X_data)[2]

  # Interception for the data
  if(length(X_data[,1]) != obs_N * obs_T){
    stop("Invalid data dimensions.")
  }
  if(length(Y_data) != obs_N * obs_T){
    stop("Invalid data dimensions.")
  }

  # Interception for NFACTORS
  if(!is.null(NFACTORS)){
    if(NFACTORS >= p){
      stop("Supplied number NFACTORS must be less then p.")
    }
    if(NFACTORS != floor(NFACTORS)){
      stop("Supplied number NFACTORS must be integer valued.")
    }
  }

  if(length(COEF_INDEX_VEC > 1)){
  Results <- vector(mode = 'list', length = length(COEF_INDEX_VEC))
  }

  #============================================================================#
  #   Estimation of number of factors K_hat (UPDATED: WITH MAIN REGRESSION)
  #============================================================================#
  # Cross-sectional averages of the regressors
  X_bar <- matrix(NA, ncol = p, nrow = obs_T)
  for(t in 1:obs_T){
    indices <- seq(t, obs_N * obs_T, by = obs_T)
    X_bar[t,] <- colMeans(X_data[indices,])
  }

  # Empirical covariance matrix and eigenstructure without the "j-th" regressor
  Cov_X_bar <- 1/obs_T * t(X_bar) %*% X_bar
  Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric = TRUE)

  # Normalize the eigenvalues by the largest one
  eigen_values <- Cov_X_bar_eigen$values/Cov_X_bar_eigen$values[1]

  if(!is.null(NFACTORS)){
    K_hat <- NFACTORS
  }else{
    K_hat <- sum((TRUNC < eigen_values))
    message(paste("Number of factors estimated given by 'K_hat' =", K_hat))
  }

  # APPEND RESULTS FOR EACH COEF INDEX TO THE RESULT LIST
  coef_counter <- 1
  for (COEF_INDEX in COEF_INDEX_VEC) {

  # Interception for COEF_INDEX
  if(COEF_INDEX > p | COEF_INDEX < 1 | floor(COEF_INDEX) != COEF_INDEX ){
    stop("COEF_INDEX must be an integer less than p. ")
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
  # Step 1 c): Computation of projection matrix
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
  # Step 2 b): Run the MAIN Lasso regression on the transformed data and get
  #            variance estimate
  #============================================================================#

    # Check for user-specified number of folds
  if(is.null(foldid)){

    if(obs_N %% NFOLDS > 0 & obs_N > NFOLDS){

      # If not divisible fill first folds with one leftover
      fold_vec <- c(rep(1:(obs_N %% NFOLDS), each = ((floor(obs_N/NFOLDS)+1) * obs_T)),
                  rep((obs_N %% NFOLDS + 1):NFOLDS, each = floor(obs_N/NFOLDS)*obs_T))

      message("Fold assignment based on NFOLDS modulo")

    }else{
      if(obs_N > NFOLDS){# If obs_N < NFOLDS use leave one out cv
        fold_vec <- rep(1:NFOLDS, each =   (obs_N/NFOLDS * obs_T))

        message("Fold assignment based on NFOLDS")

      }else{

        fold_vec <- rep(1:obs_N, each = obs_T)

        message("Leave one out cv used.")
      }
    }

    fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde,
                                        foldid = fold_vec, family = "gaussian",
                                        alpha = 1, intercept = FALSE)

    coefs_Lasso <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
    lambda_cv  <- fit_Lasso$lambda.min

  # Check if there is user-specified set of folds
   }else{
    fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde,
                                      foldid = foldid, family = "gaussian",
                                      alpha = 1, intercept = FALSE)

    coefs_Lasso <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
    lambda_cv  <- fit_Lasso$lambda.min

    message("User specified folds for CV selected.")
  }



  # Variance estimate using cross-validated penalty
  yhat_Lasso <- stats::predict(fit_Lasso, newx = X_tilde,
                        type = "response", s = "lambda.min")
  resid_Lasso <- Y_tilde - yhat_Lasso


  #============================================================================#
  # Step 2 b): Run the nodewise Lasso estimation on transformed data
  #============================================================================#

  # Check for user-specified number of folds
  if(is.null(foldid)){

    if(obs_N %% NFOLDS > 0 & obs_N > NFOLDS){

      # If not divisible fill first folds with one leftover
      foldid <- c(rep(1:(obs_N %% NFOLDS), each = ((floor(obs_N/NFOLDS)+1) * obs_T)),
                  rep((obs_N %% NFOLDS + 1):NFOLDS, each = floor(obs_N/NFOLDS)*obs_T))

      message("Fold assignment based on NFOLDS modulo.")

    }else{
      if(obs_N > NFOLDS){# If obs_N < NFOLDS use leave one out cv

        fold_vec <- rep(1:NFOLDS, each =   (obs_N/NFOLDS * obs_T))

        message("Fold assignment based on NFOLDS")

      }else{

        fold_vec <- rep(1:obs_N, each = obs_T)

        message("Leave one out cv used.")
      }
    }
    fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX],
                                   foldid = fold_vec, family = "gaussian",
                                   alpha = 1, intercept = FALSE)

    coefs_node_Lasso <- stats::coef(fit_node_Lasso, s = "lambda.min")[-1]
    kappa_cv  <- fit_node_Lasso$lambda.min



  }else{# Execute if there is user-specified set of folds
    fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX],
                                   foldid = foldid, family = "gaussian",
                                   alpha = 1, intercept = FALSE)

    coefs_node_Lasso <- stats::coef(fit_node_Lasso, s = "lambda.min")[-1]
    kappa_cv  <- fit_node_Lasso$lambda.min

    message("User specified folds for CV selected.")

  }



  # Kappa grid that cv.glmnet uses for cross-validation
  kappa_grid <- fit_node_Lasso$lambda

  # Index for the CV chosen kappa
  kappa_cv_idx <- fit_node_Lasso$index[1]

  #============================================================================#
  # Step 3: Calculate the de-sparsified estimator
  #============================================================================#

  #============================================================================#
  # Step 3 a): Choose the nodewise penalty parameter
  #============================================================================#

  # Kappa grid length
  kappa_grid_len <- length(kappa_grid)

  # Calculate the "scaled asymptotic variance"
  var_scaled <- rep(0, len = kappa_grid_len)

  # Compute Bartlett-kernel weight matrix W
  h_B <- 10 # Default Kernel Bandwidth
  W <- matrix(0, nrow = obs_T, ncol = obs_T)
  for (t_1 in 1:obs_T) {
    for (t_2 in 1:obs_T) {
      if(abs((t_1-t_2)/h_B) <= 1){
        W[t_1,t_2] <- 1 - abs((t_1-t_2)/h_B)
      }
    }
  }

 # This is the asymptotic variance estimator under heteroscedasity
    for(k in 1:kappa_grid_len){
      yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                                        type = "response", s = kappa_grid[k])
      resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso

      tmp <- numeric(obs_N)
      for (i in 1:obs_N) {
        resid_prod_components <- resid_node_Lasso[((i-1)*obs_T+1):(i*obs_T)]*resid_Lasso[((i-1)*obs_T+1):(i*obs_T)]
        tmp[i] <- t(resid_prod_components) %*% W %*% resid_prod_components
      }
      var_scaled[k] <-  (1/(t(X_tilde[,COEF_INDEX]) %*% resid_node_Lasso)^2)*sum(tmp)
    }


  # 25% increase of the scaled variance for truncation
  V_TRUNC <- 1.25 * var_scaled[kappa_cv_idx]

  # Indicator whether the scaled variance is less than the threshold
  kappa_idx <- 1
  for (l in 1:kappa_grid_len) {
    if(var_scaled[l] <= V_TRUNC){
      kappa_idx <- l
    }
    if(var_scaled[l] > V_TRUNC){
      break
    }
  }



  #============================================================================#
  # Step 3 b): Construction of the de-sparsified estimator
  #============================================================================#
  # Nodewise lasso residuals
  yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                             type = "response", s = kappa_grid[kappa_idx])
  resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso

  despar_beta <- coefs_Lasso[COEF_INDEX] + t(resid_node_Lasso) %*% resid_Lasso / t(resid_node_Lasso) %*% X_tilde[, COEF_INDEX]

  # Collect results
  Avar <- sqrt(var_scaled[kappa_idx])

  conf_band_min <- rep(despar_beta, length(alpha)) + Avar * qnorm(alpha/2)
  conf_band_max <- rep(despar_beta, length(alpha)) + Avar * qnorm(1-alpha/2)

  confidence_band <- cbind(conf_band_min, conf_band_max)

  if(length(COEF_INDEX_VEC) > 1){
  Results[[coef_counter]] <- list(coef_despar = despar_beta, Avar = Avar,
                               confidence_band = confidence_band)
  }
  else{
    Results <- list(coef_despar = despar_beta, Avar = Avar,
                    confidence_band = confidence_band)
  }

  coef_counter <- coef_counter+1
  }



  return(Results)
}





