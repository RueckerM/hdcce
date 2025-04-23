#' Estimate HD panels with IFE
#'
#' \code{hdcce_estimator} fits the high-dimensional CCE estimation procedure
#' proposed in Linton, 0., Ruecker, M., Vogt, M., Walsh, C. (2024) "Estimation
#' and Inference in High-Dimensional Panel Data Models with Interactive
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
#' @param NFACTORS Allows to set the number of factors used in estimation.
#'   Default is NULL so that the data driven choice with TRUNC = 0.01 is used.
#' @param variant Choose whether to use the lasso based estimator
#'   (variant = "Lasso" (Default)) or the least squares variant (variant = "LS").
#' @param lambda User specified lambda grid.
#' @param NFOLDS The number of folds (partitioned along cross-section n) used for
#'  cross-validation (as described in the paper). Default is NFOLDS = 10.
#'  Fold size can vary by one if obs_N is not divisible by NFOLDS.
#'@param foldid Integer vector (obs_N*obs_T - dimensional) containing a label
#'   for each sample to determine its fold for CV.
#' @param scree_plot Logical variable to indicate whether a scree plot of the
#'   eigendecomposition of \eqn{\Sigma} should be shown. The default is TRUE.
#' @param standardize Logical variable to indicate whether glmnet is called
#'   with the standardized projected data. The default is TRUE.
#' @return The function returns the estimated coefficients, the
#' number of estimated factors and the lasso penalty parameter. More
#' specifically:
#' \itemize{
#' \item If a lambda sequence was supplied or variant = "LS", then the function
#' returns a list with
#' \describe{
#'   \item{$coefs}{The coefficient estimates.}
#'   \item{$K_hat}{The estimated number of factors. If NFACTORS was
#'   set by user then this will be returned.}
#'   }
#' \item If no lambda sequence was supplied and variant = "Lasso", then the function
#' returns a list with
#'  \describe{
#'   \item{$coefs}{The estimated coefficients.}
#'   \item{$K_hat}{The estimated number of factors. If NFACTORS was
#'   set by user then this will be returned.}
#'   \item{$Lambda}{The lambda_min calculated by cross-validation.}
#'   }
#' }
#' @examples
#' # Load the data set
#' data("data_estimation")
#'
#' # Set the dimensions of the data
#' obs_N <- 50
#' obs_T <- 50
#'
#' # Do the estimation
#' estimate_model <- hdcce_estimator(data = data_estimation, obs_N = obs_N,
#'                       obs_T = obs_T, TRUNC = 0.01, NFACTORS = NULL,
#'                       variant = "Lasso", lambda = NULL, NFOLDS = 10,
#'                       foldid = NULL, scree_plot = TRUE)
#' print(estimate_model$coefs[c(1:10)])
#'
#' @references  Linton, 0., Ruecker, M., Vogt, M., Walsh, C. (2024) "Estimation
#' and Inference in High-Dimensional Panel Data Models with Interactive
#' Fixed Effects"



#' @export
hdcce_estimator <- function(data, obs_N, obs_T, TRUNC = 0.01,
                               NFACTORS = NULL, variant = "Lasso",
                               lambda = NULL, NFOLDS = 10,
                               foldid = NULL, scree_plot = TRUE,
                               standardize = TRUE){
# Initial Checks
#--------------------------------------------------------------------------

  # Check to see whether a correct variant has been supplied
  if (!((variant == "LS") | (variant == "Lasso"))){
    stop('Variant must be set to "LS" or "Lasso."')
  }
  # Interception for foldid
  if(class(foldid) == "numeric"){
    if(!all(foldid == floor(foldid))){
      stop('Provided vector for CV must contain integers only.')
    }
    if(length(foldid) != (obs_N * obs_T)){
      stop('Provided vector for CV has wrong dimension.')
    }
  }
  # Interception for trunc
    if(TRUNC > 1 | TRUNC <= 0){
      stop('Supplied truncation invalid. Must be in (0,1].')
   }

  # Interception for NFOLDS
    if(NFOLDS >= obs_N){
      stop("Supplied number NFOLDS must be less then obs_N")
    }
    if(NFOLDS != floor(NFOLDS)){
      stop("Supplied number NFOLDS must be integer valued")
    }


  # Pull out the data
  X_data <- data$x
  Y_data <- data$y
  p <- dim(X_data)[2]

  # Interception for the data
  if(length(X_data[,1]) != obs_N *obs_T){
    stop("Supplied dimensions differ.")
  }
  if(length(Y_data) != obs_N * obs_T){
    stop("Supplied dimensions differ.")
  }

  # Interception for NFACTORS
  if(!is.null(NFACTORS)){
    if(NFACTORS >= p){
      stop("Supplied number NFACTORS must be less then p")
    }
    if(NFACTORS != floor(NFACTORS)){
      stop("Supplied number NFACTORS must be integer valued")
    }
  }



  #============================================================================#
  # Step 1: Eigendecomposition of empirical covariance matrix
  #============================================================================#

  # Cross-sectional averages of the regressors
  X_bar <- matrix(NA, ncol = p, nrow = obs_T)
  for(t in 1:obs_T){
     indices <- seq(t,obs_N * obs_T, by = obs_T)
     X_bar[t,] <- colMeans(X_data[indices,])
  }

  # Empirical covariance matrix and eigenstructure
  Cov_X_bar <- 1/obs_T * t(X_bar) %*% X_bar
  Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric=TRUE)



  #============================================================================#
  # Step 2: Estimation of number of factors
  #============================================================================#

  # Normalize the eigenvalues with the largest one
  eigen_values <- Cov_X_bar_eigen$values / (Cov_X_bar_eigen$values[1])

  # Check for user-specified fixed number of factors
  if(class(NFACTORS) == "numeric"){
    if((floor(NFACTORS) == NFACTORS)){
      K_hat <- NFACTORS
      message(paste("User-supplied number of factors given by 'NFACTORS' = ",
                    NFACTORS," is used in estimation."))

      if( scree_plot == TRUE){
        graphics::plot(eigen_values, ylim=c(0,1), ylab = "Normalized Eigenvalues",
                       main=paste("Number of factors set to ", K_hat))
        graphics::points(eigen_values[1:K_hat], col="red")
      }
    } else{
      stop("Supplied numer of factors NFACTORS is not an integer.")
    }
  } else {
    # Number of normalized eigenvalues larger than TRUNC
    K_hat <- sum((TRUNC < eigen_values))
    message(paste("Number of factors estimated given by 'K_hat' =", K_hat))

    if(scree_plot == TRUE){
      graphics::plot(eigen_values, ylim = c(0,1), ylab = "Normalized Eigenvalues",
                     main = paste("Estimated number of factors = ", K_hat))
      graphics::abline(h = TRUNC, col = "red")
      graphics::legend("topright", legend = c("Truncation"), lty = c(1),
                       col = c("red"))
    }
  }



  #============================================================================#
  # Step 3: Computation of projection matrix
  #============================================================================#

  W_hat <- X_bar %*% Cov_X_bar_eigen$vectors[,1:K_hat]
  Pi_hat <- diag(obs_T) -  W_hat %*% solve(t(W_hat) %*% W_hat)  %*% t(W_hat)

  #============================================================================#
  # Step 4: Transform the data
  #============================================================================#

  Y_hat <- rep(NA, obs_T * obs_N)
  X_hat <- matrix(NA, nrow = obs_N * obs_T, ncol = p)

  for(i in 1:obs_N){
    indices <- ((i-1) * obs_T + 1):(i * obs_T)
    Y_hat[indices] <- Pi_hat %*% t(t(Y_data[indices]))
    X_hat[indices,] <- Pi_hat %*% t(t(X_data[indices,]))
  }

  #============================================================================#
  # Step 5: Estimate Lasso or OLS on the transformed data
  #============================================================================#

  # Run the LS variant if it was supplied
  if (variant == "LS"){
    # Get least squares estimates
    coef_est <- lm(Y_hat ~ X_hat - 1)$coefficients
    results <- list(coefs = coef_est, K_hat = K_hat)
    message("LS variant selected.")
  }# Close the variant == "LS" block

  # Run the Lasso variant if it was supplied
  if(variant == "Lasso"){
    # Execute if user specified lambda grid is provided
    if(is.null(lambda) == FALSE){
          fit_Lasso <- glmnet::glmnet(x = X_hat, y = Y_hat, family = "gaussian",
                                      alpha = 1, lambda = lambda/2,
                                      intercept = FALSE,
                                      standardize = standardize)

          message("User specified lambda grid selected.")
          results <- list(coefs = fit_Lasso$beta, K_hat = K_hat)

    }else {# Run the cross-validated code as a default
      # Execute if no user specific foldid
      if(is.null(foldid)){
          message(paste("User-supplied number of folds given by 'NFOLDS' = ",
                        NFOLDS," is used to create fold vector."))
          if(obs_N %% NFOLDS > 0 & obs_N > NFOLDS){ # If not divisible fill first folds with one leftover
            foldid <- c(rep(1:(obs_N %% NFOLDS),
                            each = ((floor(obs_N/NFOLDS)+1) * obs_T)),
                        rep((obs_N %% NFOLDS + 1):NFOLDS,
                            each = floor(obs_N/NFOLDS)*obs_T))

            message("Fold assignment based on NFOLDS with modulo")

          }else{
            if(obs_N > NFOLDS){# If obs_N < NFOLDS use leave one out cv
              foldid <- rep(1:NFOLDS, each =   (obs_N/NFOLDS * obs_T))

              message("Fold assignment based on NFOLDS")
            }
            else{
              foldid <- rep(1:obs_N, each = obs_T)

              message("Leave one out cv used.")
            }
          }

          fit_Lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat,
                                            foldid = foldid, family = "gaussian",
                                            alpha = 1, intercept = FALSE,
                                            standardize = standardize)

          coefs <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
          Lambda_CV  <- fit_Lasso$lambda.min

          results <- list(coefs = coefs, K_hat = K_hat, Lambda = Lambda_CV)


      # Execute if there is user-specified set of folds
      }else{
          fit_Lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat,
                                            foldid = foldid, family = "gaussian",
                                            alpha = 1, intercept = FALSE,
                                            standardize = standardize)
          coefs <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
          Lambda_CV  <- fit_Lasso$lambda.min

          message("User specified folds for CV selected.")

          results <- list(coefs = coefs, K_hat = K_hat, Lambda = Lambda_CV)

     }
  }

}# Close the variant == "Lasso" block


 # Return the estimation results
  return(results)

}# Close the estimation function





