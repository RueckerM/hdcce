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
#' @param obs_N,obs_T The number of cross-section units (obs_N) and (obs_T) in
#'   the data.
#' @param TRUNC The truncation parameter tau used to estimate the number of
#'   factors. Default is 0.01.
#' @param NFACTORS Allows to set the number of factors used in estimation.
#'   Default is NULL so that the data driven choice with TRUNC = 0.01 is used.
#' @param variant Choose whether to use the lasso based estimator
#'   (variant = "Lasso" (Default)) or the least squares variant (variant = "LS").
#' @param lambda User specified lambda grid.
#' @param NFOLDS The number of folds used for cross-validation
#'   as described in the paper. Default is NFOLDS = 10. Fold size can vary by one
#'   if obs_N is not divisible by NFOLDS.
#'@param foldid Integer vector (obs_N*obs_T - dimensional) containing a label
#'   for each sample to determine its fold for CV.
#' @param scree_plot Logical variable to indicate whether a scree plot of the
#'   eigendecomposition of \eqn{\Sigma} should be shown. The default is TRUE.
#' @return The function returns the estimated coefficients, the
#' number of estimated factors and the lasso penalty parameter. More
#' specifically:
#' \itemize{
#' \item Unless penalty_choice = "None", the function returns a list with
#' \describe{
#'   \item{$coefs}{The coefficient estimates.}
#'   \item{$K_hat}{The estimated number of factors. If NFACTORS was
#'   set by user then this will be returned.}
#'   \item{$Lambda}{The data driven lasso penalty choice. If variant ="LS",
#'   it is non-existent.}
#'   }
#' \item If the option lambda = "lambda" was set, then the function
#' returns a list with
#'  \describe{
#'   \item{$coefs}{The coefficient estimates.}
#'   \item{$K_hat}{The estimated number of factors. If NFACTORS was
#'   set by user then this will be returned.}
#'   \item{$Lambda}{The values of the lasso penalty parameters used. If
#'   variant = "LS", it is non-existent.}
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
#' estimate_model <- hdcce_estimator(data = data_example, obs_N = obs_N,
#'                       obs_T = obs_T, TRUNC = 0.01, NFACTORS = NULL,
#'                       variant = "Lasso", lambda = NULL, NFOLDS = NULL,
#'                       foldid = NULL, scree_plot = TRUE)
#'
#' @references  Linton, 0., Ruecker, M., Vogt, M., Walsh, C. (2024) "Estimation
#' and Inference in High-Dimensional Panel Data Models with Interactive
#' Fixed Effects"




hdcce_estimator <- function(data, obs_N, obs_T, TRUNC = 0.01,
                               NFACTORS = NULL, variant = "Lasso",
                               lambda = NULL, NFOLDS = NULL,
                               foldid = NULL, scree_plot = TRUE){
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
    if(TRUNC > 1 | TRUNC <=0){
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
  if(length(X_data[,1]) != obs_N *obs_T){
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
    if(is.naturalnumber(NFACTORS) == TRUE){
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

    if( scree_plot == TRUE){
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
                                      alpha = 1, lambda = lambda,
                                      intercept = FALSE)
          message("User specified lambda grid selected.")
          results <- list(coefs = fit_Lasso$beta, K_hat = K_hat,
                          Lambda = fit_Lasso$lambda)
    }else {# Run the cross-validated code as a default
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

          fit_Lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, nfolds = NFOLDS,
                                            foldid = foldid, family = "gaussian",
                                            alpha = 1, intercept = FALSE)

          coefs <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
          Lambda_CV  <- fit_Lasso.cv$lambda.min

          results <- list(coefs = coefs, K_hat = K_hat, Lambda = Lambda_CV)
      }

      # Check if there is user-specified set of folds
      } else if(is.null(NFOLDS) == FALSE & is.null(foldid) == FALSE){
          fit_Lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, nfolds = NFOLDS,
                                            foldid = foldid, family = "gaussian",
                                            alpha = 1, intercept = FALSE)
          coefs <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
          Lambda_CV  <- fit_Lasso.cv$lambda.min
          message("User specified folds for CV selected.")
          results <- list(coefs = coefs, K_hat = K_hat, Lambda = Lambda_CV)

     } else{ # Neither specified NFOLDS nor foldid
      if(obs_N %% 10 > 0){ # If not divisible fill first folds with one leftover
         foldid <- c(rep(1:(obs_N %% 10), each = ((floor(obs_N/10)+1) * obs_T)),
                       rep((obs_N %% 10 + 1):10, each = floor(obs_N/10)*obs_T))
      } else{
        foldid <- rep(1:10, each =  (obs_N/10 * obs_T))
      }

       fit_Lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, family = "gaussian",
                                         foldid = foldid, alpha = 1,
                                         intercept = FALSE)
       coefs <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
       Lambda_CV  <- fit_Lasso.cv$lambda.min
       message("10 Fold CV as explained in paper.")
       results <- list(coefs = coefs, K_hat = K_hat, Lambda = Lambda_CV)
     }
  }# Close the variant == "Lasso" block

 # Return the estimation results
  return(results)
}# Close the estimation function





