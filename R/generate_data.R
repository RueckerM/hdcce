#' Generate toy panel data
#'
#' \code{generate_data} generates data as in simulation section of Linton, O.,
#'  Ruecker, M., Vogt, M.,  Walsh, C. (2024) "Estimation and Inference in
#'  High-Dimensional Panel Data Models with Interactive Fixed Effects".
#'  The example provides the code used to create the data sets
#'  \bold{data_estimation.rda} and \bold{data_inference.rda} of the package.
#' @param obs_N Number of cross-section units.
#' @param obs_T Number of time periods.
#' @param p Number of regressors/covariates.
#' @param mu Vector for individual loadings.
#' @param RHO Pairwise correlation coefficient of the regressors.
#' @return List containing the balanced panel data for estimation
#'  (data_estimation) and inference (data_inference). data_estimation$y contains
#'  the dependent variables and data_estimation$x the regressors. data_inference
#'  is in the same spirit with the minor difference that data_inference$y has three
#'  columns (one for each signal c** = 0, 0.1, 0.2). In both cases, the panels are sorted
#'  such  that first T observations are those for unit 1, followed by the T
#'  observations for unit 2, etc.
#'@examples
#'# Simulate an example data set as in Scenario B of the paper with
#'# N = 50, T = 50, RHO = 0.25, p = 901
#'#-----------------------------------------------------------------------------
#'
#'# Set the parameters
#'#-----------------------------------------------------------------------------
#'\dontrun{
#'# Set the number of cross-sections, time periods
#'obs_N  <- 50
#'obs_T  <- 50
#'p <- 901
#'# Specify the pairwise correlation among the regressors
#'RHO <- 0.25
#'# Specify the mean for the factor loadings
#'mu <- c(1, 1, 1, rep(1, (p-1)))
#'
#'# Simulate the data
#'set.seed(2024)
#'data <- generate_data(obs_N = obs_N, obs_T = obs_T, p,  mu = mu,
#'                              RHO = RHO)
#'data_estimation <- data$data_estimation
#'data_inference <- data$data_inference
#'save(data_estimation, file = "data_estimation.rda", compress = TRUE)
#'save(data_inference, file = "data_inference.rda", compress = TRUE)
#'}

generate_data <- function(obs_N, obs_T, p, mu, RHO){

  # K = 3 Factors
  K <- 3

  # Generate beta
  if((p-1) %% 3 > 0){
    stop("beta can not be generated with group structure.")
  }
  if(p < 10){
    beta <- rep(1, p)
    group_size <- (p-1)/3
    message( beta)
  }else{
   group_size <- (p-1)/3
   beta_group <- c(1,1,1, rep(0, (group_size-3)))
   beta <- c(1, rep(beta_group, 3))
   message( beta)
  }

  X <- numeric(p)
  Y <- numeric(1)

  Y_INFERENCE <- numeric(3) # Generate response variable for inference


  # Generate factors
  F_1 <- as.vector(stats::arima.sim(model = list(ar = 0.5), n = obs_T,
                             innov = rnorm(n = obs_T, mean = 0, sd = sqrt(0.75))))
  F_2 <- as.vector(stats::arima.sim(model = list(ar = 0.5), n = obs_T,
                             innov = rnorm(n = obs_T, mean = 0, sd = sqrt(0.75))))
  F_3 <- as.vector(stats::arima.sim(model = list(ar = 0.5), n = obs_T,
                             innov = rnorm(n = obs_T, mean = 0, sd = sqrt(0.75))))
  F_matrix <- cbind(F_1, F_2, F_3)

  # Oracle projection matrix
  Pi <- diag(1, obs_T, obs_T) - F_matrix %*% solve((t(F_matrix) %*% F_matrix)) %*% t(F_matrix)

  # Covariance matrix of G_i
  Sigma_G <- matrix(RHO, nrow = (K + p - 1 ),
                     ncol = (K + p - 1)) + diag(1 - RHO, nrow = (K + p - 1 ),
                                                      ncol = (K + p - 1))

  # Idiosyncratic error iid across i and t
  eps <- rnorm(n = obs_N * obs_T, mean = 0, sd = 1)

  # Nodewise errors iid across i and t
  u <- rnorm(n = obs_N * obs_T, mean = 0, sd = 1)

  # Idiosyncratic regressors
  Z <- mvtnorm::rmvnorm(n = obs_N * obs_T, mean = rep(0, (p - 1),
                                                        sd = diag(1, nrow = (p - 1),
                                                                  ncol = (p - 1) )))

  G_i <- mvtnorm::rmvnorm(obs_N, mean = mu, sigma = Sigma_G)



  for (i in 1:obs_N) {################# Start iteration i #################

    gamma_i <- G_i[i,1:K]

    Gamma_i <- cbind(c(G_i[i, (K + 1):(K + group_size)], rep(0, 2 * group_size)),
                      c(rep(0, group_size), G_i[i, (K + group_size + 1):(K + 2 * group_size)], rep(0, group_size)),
                      c(rep(0, 2 * group_size), G_i[i, (K + 2 * group_size + 1):(K + 3 * group_size)]))


    # Regressors generated with factor structure
    X_i <- F_matrix %*% t(Gamma_i) + Z[((i-1) * obs_T + 1) : (i * obs_T),]

    # First regressor generated with nodewise structure
    X_i1 <- X_i[, 1] * sqrt(2/3) + u[((i-1) * obs_T + 1) : (i * obs_T)]


    # Collect the design X_i in a first step
    X_i <- cbind(X_i1, X_i)


    # Generate response variable Y
    Y_i <- X_i %*% beta + F_matrix %*% gamma_i + eps[((i-1) * obs_T + 1) : (i * obs_T)]


    # Generate response Y for INFERENCE
    Y_INFERENCE_i <- matrix(nrow = obs_T, ncol = 3)


    k <- 1
    for (signal in c(0, 0.1, 0.2)) {
      beta_INFERENCE <- beta
      beta_INFERENCE[1] <- signal
      Y_INFERENCE_i[,k] <- X_i %*% beta_INFERENCE + F_matrix %*% gamma_i + eps[((i-1) * obs_T + 1) : (i * obs_T)]

      k <- k + 1
    }


    # Collect the data
    X <- rbind(X, X_i)
    Y <- rbind(Y, Y_i)

    Y_INFERENCE <- rbind(Y_INFERENCE, Y_INFERENCE_i)

  } ################# End iteration i #################

  return(list(data_estimation = list(x = X[-1,], y = Y[-1]),
         data_inference = list(x = X[-1,], y = Y_INFERENCE[-1,])))
}












