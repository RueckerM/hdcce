#' Simulated panel data for estimation
#'
#' A balanced panel simulated with \code{\link{generate_data}} from the design
#' of the simulation study in Rücker et al. (2025), used to
#' illustrate \code{\link{hdcce_estimator}}. The panel is sorted such that the
#' first \code{obs_T = 20} rows belong to unit 1, the next 20 to unit 2, etc.
#'
#' @format A list with two components:
#' \describe{
#'   \item{x}{numeric matrix of dimension 400 x 61 containing the regressors.}
#'   \item{y}{numeric vector of length 400 containing the response.}
#' }
#' @source Simulated with \code{generate_data(obs_N = 20, obs_T = 20, p = 61,
#'   mu = rep(1, 63), RHO = 0.25)} after \code{set.seed(2024)}.
#' @references  Rücker, M., Vogt, M., Linton, O. and Walsh, C. (2025).
#'   Estimation and inference in high-dimensional panel data models with
#'   interactive fixed effects. \emph{Quantitative Economics}, 16(4),
#'   1457--1509. \doi{10.3982/QE2308}.
#'
#'   Rücker, M., Vogt, M. and Linton, O. (2026) "High-Dimensional Panel Data Models with
#'   Interactive Fixed Effects: Beyond the Linear Case" \url{https://arxiv.org/pdf/2608.02055}.
"data_estimation"


#' Simulated panel data for inference
#'
#' As \code{\link{data_estimation}}, but with three response vectors differing
#' in the coefficient of the first regressor; used to illustrate
#' \code{\link{hdcce_inference}}.
#'
#' @format A list with two components:
#' \describe{
#'   \item{x}{numeric matrix of dimension 400 x 61 containing the regressors.}
#'   \item{y}{numeric matrix of dimension 400 x 3; column \eqn{k} is generated
#'     with coefficient of the first regressor equal to c** = 0, 0.1 and 0.2.}
#' }
#' @source Simulated with \code{generate_data(obs_N = 20, obs_T = 20, p = 61,
#'   mu = rep(1, 63), RHO = 0.25)} after \code{set.seed(2024)}.
#' @references Rücker, M., Vogt, M., Linton, O. and Walsh, C. (2025).
#'   Estimation and inference in high-dimensional panel data models with
#'   interactive fixed effects. \emph{Quantitative Economics}, 16(4),
#'   1457--1509. \doi{10.3982/QE2308}
#'
#'   Rücker, M., Vogt, M. and Linton, O. (2026) "High-Dimensional Panel Data Models with
#'   Interactive Fixed Effects: Beyond the Linear Case" <https://arxiv.org/pdf/2608.02055>.
"data_inference"
