#===============================================================================
# Internal helpers. Not exported.
#===============================================================================

# Kernels with compact support [-1, 1] -----------------------------------------
.hdcce_kernel <- function(kernel){
  if(is.function(kernel)) return(kernel)
  switch(kernel,
    epanechnikov = function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
    biweight     = function(u) ifelse(abs(u) <= 1, (15/16) * (1 - u^2)^2, 0),
    stop("kernel must be 'epanechnikov', 'biweight', or a function.")
  )
}

# Subtract the cross-sectional (over i) mean at each time point -----------------
# M is n_obs x q, sorted unit-major: rows (i-1)*obs_T + 1, ..., i*obs_T.
.hdcce_center_units <- function(M, obs_N, obs_T){
  M   <- as.matrix(M)
  out <- M
  for(t in 1:obs_T){
    idx <- seq(t, obs_N * obs_T, by = obs_T)
    out[idx, ] <- M[idx, , drop = FALSE] - matrix(colMeans(M[idx, , drop = FALSE]),
                                                  nrow = obs_N, ncol = ncol(M), byrow = TRUE)
  }
  out
}

# Apply a T x T matrix to every unit block --------------------------------------
.hdcce_project_units <- function(M, Pi, obs_N, obs_T){
  M   <- as.matrix(M)
  out <- M
  for(q in seq_len(ncol(M))){
    out[, q] <- as.vector(Pi %*% matrix(M[, q], nrow = obs_T, ncol = obs_N))
  }
  out
}
