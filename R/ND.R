#' Network Deconvolution
#'
#' Clean up a network adjacency matrix, filter the false-positive edges.
#'
#' Feizi, S.;  Marbach, D.;  Médard, M.; Kellis, M., Network deconvolution as a general method to distinguish direct dependencies in networks. Nature Biotechnology 2013, 31, 726–733.
#'
#' @param mat Input matrix, if it is a square matrix, the program assumes Input matrix,
#' if it is a square matrix, the program assumes between nodes i and j. Elements of matrix should be
#' non-negative.
#' @param beta Scaling parameter, the program maps the largest absolute eigenvalue
#' Scaling parameter, the program maps the largest absolute eigenvalue between 0 and 1.
#' @param alpha fraction of edges of the observed dependency matrix to be kept in deconvolution process.
#' @param control If FALSE, displaying direct weights for observed interactions, if 1, displaying direct
#' weights for both observed and non-observed interactions.
#' @param linear_mapping_before If TRUE, mat will be linearly mapped to be between 0 and 1 before deconvolution.
#' @param linear_mapping_after If TRUE, result will be linearly mapped to be between 0 and 1 after deconvolution.
#'
#' @return mat_nd, Output deconvolved matrix (direct dependency matrix). Its components
#' represent direct edge weights of observed interactions. Choosing top direct interactions (a cut-off) depends on the application and
#' is not implemented in this code.
#'
#' @export
#'
#' @examples
#' a <- matrix(1:9, nrow = 3)
#' ND(a)
#'
ND <- function(mat, beta = 0.99, alpha = 1, control = FALSE, linear_mapping_before = TRUE, linear_mapping_after = TRUE) {
  stopifnot(is.matrix(mat))
  stopifnot(nrow(mat) == ncol(mat))
  stopifnot(all(mat >= 0))
  stopifnot(beta > 0 & beta < 1)
  stopifnot(alpha > 0 & alpha <= 1)
  mat_max = max(mat)
  mat_min = min(mat)
  if (mat_max == mat_min) {
    stop("the input matrix is a constant matrix")
  }

  ################ preprocessing the input matrix ################

  # linearly mapping the input matrix to be between 0 and 1
  if (linear_mapping_before) {
    mat = (mat - mat_min) / (mat_max - mat_min)
  }
  # diagonal values are filtered, as 0
  diag(mat) <- 0.
  # filtered the edges
  y <- quantile(mat, 1 - alpha)
  mat_th <- mat
  mat_th[mat_th < y] <- 0.
  # making the matrix symetric if already not
  mat_th <- (mat_th + t(mat_th)) / 2

  ################ eigen decomposition ################
  eigen_fit <- eigen(mat_th)
  U <- eigen_fit$vectors
  D <- eigen_fit$values

  lam_n <- min(c(min(D), 0))
  lam_p <- max(c(max(D), 0))
  m1 <- lam_p * (1 - beta) / beta
  m2 <- lam_n * (1 + beta) / beta
  m <- max(m1, m2)

  ################ network deconvolution ################
  D <- D / (D + m)
  mat_new1 <- U %*% diag(D) %*% solve(U)

  ################ displying direct weights ################
  if (control) {
    m2 <- min(mat_new1)
    mat_new2 <- mat_new1 + max(-m2, 0)
  } else {
    ind_nonedges <- mat_th == 0.
    m1 <- max(mat * ind_nonedges)
    m2 <- min(mat_new1)
    mat_new2 <- mat_new1 + max(m1 - m2, 0)
    mat_new2[ind_nonedges] <- mat[ind_nonedges]
  }

  if (linear_mapping_after) {
    m1 <- min(mat_new2)
    m2 <- max(mat_new2)
    mat_nd <- (mat_new2 - m1) / (m2 - m1)
  } else {
    mat_nd <- mat_new2
  }

  mat_nd
}
