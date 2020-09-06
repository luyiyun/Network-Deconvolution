#' Network deconvolution for gene regulatory networks of DREAM5
#'
#' This function can suit for non-square matrix
#'
#' Feizi, S.;  Marbach, D.;  Médard, M.; Kellis, M., Network deconvolution as a general method to distinguish direct dependencies in networks. Nature Biotechnology 2013, 31, 726–733.
#'
#'
#' @param mat Input matrix, it is a n_tf by n matrix where first n_tf genes are TFs.
#' Elements of the input matrix are nonnegative.
#' @param beta Scaling parameter, the program maps the largest absolute eigenvalue of the direct dependency matrix to beta. It should be
#' between 0 and 1. You should skip this scaling step if you know eigenvalues of your matrix satisfy ND conditions.
#' @param alpha fraction of edges of the observed dependency matrix to be kept in deconvolution process.
#' @param linear_mapping_before If TRUE, mat will be linearly mapped to be between 0 and 1 before deconvolution.
#' @param linear_mapping_after If TRUE, result will be linearly mapped to be between 0 and 1 after deconvolution.
#' @param control_p If set to TRUE, it perturbs input networks slightly to have stable results in case of non-diagonalizable matrices.
#' If FALSE, it checks some sufficient condition and then add a small perturbation (this may be slower). Default is FALSE.
#'
#' @return mat_nd, Output deconvolved matrix (direct dependency matrix). Its components
#' represent direct edge weights of observed interactions. Choosing top direct interactions (a cut-off) depends on the application and
#' is not implemented in this code.
#'
#' @export
#'
#' @examples
#' aa <- matrix(1:12, nrow=3)
#' ND_regulatory(aa)
ND_regulatory <- function(mat, beta = 0.5, alpha = 0.1, linear_mapping_before = TRUE, linear_mapping_after = TRUE, control_p = FALSE) {
  stopifnot(is.matrix(mat))
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

  # making the TF-TF network symmetric
  # note some algorithms only output one-directional edges
  # the other direction is added in that case
  rr <- nrow(mat)
  cc <- ncol(mat)
  if (rr > cc) {
    mat <- t(mat)
    temp <- rr
    cc <- rr
    rr <- temp
  }

  # 将tf-tf的网络变成无向网。
  #   如果只有一条有方向的边，则再复制一遍其反方向的边即可
  #   如果两个方向都有边，则平均一下。
  tf_net <- mat[, 1:rr]
  ind <- (tf_net == 0) & (t(tf_net) != 0)
  tf_net[ind] <- tf_net[t(ind)]
  tf_net <- (tf_net + t(tf_net)) / 2
  mat[, 1:rr] <- tf_net

  # filtered the edges
  y <- quantile(mat, 1 - alpha)
  mat_th <- mat
  mat_th[mat_th < y] <- 0.

  # 因为滤过了一些边，所以还需要再进行一次对称化
  mat_th[, 1:rr] <- (mat_th[, 1:rr] + t(mat_th[, 1:rr])) / 2
  temp_net <- mat_th > 0
  temp_net_remain <- mat_th == 0.
  mat_th_remain <- mat[temp_net_remain]
  m11 <- max(mat_th_remain)

  ################ network deconvolution ################

  # check if matrix is diagonalizable
  if (!control_p) {
    mat1 <- rbind(mat, matrix(0, nrow = cc - rr, ncol = cc))
    eigen_fit <- eigen(mat1)
    U <- eigen_fit$vectors
    D <- eigen_fit$values
    if ((rcond(U)) < 1e10) {
      control_p <- TRUE
    }
  }

  # if matrix is not diagonalizable,
  # add random perturbation to make it diagonalizable
  if (control_p) {
    r_p <- 0.001
    rand_tf <- r_p * matrix(runif(rr*rr), nrow = rr)
    rand_tf <- (rand_tf + t(rand_tf)) / 2
    diag(rand_tf) <- 0
    rand_target <- r_p * matrix(runif(rr*(cc-rr)), nrow = rr)
    mat_rand <- cbind(rand_tf, rand_target)
    mat_th <- mat_th + mat_rand

    mat1 <- rbind(mat, matrix(0, nrow = cc - rr, ncol = cc))
    eigen_fit <- eigen(mat1)
    U <- eigen_fit$vectors
    D <- eigen_fit$values
  }

  lam_n <- min(c(min(D), 0))
  lam_p <- max(c(max(D), 0))
  m1 <- lam_p * (1 - beta) / beta
  m2 <- lam_n * (1 + beta) / beta
  m <- max(m1, m2)

  ################ network deconvolution ################
  D <- D / (D + m)
  mat_new1 <- U %*% diag(D) %*% solve(U)

  ################ displying direct weights ################
  # adding remaining edges
  mat_new2 <- mat_new1[1:rr, ]
  m2 <- min(mat_new2)
  mat_new3 <- mat_new2 + max(m11 - m2, 0)
  mat_new3[temp_net_remain] <- mat_th_remain

  mat_new3
}
