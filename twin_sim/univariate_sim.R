################################################################################
# simulate_twins: Generate twin data for a univariate model
################################################################################
library(mvtnorm)
simulate_twins <- function(
    n_MZ, n_DZ,
    seed = NULL,
    s = 0,
    A_diag = 0.5, C_diag = 0.3, D_diag = 0, E_diag = 0.2,
    mean = 0
) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Package 'mvtnorm' is required.")
  if (!is.null(seed)) set.seed(seed)
  
  # resolve inputs
  A_diag <- as.numeric(A_diag)
  C_diag <- as.numeric(C_diag)
  D_diag <- as.numeric(D_diag)
  E_diag <- as.numeric(E_diag)
  tot <- A_diag + C_diag + D_diag + E_diag
  
  A_diag <- A_diag / tot
  C_diag <- C_diag / tot
  D_diag <- D_diag / tot
  E_diag <- E_diag / tot
  
  mean <- as.numeric(mean)
  s <- as.numeric(s)
  
  # latent A, C, D, E
  A <- A_diag
  C <- C_diag
  D <- D_diag
  E <- E_diag
  
  # observed-scale (phenotypic) variance
  V  <- A + C + D + E
  cMZ <- A + D + C
  cDZ <- 0.5 * A + 0.25 * D + C
  
  # 2x2 covariance matrices for MZ and DZ
  Sigma_MZ <- matrix(c(
    V,    cMZ,
    cMZ,  V
  ), 2, 2, byrow = TRUE)
  
  Sigma_DZ <- matrix(c(
    V,    cDZ,
    cDZ,  V
  ), 2, 2, byrow = TRUE)
  
  # means
  mu2 <- c(mean, mean)
  
  # apply sibling interaction transformation 
  if (s != 0) {
    I2 <- diag(2)
    Bsib <- matrix(c(
      0, s,
      s, 0
    ), 2, 2, byrow = TRUE)
    
    I2_Bsib <- I2 - Bsib
    if (abs(det(I2_Bsib)) <= .Machine$double.eps) stop("I2 - Bsib is not invertible.")
    G_sib <- solve(I2_Bsib)
    
    # transform covariances and means
    Sigma_MZ <- G_sib %*% Sigma_MZ %*% t(G_sib)
    Sigma_DZ <- G_sib %*% Sigma_DZ %*% t(G_sib)
    mu2 <- as.vector(G_sib %*% mu2)
  } else {
    Bsib <- matrix(0, 2, 2)
    G_sib <- diag(2)
  }
  
  # draws
  draw <- function(n, mu2, Sigma2) {
    if (n <= 0) return(NULL)
    mvtnorm::rmvnorm(n, mean = mu2, sigma = Sigma2)
  }
  to_df <- function(X, z) {
    if (is.null(X)) return(NULL)
    df <- data.frame(X)
    names(df) <- c("p1_t1", "p1_t2")
    df$zyg <- z
    df
  }
  
  X_MZ <- draw(n_MZ, mu2, Sigma_MZ)
  X_DZ <- draw(n_DZ, mu2, Sigma_DZ)
  
  dat <- do.call(rbind, list(
    to_df(X_MZ, 1L),
    to_df(X_DZ, 2L)
  ))
  dat <- dat[, c("p1_t1", "p1_t2", "zyg")]
  
  # truth object with univariate parameters
  truth <- list(
    flags = list(),
    varparam = list(
      va11 = A,
      vc11 = C,
      vd11 = D,
      ve11 = E
    ),
    realized = list(
      s1 = s
    ),
    matrices = list(
      A = A, C = C, D = D, E = E,
      V = V,
      Bsib = Bsib,
      G_sib = G_sib
    )
  )
  
  list(data = dat, truth = truth)
}
