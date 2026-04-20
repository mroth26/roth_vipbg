################################################################################
# simulate_twins: Generate twin data for a bivariate ACE + DoC + sibling-interaction model
#
# CONVENTIONS:
#   - Column order (phenotype-major): p1_t1, p1_t2, p2_t1, p2_t2
#   - Zygosity codes: zyg {MZ=1, DZ=2}
#   - DoC transform: Tm = (I_2 - B)^(-1)
#   - Sibling interaction transform: G_sib = (I_4 - Bsib)^(-1)
################################################################################

simulate_twins <- function(
    n_MZ, n_DZ,
    seed = NULL,
    rA = 0, rC = 0, rE = 0,
    rA12 = NULL, rC12 = NULL, rE12 = NULL,
    b12 = 0, b21 = 0,
    s = c(0, 0),
    A_diag = c(0.5, 0.5),
    C_diag = c(0.3, 0.3),
    E_diag = c(0.2, 0.2),
    mean = c(0, 0)
) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required.")
  }
  if (!is.null(seed)) set.seed(seed)

  s <- as.numeric(s)
  if (length(s) == 1) s <- rep(s, 2)
  if (length(s) != 2) stop("s must be a scalar or length-2 numeric vector.")

  A_diag <- as.numeric(A_diag)
  C_diag <- as.numeric(C_diag)
  E_diag <- as.numeric(E_diag)
  mean   <- as.numeric(mean)

  if (length(A_diag) != 2) stop("A_diag must be length 2.")
  if (length(C_diag) != 2) stop("C_diag must be length 2.")
  if (length(E_diag) != 2) stop("E_diag must be length 2.")
  if (length(mean)   != 2) stop("mean must be length 2.")

  tot <- A_diag + C_diag + E_diag
  if (any(!is.finite(tot)) || any(tot <= 0)) {
    stop("A_diag + C_diag + E_diag must be positive for each trait.")
  }

  # trait-wise normalization
  A_diag <- A_diag / tot
  C_diag <- C_diag / tot
  E_diag <- E_diag / tot

  make_R2 <- function(r_default, r12) {
    r_eff <- if (!is.null(r12)) as.numeric(r12) else as.numeric(r_default)
    R <- matrix(c(1, r_eff, r_eff, 1), 2, 2, byrow = TRUE)
    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) {
      stop("Correlation matrix is not positive definite. Check rA/rC/rE inputs.")
    }
    R
  }

  R_A <- make_R2(rA, rA12)
  R_C <- make_R2(rC, rC12)
  R_E <- make_R2(rE, rE12)

  make_cov <- function(diagv, R) {
    Dm <- diag(sqrt(diagv))
    Dm %*% R %*% Dm
  }

  A <- make_cov(A_diag, R_A)
  C <- make_cov(C_diag, R_C)
  E <- make_cov(E_diag, R_E)

  # DoC transform
  B <- matrix(c(
    0,   b12,
    b21, 0
  ), 2, 2, byrow = TRUE)

  I2 <- diag(2)
  IB <- I2 - B
  if (abs(det(IB)) <= .Machine$double.eps) {
    stop("I_2 - B is not invertible.")
  }
  Tm <- solve(IB)

  # within-twin and cross-twin blocks in trait space
  V    <- Tm %*% (A + C + E) %*% t(Tm)
  X_MZ <- Tm %*% (A + C)     %*% t(Tm)
  X_DZ <- Tm %*% (0.5 * A + C) %*% t(Tm)

  # build 4x4 family covariance directly in phenotype-major order:
  # [p1_t1, p1_t2, p2_t1, p2_t2]
  Sigma_interleaved <- function(V, X) {
    matrix(c(
      V[1,1], X[1,1], V[1,2], X[1,2],
      X[1,1], V[1,1], X[1,2], V[1,2],
      V[2,1], X[2,1], V[2,2], X[2,2],
      X[2,1], V[2,1], X[2,2], V[2,2]
    ), 4, 4, byrow = TRUE)
  }

  Sigma_MZ <- Sigma_interleaved(V, X_MZ)
  Sigma_DZ <- Sigma_interleaved(V, X_DZ)

  # means: transform through DoC, then expand to phenotype-major family mean
  mu_obs <- as.vector(Tm %*% mean)
  mu4 <- c(mu_obs[1], mu_obs[1], mu_obs[2], mu_obs[2])

  # sibling interaction in phenotype-major order
  Bsib <- matrix(c(
    0,    s[1], 0,    0,
    s[1], 0,    0,    0,
    0,    0,    0,    s[2],
    0,    0,    s[2], 0
  ), 4, 4, byrow = TRUE)

  I4 <- diag(4)
  I4_Bsib <- I4 - Bsib
  if (abs(det(I4_Bsib)) <= .Machine$double.eps) {
    stop("I_4 - Bsib is not invertible.")
  }
  G_sib <- solve(I4_Bsib)

  Sigma_MZ <- G_sib %*% Sigma_MZ %*% t(G_sib)
  Sigma_DZ <- G_sib %*% Sigma_DZ %*% t(G_sib)
  mu4 <- as.vector(G_sib %*% mu4)

  check_sigma <- function(S, name) {
    if (!is.matrix(S) || any(dim(S) != c(4, 4))) {
      stop(paste0(name, " must be a 4x4 matrix."))
    }
    if (any(!is.finite(S))) {
      stop(paste0(name, " contains non-finite values."))
    }
    if (max(abs(S - t(S))) > 1e-10) {
      stop(paste0(name, " is not symmetric."))
    }
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-10) {
      stop(paste0(name, " is not positive definite (min eigenvalue = ",
                  signif(min(ev), 6), ")."))
    }
    TRUE
  }

  check_sigma(Sigma_MZ, "Sigma_MZ")
  check_sigma(Sigma_DZ, "Sigma_DZ")

  draw <- function(n, mu4, Sigma4) {
    if (n <= 0) return(NULL)
    mvtnorm::rmvnorm(n, mean = mu4, sigma = Sigma4)
  }

  to_df <- function(X, z) {
    if (is.null(X)) return(NULL)
    df <- data.frame(X)
    names(df) <- c("p1_t1", "p1_t2", "p2_t1", "p2_t2")
    df$zyg <- z
    df
  }

  X_MZ_draw <- draw(n_MZ, mu4, Sigma_MZ)
  X_DZ_draw <- draw(n_DZ, mu4, Sigma_DZ)

  dat <- do.call(rbind, list(
    to_df(X_MZ_draw, 1L),
    to_df(X_DZ_draw, 2L)
  ))
  dat <- dat[, c("p1_t1", "p1_t2", "p2_t1", "p2_t2", "zyg")]

  truth <- list(
    flags = list(),
    varparam = list(
      va11 = A[1,1],
      va22 = A[2,2],
      va12 = A[1,2],
      vc11 = C[1,1],
      vc22 = C[2,2],
      vc12 = C[1,2],
      ve11 = E[1,1],
      ve22 = E[2,2],
      ve12 = E[1,2],
      b12 = B[1,2],
      b21 = B[2,1]
    ),
    realized = list(
      rA12 = R_A[1,2],
      rC12 = R_C[1,2],
      rE12 = R_E[1,2],
      b12 = B[1,2],
      b21 = B[2,1],
      s1 = s[1],
      s2 = s[2]
    ),
    matrices = list(
      A = A,
      C = C,
      E = E,
      B = B,
      Tm = Tm,
      RA = R_A,
      RC = R_C,
      RE = R_E,
      V = V,
      X_MZ = X_MZ,
      X_DZ = X_DZ,
      Sigma_MZ = Sigma_MZ,
      Sigma_DZ = Sigma_DZ,
      Bsib = Bsib,
      G_sib = G_sib
    )
  )

  list(data = dat, truth = truth)
}
