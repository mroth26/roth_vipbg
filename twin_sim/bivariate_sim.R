################################################################################
# simulate_twins: Generate twin data for a bivariate ACE + DoC + sibling
#                 interaction model
#
# CONVENTIONS:
#   - Two phenotypes, each with latent A, C, E components
#   - Column order (phenotype-major): p1_t1, p1_t2, p2_t1, p2_t2
#   - Zygosity codes: zyg {MZ = 1, DZ = 2}
#   - DoC transform:           T     = (I_2 - B)^(-1)
#   - Sibling interaction:     G_sib = (I_4 - B_sib)^(-1)
#     where B_sib acts on the 4-vector [p1_t1, p1_t2, p2_t1, p2_t2], coupling
#     the two twins within each trait via s1 (trait 1) and s2 (trait 2).
#
# MODEL SPECIFICATION:
#   A = D_A R_A D_A,   C = D_C R_C D_C,   E = D_E R_E D_E   (each 2x2)
#   B = [[0, b12], [b21, 0]]
#   V    = T (A + C + E) T'
#   X_MZ = T (A + C) T'
#   X_DZ = T (0.5 A + C) T'
#   Sigma_MZ, Sigma_DZ are 4x4 in phenotype-major order, then transformed by
#   the sibling interaction:
#     Sigma_*_sib = G_sib Sigma_* G_sib'
#     mu_4_sib    = G_sib mu_4
################################################################################

simulate_twins <- function(
    n_MZ, n_DZ,
    seed = NULL,
    rA = 0, rC = 0, rE = 0,
    rA12 = NULL, rC12 = NULL, rE12 = NULL,
    b12 = 0, b21 = 0,
    s = c(0, 0),
    A_diag = c(0.5, 0.5), C_diag = c(0.3, 0.3), E_diag = c(0.2, 0.2),
    mean = c(0, 0)
) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Package 'mvtnorm' is required.")
  if (!is.null(seed)) set.seed(seed)

  # process sibling interaction parameter
  s <- as.numeric(s)
  if (length(s) == 1) s <- rep(s, 2)
  if (length(s) != 2) stop("s must be a scalar or length-2 numeric vector.")

  # resolve inputs
  A_diag <- as.numeric(A_diag)
  C_diag <- as.numeric(C_diag)
  E_diag <- as.numeric(E_diag)
  mean   <- as.numeric(mean)

  if (length(A_diag) != 2) stop("A_diag must be length 2.")
  if (length(C_diag) != 2) stop("C_diag must be length 2.")
  if (length(E_diag) != 2) stop("E_diag must be length 2.")
  if (length(mean)   != 2) stop("mean must be length 2.")

  # trait-wise normalization: each trait's latent variance sums to 1
  tot <- A_diag + C_diag + E_diag
  if (any(tot <= 0)) stop("A_diag + C_diag + E_diag must be positive for each trait.")
  A_diag <- A_diag / tot
  C_diag <- C_diag / tot
  E_diag <- E_diag / tot

  # correlation matrices (bivariate: one off-diagonal per matrix)
  rA_eff <- if (!is.null(rA12)) rA12 else rA
  rC_eff <- if (!is.null(rC12)) rC12 else rC
  rE_eff <- if (!is.null(rE12)) rE12 else rE

  R_A <- matrix(c(1, rA_eff, rA_eff, 1), 2, 2)
  R_C <- matrix(c(1, rC_eff, rC_eff, 1), 2, 2)
  R_E <- matrix(c(1, rE_eff, rE_eff, 1), 2, 2)

  # DoC matrix: 2x2 with directional paths b12, b21
  B <- matrix(c(
    0,   b12,
    b21, 0
  ), 2, 2, byrow = TRUE)

  I2 <- diag(2)
  IB <- I2 - B
  if (abs(det(IB)) <= .Machine$double.eps^0.5) stop("I_2 - B is (near-)singular; check b12 and b21.")
  Tm <- solve(IB)

  # latent A, C, E (diagonal SD matrices, not variance matrices)
  D_A <- diag(sqrt(A_diag))
  D_C <- diag(sqrt(C_diag))
  D_E <- diag(sqrt(E_diag))
  A <- D_A %*% R_A %*% D_A
  C <- D_C %*% R_C %*% D_C
  E <- D_E %*% R_E %*% D_E

  # observed-scale (phenotypic) blocks
  V    <- Tm %*% (A + C + E)       %*% t(Tm)
  X_MZ <- Tm %*% (A + C)           %*% t(Tm)
  X_DZ <- Tm %*% (0.5 * A + C)     %*% t(Tm)

  # means on observed scale; phenotype-major 4-vector
  mu_obs <- as.vector(Tm %*% mean)
  mu_4   <- c(mu_obs[1], mu_obs[1], mu_obs[2], mu_obs[2])

  # 4x4 covariance built directly in phenotype-major order:
  # [p1_t1, p1_t2, p2_t1, p2_t2]
  Sigma_phenmajor <- function(V, X) {
    matrix(c(
      V[1,1], X[1,1], V[1,2], X[1,2],
      X[1,1], V[1,1], X[2,1], V[1,2],
      V[2,1], X[1,2], V[2,2], X[2,2],
      X[2,1], V[2,1], X[2,2], V[2,2]
    ), 4, 4, byrow = TRUE)
  }

  Sigma_MZ <- Sigma_phenmajor(V, X_MZ)
  Sigma_DZ <- Sigma_phenmajor(V, X_DZ)

  # sibling interaction: 4x4 in phenotype-major order, coupling twins within trait
  I4   <- diag(4)
  Bsib <- matrix(0, 4, 4)
  # trait 1: positions (1,2) = (p1_t1, p1_t2)
  Bsib[1, 2] <- s[1]
  Bsib[2, 1] <- s[1]
  # trait 2: positions (3,4) = (p2_t1, p2_t2)
  Bsib[3, 4] <- s[2]
  Bsib[4, 3] <- s[2]

  if (any(s != 0)) {
    I4_Bsib <- I4 - Bsib
    if (abs(det(I4_Bsib)) <= .Machine$double.eps^0.5) stop("I_4 - B_sib is (near-)singular; check s1 and s2.")
    G_sib <- solve(I4_Bsib)

    # transform covariances and means
    Sigma_MZ <- G_sib %*% Sigma_MZ %*% t(G_sib)
    Sigma_DZ <- G_sib %*% Sigma_DZ %*% t(G_sib)
    mu_4     <- as.vector(G_sib %*% mu_4)
  } else {
    G_sib <- diag(4)
  }

  # symmetrize (guard against tiny numerical asymmetry from matrix products)
  Sigma_MZ <- (Sigma_MZ + t(Sigma_MZ)) / 2
  Sigma_DZ <- (Sigma_DZ + t(Sigma_DZ)) / 2

  # positive-definiteness checks before drawing
  check_pd <- function(S, label) {
    if (!isSymmetric(S, tol = 1e-8)) stop(sprintf("%s is not symmetric.", label))
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev <= -1e-8)) stop(sprintf("%s is not positive semi-definite (min eigenvalue = %g).",
                                       label, min(ev)))
  }
  check_pd(Sigma_MZ, "Sigma_MZ")
  check_pd(Sigma_DZ, "Sigma_DZ")

  # draws
  draw <- function(n, mu4, Sigma4) {
    if (n <= 0) return(NULL)
    mvtnorm::rmvnorm(n, mean = mu4, sigma = Sigma4)
  }
  to_df <- function(X, z) {
    if (is.null(X)) return(NULL)
    df <- data.frame(X)
    # names in phenotype-major order (matching the underlying simulation order)
    names(df) <- c("p1_t1", "p1_t2", "p2_t1", "p2_t2")
    df$zyg <- z
    df
  }

  Y_MZ <- draw(n_MZ, mu_4, Sigma_MZ)
  Y_DZ <- draw(n_DZ, mu_4, Sigma_DZ)

  dat <- do.call(rbind, list(
    to_df(Y_MZ, 1L),
    to_df(Y_DZ, 2L)
  ))
  dat <- dat[, c("p1_t1", "p1_t2", "p2_t1", "p2_t2", "zyg")]

  # truth object with bivariate parameters
  truth <- list(
    flags = list(),
    varparam = list(
      va11 = A[1,1], va22 = A[2,2], va12 = A[1,2],
      vc11 = C[1,1], vc22 = C[2,2], vc12 = C[1,2],
      ve11 = E[1,1], ve22 = E[2,2], ve12 = E[1,2],
      b12  = B[1,2], b21  = B[2,1]
    ),
    realized = list(
      rA12 = R_A[1,2], rC12 = R_C[1,2], rE12 = R_E[1,2],
      b12  = B[1,2],   b21  = B[2,1],
      s1   = s[1],     s2   = s[2]
    ),
    matrices = list(
      A = A, C = C, E = E,
      B = B, Tm = Tm,
      RA = R_A, RC = R_C, RE = R_E,
      V = V,
      X_MZ = X_MZ, X_DZ = X_DZ,
      Sigma_MZ = Sigma_MZ, Sigma_DZ = Sigma_DZ,
      Bsib = Bsib, G_sib = G_sib
    )
  )

  list(data = dat, truth = truth)
}
