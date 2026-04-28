# ------------------------------------------------------------------------------
# Program: AllSourcePipeline.R
# Author: Matthew Roth
# Date: October 9, 2025
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1) Libraries and options
# ------------------------------------------------------------------------------
require(OpenMx)
require(psych)
require(mvtnorm)  # For multivariate normal random number generation
options(width = 245)
mxVersion()
mxOption(NULL, "Default optimizer", "NPSOL")

# ------------------------------------------------------------------------------
# Data simulation function for univariate models
# ------------------------------------------------------------------------------
source('models/univariate_sim.R')

# ------------------------------------------------------------------------------
# Single Model Univariate Selection Algorithm
# ------------------------------------------------------------------------------

single_model_select_sib <- function(fitACE, fitACES, fits_ACE, fits_ACES, alpha = 0.05, silent = TRUE) {
  # --- helpers (kept inside function for portability) ---
  get_aic <- function(m) {
    tryCatch({
      ic <- summary(m)$informationCriteria
      if (is.data.frame(ic)) {
        as.numeric(ic$par[1])
      } else if (is.matrix(ic)) {
        as.numeric(ic[1, 1])
      } else {
        as.numeric(ic[1])
      }
    }, error = function(e) NA_real_)
  }
  
  safe_mxCompare <- function(m1, m2) {
    tryCatch(as.data.frame(mxCompare(m1, m2)), error = function(e) NULL)
  }
  
  safe_pval_extract <- function(cmp_df, row_idx) {
    tryCatch({
      pval <- as.numeric(cmp_df[row_idx, "p"])
      if (is.finite(pval)) pval else NA_real_
    }, error = function(e) NA_real_)
  }
  
  # --- 0) saturated model ---
  sat0 <- mxRefModels(fitACE)$Saturated
  sat  <- mxRun(sat0, silent = TRUE)
  
  # --- 1) absolute fit screening: ACE vs ACES ---
  cmp_ace <- safe_mxCompare(sat, fitACE)
  cmp_aces <- safe_mxCompare(sat, fitACES)
  
  compare_errors <- list()
  
  if (is.null(cmp_ace) || nrow(cmp_ace) < 2) {
    compare_errors$ACE_vs_sat <- "mxCompare(sat, fitACE) failed"
    ace_pass <- NA
    ace_aic <- NA
    ace_p <- NA
    ace_df <- NA
    ace_diffLL <- NA
  } else {
    ace_aic <- get_aic(fitACE)
    ace_p <- safe_pval_extract(cmp_ace, 2)
    ace_df <- as.numeric(cmp_ace[2, "df"])
    ace_diffLL <- as.numeric(cmp_ace[2, "diffLL"])
    ace_pass <- is.finite(ace_p) && ace_p >= alpha
  }
  
  if (is.null(cmp_aces) || nrow(cmp_aces) < 2) {
    compare_errors$ACES_vs_sat <- "mxCompare(sat, fitACES) failed"
    aces_pass <- NA
    aces_aic <- NA
    aces_p <- NA
    aces_df <- NA
    aces_diffLL <- NA
  } else {
    aces_aic <- get_aic(fitACES)
    aces_p <- safe_pval_extract(cmp_aces, 2)
    aces_df <- as.numeric(cmp_aces[2, "df"])
    aces_diffLL <- as.numeric(cmp_aces[2, "diffLL"])
    aces_pass <- is.finite(aces_p) && aces_p >= alpha
  }
  
  abs_fit_table <- data.frame(
    model = c("ACE", "ACES"),
    AIC = c(ace_aic, aces_aic),
    p = c(ace_p, aces_p),
    df = c(ace_df, aces_df),
    diffLL = c(ace_diffLL, aces_diffLL),
    pass = c(ace_pass, aces_pass),
    stringsAsFactors = FALSE
  )
  
  # --- 2) family choice logic ---
  n_pass <- sum(c(ace_pass, aces_pass), na.rm = TRUE)
  
  if (n_pass == 1 && ace_pass) {
    chosen_family <- "ACE"
    family_reason <- "only_ACE_passed"
    full_model <- fitACE
    fits_chosen <- fits_ACE
  } else if (n_pass == 1 && aces_pass) {
    chosen_family <- "ACES"
    family_reason <- "only_ACES_passed"
    full_model <- fitACES
    fits_chosen <- fits_ACES
  } else if (n_pass == 2) {
    # both passed: choose by AIC
    if (ace_aic < aces_aic) {
      chosen_family <- "ACE"
      family_reason <- "both_passed_AIC"
      full_model <- fitACE
      fits_chosen <- fits_ACE
    } else {
      chosen_family <- "ACES"
      family_reason <- "both_passed_AIC"
      full_model <- fitACES
      fits_chosen <- fits_ACES
    }
  } else {
    # both failed or both NA: choose by AIC
    if (is.finite(ace_aic) && is.finite(aces_aic)) {
      if (ace_aic < aces_aic) {
        chosen_family <- "ACE"
        family_reason <- "both_failed_AIC"
        full_model <- fitACE
        fits_chosen <- fits_ACE
      } else {
        chosen_family <- "ACES"
        family_reason <- "both_failed_AIC"
        full_model <- fitACES
        fits_chosen <- fits_ACES
      }
    } else if (is.finite(ace_aic)) {
      chosen_family <- "ACE"
      family_reason <- "both_failed_AIC"
      full_model <- fitACE
      fits_chosen <- fits_ACE
    } else {
      chosen_family <- "ACES"
      family_reason <- "both_failed_AIC"
      full_model <- fitACES
      fits_chosen <- fits_ACES
    }
  }
  
  if (!silent) {
    cat("Chosen family:", chosen_family, "\n")
    cat("Family reason:", family_reason, "\n")
  }
  
  # --- 3) nested reduction within chosen family ---
  reduction_table <- data.frame(
    model = character(),
    AIC = numeric(),
    p = numeric(),
    df = numeric(),
    diffLL = numeric(),
    pass = logical(),
    compare_status = character(),
    stringsAsFactors = FALSE
  )
  
  retained_models <- character(0)
  
  for (child_nm in names(fits_chosen)) {
    child_model <- fits_chosen[[child_nm]]
    cmp <- safe_mxCompare(full_model, child_model)
    
    if (is.null(cmp) || nrow(cmp) < 2) {
      child_aic <- get_aic(child_model)
      reduction_table <- rbind(reduction_table, data.frame(
        model = child_nm,
        AIC = child_aic,
        p = NA_real_,
        df = NA_real_,
        diffLL = NA_real_,
        pass = NA,
        compare_status = "compare_error",
        stringsAsFactors = FALSE
      ))
      if (length(compare_errors) == 0 || !("reduction_compare_errors" %in% names(compare_errors))) {
        compare_errors[[paste0("reduction_", child_nm)]] <- "mxCompare failed"
      }
      next
    }
    
    pval <- safe_pval_extract(cmp, 2)
    dfv <- as.numeric(cmp[2, "df"])
    dll <- as.numeric(cmp[2, "diffLL"])
    child_aic <- get_aic(child_model)
    keep_child <- is.finite(pval) && pval >= alpha
    
    reduction_table <- rbind(reduction_table, data.frame(
      model = child_nm,
      AIC = child_aic,
      p = pval,
      df = dfv,
      diffLL = dll,
      pass = keep_child,
      compare_status = if (keep_child) "retained" else "rejected",
      stringsAsFactors = FALSE
    ))
    
    if (keep_child) {
      retained_models <- c(retained_models, child_nm)
    }
  }
  
  if (!silent) {
    cat("Reduced models retained:", 
        if (length(retained_models) > 0) paste(retained_models, collapse = ", ") else "none", "\n")
  }
  
  # --- 4) final single-model choice ---
  if (length(retained_models) > 0) {
    # multiple or single retained: choose by lowest AIC among retained
    retained_aics <- setNames(
      reduction_table$AIC[reduction_table$model %in% retained_models],
      retained_models
    )
    selected_name <- names(which.min(retained_aics))
    single_fit <- fits_chosen[[selected_name]]
    selection_method <- if (length(retained_models) == 1) "single_retained" else "lowest_AIC_among_retained"
  } else {
    # none retained: fallback to full model
    selected_name <- chosen_family
    single_fit <- full_model
    selection_method <- "fallback_to_full"
  }
  
  if (!silent) {
    cat("Selected model:", selected_name, "\n")
    cat("Selection method:", selection_method, "\n")
  }
  
  list(
    single_fit = single_fit,
    selected_name = selected_name,
    chosen_family = chosen_family,
    family_reason = family_reason,
    abs_fit_table = abs_fit_table,
    reduction_table = reduction_table,
    compare_errors = if (length(compare_errors) > 0) compare_errors else NULL
  )
}



# ------------------------------------------------------------------------------
# Univariate Model Fuction
# ------------------------------------------------------------------------------
fit_univariate <- function(
    data,
    topName = "uni",
    selVars = c("p1_t1","p1_t2"),
    nv = 1,
    lboundVar = 1e-5,  
    
    means = list(me1 = c(TRUE, 0)),
    
    A_diag = list(a11 = c(FALSE, 0)),
    C_diag = list(c11 = c(FALSE, 0)),
    D_diag = list(d11 = c(FALSE, 0)),
    E_diag = list(e11 = c(FALSE, 0)),
    
    S = list(
      s1 = c(FALSE, 0)
    ),
    
    ci = character(0),
    tryHard = TRUE,
    extraTries = 15,
    intervals = TRUE
) {
  stopifnot(nv == 1)
  
  dataMZ <- subset(data, zyg==1)
  dataDZ <- subset(data, zyg==2)
  
  get_fs <- function(x, nm) {
    if (is.null(x) || length(x) != 2) stop("Parameter '", nm, "' must be c(free, start).")
    list(free = as.logical(x[[1]]), start = as.numeric(x[[2]]))
  }
  
  # ----------------------------
  # Means
  # ----------------------------
  me1 <- get_fs(means$me1, "me1")
  
  M <- mxMatrix("Full", 1, 1,
                free   = me1$free,
                values = me1$start,
                labels = "me1",
                name   = "M")
  
  expMean <- mxAlgebra(
    cbind(M[1,1], M[1,1]),
    name = "expMean"
  )
  
  # ----------------------------
  # Sibling interaction (s1)
  # ----------------------------
  if (is.null(S$s1)) S$s1 <- c(FALSE, 0)
  s1 <- get_fs(S$s1, "s1")
  
  Bsib_free   <- matrix(FALSE, 2, 2)
  Bsib_values <- matrix(0,     2, 2)
  Bsib_labels <- matrix(NA_character_, 2, 2)
  
  setSib <- function(i, j, spec, nm) {
    Bsib_free[i,j]   <<- spec$free
    Bsib_values[i,j] <<- spec$start
    if (spec$free) Bsib_labels[i,j] <<- nm
  }
  
  setSib(1, 2, s1, "s1")
  setSib(2, 1, s1, "s1")
  
  Bsib <- mxMatrix("Full", 2, 2,
                   free   = as.vector(t(Bsib_free)),
                   values = as.vector(t(Bsib_values)),
                   labels = as.vector(t(Bsib_labels)),
                   lbound = as.vector(t(ifelse(Bsib_free, -.50, NA))),
                   ubound = as.vector(t(ifelse(Bsib_free,  .50, NA))),
                   byrow  = TRUE,
                   name   = "Bsib")
  
  I2 <- mxMatrix("Iden", 2, 2, name="I2")
  G_sib <- mxAlgebra(solve(I2 - Bsib), name="G_sib")
  
  # ----------------------------
  # Variance components
  # ----------------------------
  a11 <- get_fs(A_diag$a11, "A_diag$a11")
  c11 <- get_fs(C_diag$c11, "C_diag$c11")
  d11 <- get_fs(D_diag$d11, "D_diag$d11")
  e11 <- get_fs(E_diag$e11, "E_diag$e11")
  
  mkDiagVar <- function(spec, label, name) {
    free  <- spec$free
    val_v <- spec$start
    if (!is.finite(val_v)) stop(name, " start variance must be finite.")
    
    lab <- if (free) label else NA_character_
    
    
    if (free && !is.na(lboundVar)) {
      mxMatrix("Diag", 1, 1,
               free   = free,
               values = val_v,
               lbound = lboundVar,
               labels = lab,
               name   = name)
    } else {
      mxMatrix("Diag", 1, 1,
               free   = free,
               values = val_v,
               labels = lab,
               name   = name)
    }
  }
  
  aV <- mkDiagVar(a11, "va11", "aV")
  cV <- mkDiagVar(c11, "vc11", "cV")
  dV <- mkDiagVar(d11, "vd11", "dV")
  eV <- mkDiagVar(e11, "ve11", "eV")
  
  A <- mxAlgebra(aV, name="A")
  C <- mxAlgebra(cV, name="C")
  D <- mxAlgebra(dV, name="D")
  E <- mxAlgebra(eV, name="E")
  
  V   <- mxAlgebra(A + C + D + E, name="V")
  cMZ <- mxAlgebra(A + C + D, name="cMZ")
  cDZ <- mxAlgebra(0.5 * A + C + 0.25 * D, name="cDZ")
  
  expCovMZ <- mxAlgebra(rbind(
    cbind( V[1,1],   cMZ[1,1] ),
    cbind( cMZ[1,1], V[1,1]   )
  ), name="expCovMZ")
  
  expCovDZ <- mxAlgebra(rbind(
    cbind( V[1,1],   cDZ[1,1] ),
    cbind( cDZ[1,1], V[1,1]   )
  ), name="expCovDZ")
  
  expCovMZ_sib <- mxAlgebra(G_sib %*% expCovMZ %*% t(G_sib), name="expCovMZ_sib")
  expCovDZ_sib <- mxAlgebra(G_sib %*% expCovDZ %*% t(G_sib), name="expCovDZ_sib")
  expMean_sib  <- mxAlgebra(expMean %*% t(G_sib), name="expMean_sib")
  
  subMZ <- mxModel("MZ",
                   mxData(dataMZ[, selVars, drop=FALSE], type="raw"),
                   mxExpectationNormal(paste0(topName,".expCovMZ_sib"),
                                       paste0(topName,".expMean_sib"),
                                       dimnames=selVars),
                   mxFitFunctionML())
  
  subDZ <- mxModel("DZ",
                   mxData(dataDZ[, selVars, drop=FALSE], type="raw"),
                   mxExpectationNormal(paste0(topName,".expCovDZ_sib"),
                                       paste0(topName,".expMean_sib"),
                                       dimnames=selVars),
                   mxFitFunctionML())
  
  top <- mxModel(topName,
                 M, expMean,
                 Bsib, I2, G_sib,
                 aV, cV, dV, eV,
                 A, C, D, E,
                 V, cMZ, cDZ,
                 expCovMZ, expCovDZ,
                 expCovMZ_sib, expCovDZ_sib, expMean_sib,
                 subMZ, subDZ,
                 mxFitFunctionMultigroup(c("MZ","DZ")))
  
  if (length(ci)) top <- mxModel(top, mxCI(ci))
  
  mxOption(NULL, "Calculate Hessian", "Yes")
  mxOption(NULL, "Standard Errors",  "Yes")
  
  if (isTRUE(tryHard)) {
    fit <- mxTryHard(top, greenOK=TRUE, extraTries=extraTries, intervals=intervals)
  } else {
    fit <- mxRun(top, intervals=intervals)
  }
  
  fit
}





# ------------------------------------------------------------------------------
# Attach standardized-proportion algebra for mxModelAverage
# ------------------------------------------------------------------------------

add_prop_block_uni <- function(fit) {
  add_if_missing <- function(model, name, algebra_or_matrix) {
    if (!(name %in% names(model$algebras)) && !(name %in% names(model$matrices))) {
      model <- mxModel(model, algebra_or_matrix)
    }
    model
  }
  
  # Require A, E at the top level
  req <- c("A", "E")
  missing_req <- setdiff(req, c(names(fit$algebras), names(fit$matrices)))
  if (length(missing_req)) stop("downstream: missing required objects: ", paste(missing_req, collapse=", "))
  
  # Ensure C exists; if absent, add 1x1 Zero
  if (!("C" %in% names(fit$algebras) || "C" %in% names(fit$matrices))) {
    fit <- mxModel(fit, mxMatrix(type="Zero", nrow=1, ncol=1, name="C"))
  }
  
  # Ensure D exists; if absent, add 1x1 Zero
  if (!("D" %in% names(fit$algebras) || "D" %in% names(fit$matrices))) {
    fit <- mxModel(fit, mxMatrix(type="Zero", nrow=1, ncol=1, name="D"))
  }
  
  # Ensure V exists; if absent, build from components
  hasV <- "V" %in% c(names(fit$algebras), names(fit$matrices))
  if (!hasV) {
    fit <- mxModel(fit, mxAlgebra(A + C + D + E, name="V"))
  }
  
  # Inverse SD matrix from phenotypic V
  fit <- add_if_missing(fit, "iSD", mxAlgebra(vec2diag(1 / sqrt(diag2vec(V))), name="iSD"))
  
  # Proportions for univariate: Pva11, Pvd11, Pvc11, Pve11, totV
  if (!("Prop" %in% names(fit$algebras))) {
    estProp <- mxAlgebra(
      cbind(
        A[1,1]/V[1,1],
        D[1,1]/V[1,1],
        C[1,1]/V[1,1],
        E[1,1]/V[1,1],
        V[1,1]
      ),
      name="Prop",
      dimnames=list("pro", c("Pva11","Pvd11","Pvc11","Pve11","totV"))
    )
    fit <- mxModel(fit, estProp)
  }
  
  fit
}




# Function to flatten nested lists and create column names with model prefix
flatten_output <- function(output_list, model_name) {
  result <- list()
  
  # Extract matrices
  if (!is.null(output_list$matrices)) {
    for (mat_name in names(output_list$matrices)) {
      mat <- output_list$matrices[[mat_name]]
      # Flatten matrix to vector with indexed column names
      for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
          col_name <- paste0(model_name, "_", mat_name, "_", i, "_", j)
          result[[col_name]] <- mat[i, j]
        }
      }
    }
  }
  
  # Extract algebras (excluding matrices and fit functions with attributes)
  if (!is.null(output_list$algebras)) {
    for (alg_name in names(output_list$algebras)) {
      alg <- output_list$algebras[[alg_name]]
      
      # Skip complex objects with attributes
      if (alg_name %in% c("MZ.fitfunction", "DZ.fitfunction")) {
        result[[paste0(model_name, "_", alg_name, "_value")]] <- as.numeric(alg)
        next
      }
      
      # Handle matrices
      if (is.matrix(alg)) {
        for (i in 1:nrow(alg)) {
          for (j in 1:ncol(alg)) {
            col_name <- paste0(model_name, "_", alg_name, "_", i, "_", j)
            result[[col_name]] <- alg[i, j]
          }
        }
      } else if (is.numeric(alg) && length(alg) == 1) {
        result[[paste0(model_name, "_", alg_name)]] <- alg
      }
    }
  }
  
  # Extract key scalar values
  scalar_fields <- c("fit", "Minus2LogLikelihood", "minimum", 
                     "maxRelativeOrdinalError", "iterations", "evaluations")
  
  for (field in scalar_fields) {
    if (!is.null(output_list[[field]])) {
      result[[paste0(model_name, "_", field)]] <- output_list[[field]]
    }
  }
  
  # Extract status
  if (!is.null(output_list$status)) {
    result[[paste0(model_name, "_status_code")]] <- output_list$status$code
    result[[paste0(model_name, "_status_status")]] <- output_list$status$status
  }
  
  # Extract timing information
  if (!is.null(output_list$wallTime)) {
    result[[paste0(model_name, "_wallTime_seconds")]] <- as.numeric(output_list$wallTime, units = "secs")
  }
  
  return(result)
}

# Main function to process all models into a single row
process_all_models <- function(fits) {
  # Extract output for each model and combine into single list
  combined_results <- list()
  
  for (model_name in names(fits)) {
    model_output <- fits[[model_name]]$output
    flattened <- flatten_output(model_output, model_name)
    # Merge into combined results
    combined_results <- c(combined_results, flattened)
  }
  
  # Convert to single-row dataframe
  df <- as.data.frame(combined_results, stringsAsFactors = FALSE)
  
  return(df)
}




run_simulation_and_analysis <- function(
    n_MZ = 500, n_DZ = 500,
    seed = 1L,
    output_prefix = "",
    s = 0,
    A_diag = .5, C_diag = .3, D_diag = 0, E_diag = 0.2,
    mean = 0
) {
  
  if (identical(output_prefix, "")) {
    output_prefix <- sprintf("simulation_results_(A=0.70, C=.20, E=0.10, S=0.30)/run_%04d_", as.integer(seed))
    dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
  }
  
  
  
  safe_div <- function(x, y) {
    out <- rep(NA_real_, length.out = max(length(x), length(y)))
    ok <- is.finite(x) & is.finite(y) & (y != 0)
    out[ok] <- x[ok] / y[ok]
    out
  }
  
  mat2 <- function(a11, a22, a12 = 0) matrix(c(a11, a12, a12, a22), nrow = 2, byrow = TRUE)
  
  
  sim_results <- simulate_twins(
    n_MZ = 1000,
    n_DZ = 1000,
    seed = seed,
    s = 0.3 + runif(1, -0.01, 0.01),
    A_diag = pmax(0.70 + runif(1, -0.01, 0.01), 0),
    C_diag = pmax(0.20 + runif(1, -0.01, 0.01), 0),
    D_diag = pmax(0 + runif(1, 0, 0), 0),
    E_diag = pmax(0.10 + runif(1, -0.01, 0.01), 0),
    mean = 0
  )
  
  
  
  sim_data <- sim_results$data
  truth    <- sim_results$truth
  
  vars <- c("p1_t1","p1_t2")
  
  # Covariance matrices by zygosity
  S_mz <- cov(sim_data[sim_data$zyg==1, vars], use="pairwise.complete.obs")
  S_dz <- cov(sim_data[sim_data$zyg==2, vars], use="pairwise.complete.obs")
  
  utils::write.csv(sim_data, paste0(output_prefix, "sim_data.csv"), row.names = FALSE)
  # ========================================================================
  # Precomputation of Start Values
  # ========================================================================
  
  vars4 <- c("p1_t1","p1_t2")
  
  dat_all <- sim_data[, vars4]
  dat_mz  <- sim_data[sim_data$zyg == 1, vars4, drop=FALSE]
  dat_dz  <- sim_data[sim_data$zyg == 2, vars4, drop=FALSE]
  
  Z_all <- scale(dat_all)
  Z_mz  <- scale(dat_mz)
  Z_dz  <- scale(dat_dz)
  
  cov_all <- cov(Z_all, use="pairwise.complete.obs")
  cov_mz  <- cov(Z_mz,  use="pairwise.complete.obs")
  cov_dz  <- cov(Z_dz,  use="pairwise.complete.obs")
  
  clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  
  
  rMZ_p1 <- cov_mz["p1_t1","p1_t2"]
  rDZ_p1 <- cov_dz["p1_t1","p1_t2"]
  
  h2_1 <- clip(2 * (rMZ_p1 - rDZ_p1), 0.05, 0.95)
  c2_1 <- clip(rMZ_p1 - h2_1,         0.00, 0.90)
  e2_1 <- clip(1 - h2_1 - c2_1,       0.05, 0.95)
  a2_1 <- clip(4 * rDZ_p1 - rMZ_p1,       0.05, 0.95)
  
  
  sv_va11 <- h2_1; sv_vc11 <- c2_1; sv_ve11 <- e2_1
  
  sib_cross_all <- cov_all["p1_t1","p1_t2"]
  sv_s1 <- clip(sib_cross_all - rDZ_p1, -0.95, 0.95)
  
  starts <- c(
    sv_va11 = sv_va11, sv_vc11 = sv_vc11, sv_ve11 = sv_ve11,
    sv_s1   = sv_s1
  )
  
  # ========================================================================
  # METHOD 0 – Quality Control if DGP were known
  # ========================================================================
  
  fitACE   <- fit_univariate( data = sim_data, topName = "fitACE",
                              A_diag = list(a11 = c(TRUE, sv_va11)),
                              C_diag = list(c11 = c(TRUE, sv_vc11)),
                              E_diag = list(e11 = c(TRUE, sv_ve11))
  )
  fitAE   <- fit_univariate( data = sim_data, topName = "fitAE",
                             A_diag = list(a11 = c(TRUE, sv_va11)),
                             E_diag = list(e11 = c(TRUE, sv_ve11))
  )
  fitCE   <- fit_univariate( data = sim_data, topName = "fitCE",
                             C_diag = list(c11 = c(TRUE, sv_vc11)),
                             E_diag = list(e11 = c(TRUE, sv_ve11))
  )
  fitE   <- fit_univariate( data = sim_data, topName = "fitE",
                            E_diag = list(e11 = c(TRUE, sv_ve11))
  )
  
  
  fitACES   <- fit_univariate( data = sim_data, topName = "fitACES",
                               A_diag = list(a11 = c(TRUE, sv_va11)),
                               C_diag = list(c11 = c(TRUE, sv_vc11)),
                               E_diag = list(e11 = c(TRUE, sv_ve11)),
                               S = list(s1 = c(TRUE, sv_s1))
  )
  fitAES   <- fit_univariate( data = sim_data, topName = "fitAES",
                              A_diag = list(a11 = c(TRUE, sv_va11)),
                              E_diag = list(e11 = c(TRUE, sv_ve11)),
                              S = list(s1 = c(TRUE, sv_s1))
  )
  fitES   <- fit_univariate( data = sim_data, topName = "fitES",
                             E_diag = list(e11 = c(TRUE, sv_ve11)),
                             S = list(s1 = c(TRUE, sv_s1))
  )
  
  fits <- list(
    ACES = fitACES
  )
  
  single_fit <- fits[[1]]
  results_df <- process_all_models(fits)
  
  write.csv(results_df,paste0(output_prefix, "models.csv"), row.names = FALSE)
  
  use_vp <- if (isTRUE(truth$flags$sex_differences)) truth$varparam_m else truth$varparam
  if (is.null(use_vp)) use_vp <- truth$varparam
  
  truth_raw <- c(
    va11 = use_vp$va11,
    vc11 = if (!is.null(use_vp$vc11)) use_vp$vc11 else 0,
    ve11 = use_vp$ve11
  )
  
  totV_true <- truth_raw["va11"] + truth_raw["vc11"]  + truth_raw["ve11"]
  
  Pva11_true <- safe_div(truth_raw["va11"], totV_true)
  Pvc11_true <- safe_div(truth_raw["vc11"], totV_true)
  Pve11_true <- safe_div(truth_raw["ve11"], totV_true)
  
  s1_truth <- NA_real_
  if (!is.null(truth$realized) && !is.null(truth$realized$s1)) {
    s1_truth <- truth$realized$s1
  } else if (!is.null(truth$matrices$Bsib) && nrow(truth$matrices$Bsib) >= 1 && ncol(truth$matrices$Bsib) >= 2) {
    s1_truth <- truth$matrices$Bsib[1, 2]
  }
  
  truth_raw["s1"] <- s1_truth
  
  truth_std <- c(
    Pva11 = Pva11_true,
    Pvc11 = Pvc11_true,
    Pve11 = Pve11_true,
    totV  = totV_true
  )
  
  map_prop_names <- c(
    "Prop[1,1]"="Pva11",
    "Prop[1,2]"="Pvc11",
    "Prop[1,3]"="Pve11",
    "Prop[1,4]"="totV"
  )
  
  # Extract fitted parameters
  fit_params <- single_fit$output$estimate
  coef_vals <- coef(single_fit)
  vcov_mat <- vcov(single_fit)
  
  ### SIBLING INTERACTION ADD-ON ###
  s1_est <- if ("s1" %in% names(coef_vals)) coef_vals["s1"] else 0
  ### END SIBLING INTERACTION ADD-ON ###
  
  # Helper function to compute Prop vector from parameter vector
  compute_prop_from_params <- function(param_vals) {
    # Extract variance components
    va11_p <- if ("va11" %in% names(param_vals)) param_vals["va11"] else NA_real_
    vc11_p <- if ("vc11" %in% names(param_vals)) param_vals["vc11"] else 0
    ve11_p <- if ("ve11" %in% names(param_vals)) param_vals["ve11"] else NA_real_
    
    va11_p <- if (va11_p > 0) va11_p else 0
    vc11_p <- if (vc11_p > 0) vc11_p else 0
    ve11_p <- if (ve11_p > 0) ve11_p else 0
    
    # Compute total variance
    totV_p <- va11_p + vc11_p  + ve11_p
    
    # Compute proportions
    Pva11_p <- safe_div(va11_p, totV_p)
    Pvc11_p <- safe_div(vc11_p, totV_p)
    Pve11_p <- safe_div(ve11_p, totV_p)
    
    c(Pva11_p, Pvc11_p, Pve11_p, totV_p)
  }
  
  # Compute Prop vector
  prop_vec <- compute_prop_from_params(fit_params)
  
  # Compute delta-method SEs
  eps_base <- 1e-6
  
  if (is.null(coef_vals) || length(coef_vals) == 0 || 
      is.null(vcov_mat) || !all(is.finite(vcov_mat))) {
    prop_ses <- rep(NA_real_, length(prop_vec))
  } else {
    n_params <- length(coef_vals)
    n_props <- length(prop_vec)
    gradients <- matrix(0, nrow = n_props, ncol = n_params)
    
    for (j in seq_len(n_params)) {
      pj <- coef_vals[j]
      h <- eps_base * (abs(pj) + 1)
      
      # Forward perturbation
      coef_plus <- coef_vals
      coef_plus[j] <- pj + h
      prop_plus <- tryCatch(compute_prop_from_params(coef_plus), 
                            error = function(e) rep(NA_real_, n_props))
      
      # Backward perturbation
      coef_minus <- coef_vals
      coef_minus[j] <- pj - h
      prop_minus <- tryCatch(compute_prop_from_params(coef_minus), 
                             error = function(e) rep(NA_real_, n_props))
      
      # Central difference gradient
      for (k in seq_len(n_props)) {
        if (is.finite(prop_plus[k]) && is.finite(prop_minus[k])) {
          gradients[k, j] <- (prop_plus[k] - prop_minus[k]) / (2 * h)
        } else {
          gradients[k, j] <- 0
        }
      }
    }
    
    # Compute delta-method SEs: sqrt(g' * vcov * g)
    prop_ses <- apply(gradients, 1, function(g) {
      if (all(is.finite(g)) && any(g != 0)) {
        se_sq <- as.numeric(g %*% vcov_mat %*% g)
        if (se_sq >= 0) sqrt(se_sq) else NA_real_
      } else {
        NA_real_
      }
    })
  }
  
  # Build estimates table
  prop_rownames <- names(map_prop_names)
  estimates_df <- data.frame(
    param    = prop_rownames,
    estimate = prop_vec,
    SE       = prop_ses,
    source   = "standardized",
    stringsAsFactors = FALSE
  )
  estimates_df$lcl95 <- estimates_df$estimate - 1.96 * estimates_df$SE
  estimates_df$ucl95 <- estimates_df$estimate + 1.96 * estimates_df$SE
  
  # Map param names and add scale/truth
  estimates_df$param <- unname(map_prop_names[estimates_df$param])
  
  scale_tag <- function(p) {
    if (p %in% names(truth_raw)) "raw" else "standardized"
  }
  
  estimates_df$scale <- vapply(estimates_df$param, scale_tag, FUN.VALUE = character(1))
  estimates_df$truth <- mapply(function(p, sc) {
    if (sc == "raw") {
      if (p %in% names(truth_raw)) truth_raw[[p]] else NA_real_
    } else {
      if (p %in% names(truth_std)) truth_std[[p]] else NA_real_
    }
  }, estimates_df$param, estimates_df$scale)
  
  ### SIBLING INTERACTION ADD-ON ###
  s1_se <- if (!is.null(vcov_mat) && "s1" %in% rownames(vcov_mat)) {
    sqrt(vcov_mat["s1", "s1"])
  } else {
    NA_real_
  }
  
  sib_df <- data.frame(
    param    = "s1",
    estimate = s1_est,
    SE       = s1_se,
    source   = "raw",
    lcl95    = s1_est - 1.96 * s1_se,
    ucl95    = s1_est + 1.96 * s1_se,
    scale    = "raw",
    truth    = s1_truth,
    stringsAsFactors = FALSE
  )
  
  # Bind sibling row to estimates_df
  estimates_df <- rbind(estimates_df, sib_df)
  rownames(estimates_df) <- NULL
  ### END SIBLING INTERACTION ADD-ON ###
  
  # Compute performance metrics
  performance_metrics <- within(estimates_df, {
    bias      <- estimate - truth
    MAE       <- abs(bias)
    MSE       <- bias^2
    covered95 <- ifelse(!is.na(truth) & (truth >= lcl95) & (truth <= ucl95), 1L, 0L)
    typeI     <- ifelse(!is.na(truth) & truth == 0 & (lcl95 > 0 | ucl95 < 0), 1L, 0L)
    typeII    <- ifelse(!is.na(truth) & truth != 0 & (lcl95 <= 0 & ucl95 >= 0), 1L, 0L)
  })
  
  # Write CSV outputs
  utils::write.csv(estimates_df,
                   paste0(output_prefix, "model_estimates_m0", ".csv"),
                   row.names = FALSE)
  
  utils::write.csv(performance_metrics,
                   paste0(output_prefix, "performance_metrics_m0.csv"),
                   row.names = FALSE)
  
  # Return results
  res_m0 <- list(estimates = estimates_df, metrics = performance_metrics)
  
  # ========================================================================
  # METHOD 1 –  Single Model Selection  
  # ========================================================================
  
  fits_ACE <- list(ACE = fitACE, AE = fitAE, CE = fitCE, E = fitE)
  
  fits_ACES <- list(ACES = fitACES, AES = fitAES, ES = fitES)
  
  
  fits <- list(fit = single_model_select_sib(fitACE, fitACES, fits_ACE, fits_ACES)$single_fit)
  
  single_fit <- fits[[1]]
  results_df <- process_all_models(fits)
  
  use_vp <- if (isTRUE(truth$flags$sex_differences)) truth$varparam_m else truth$varparam
  if (is.null(use_vp)) use_vp <- truth$varparam
  
  truth_raw <- c(
    va11 = use_vp$va11,
    vc11 = if (!is.null(use_vp$vc11)) use_vp$vc11 else 0,
    ve11 = use_vp$ve11
  )
  
  totV_true <- truth_raw["va11"] + truth_raw["vc11"]  + truth_raw["ve11"]
  
  Pva11_true <- safe_div(truth_raw["va11"], totV_true)
  Pvc11_true <- safe_div(truth_raw["vc11"], totV_true)
  Pve11_true <- safe_div(truth_raw["ve11"], totV_true)
  
  s1_truth <- NA_real_
  if (!is.null(truth$realized) && !is.null(truth$realized$s1)) {
    s1_truth <- truth$realized$s1
  } else if (!is.null(truth$matrices$Bsib) && nrow(truth$matrices$Bsib) >= 1 && ncol(truth$matrices$Bsib) >= 2) {
    s1_truth <- truth$matrices$Bsib[1, 2]
  }
  
  truth_raw["s1"] <- s1_truth
  
  truth_std <- c(
    Pva11 = Pva11_true,
    Pvc11 = Pvc11_true,
    Pve11 = Pve11_true,
    totV  = totV_true
  )
  
  map_prop_names <- c(
    "Prop[1,1]"="Pva11",
    "Prop[1,2]"="Pvc11",
    "Prop[1,3]"="Pve11",
    "Prop[1,4]"="totV"
  )
  
  # Extract fitted parameters
  fit_params <- single_fit$output$estimate
  coef_vals <- coef(single_fit)
  vcov_mat <- vcov(single_fit)
  
  s1_est <- if ("s1" %in% names(coef_vals)) coef_vals["s1"] else 0

  # Helper function to compute Prop vector from parameter vector
  compute_prop_from_params <- function(param_vals) {
    # Extract variance components
    va11_p <- if ("va11" %in% names(param_vals)) param_vals["va11"] else NA_real_
    vc11_p <- if ("vc11" %in% names(param_vals)) param_vals["vc11"] else 0
    ve11_p <- if ("ve11" %in% names(param_vals)) param_vals["ve11"] else NA_real_
    
    va11_p <- if (va11_p > 0) va11_p else 0
    vc11_p <- if (vc11_p > 0) vc11_p else 0
    ve11_p <- if (ve11_p > 0) ve11_p else 0
    
    # Compute total variance
    totV_p <- va11_p + vc11_p  + ve11_p
    
    # Compute proportions
    Pva11_p <- safe_div(va11_p, totV_p)
    Pvc11_p <- safe_div(vc11_p, totV_p)
    Pve11_p <- safe_div(ve11_p, totV_p)
    
    c(Pva11_p, Pvc11_p, Pve11_p, totV_p)
  }
  
  # Compute Prop vector
  prop_vec <- compute_prop_from_params(fit_params)
  
  # Compute delta-method SEs
  eps_base <- 1e-6
  
  if (is.null(coef_vals) || length(coef_vals) == 0 || 
      is.null(vcov_mat) || !all(is.finite(vcov_mat))) {
    prop_ses <- rep(NA_real_, length(prop_vec))
  } else {
    n_params <- length(coef_vals)
    n_props <- length(prop_vec)
    gradients <- matrix(0, nrow = n_props, ncol = n_params)
    
    for (j in seq_len(n_params)) {
      pj <- coef_vals[j]
      h <- eps_base * (abs(pj) + 1)
      
      # Forward perturbation
      coef_plus <- coef_vals
      coef_plus[j] <- pj + h
      prop_plus <- tryCatch(compute_prop_from_params(coef_plus), 
                            error = function(e) rep(NA_real_, n_props))
      
      # Backward perturbation
      coef_minus <- coef_vals
      coef_minus[j] <- pj - h
      prop_minus <- tryCatch(compute_prop_from_params(coef_minus), 
                             error = function(e) rep(NA_real_, n_props))
      
      # Central difference gradient
      for (k in seq_len(n_props)) {
        if (is.finite(prop_plus[k]) && is.finite(prop_minus[k])) {
          gradients[k, j] <- (prop_plus[k] - prop_minus[k]) / (2 * h)
        } else {
          gradients[k, j] <- 0
        }
      }
    }
    
    # Compute delta-method SEs: sqrt(g' * vcov * g)
    prop_ses <- apply(gradients, 1, function(g) {
      if (all(is.finite(g)) && any(g != 0)) {
        se_sq <- as.numeric(g %*% vcov_mat %*% g)
        if (se_sq >= 0) sqrt(se_sq) else NA_real_
      } else {
        NA_real_
      }
    })
  }
  
  # Build estimates table
  prop_rownames <- names(map_prop_names)
  estimates_df <- data.frame(
    param    = prop_rownames,
    estimate = prop_vec,
    SE       = prop_ses,
    source   = "standardized",
    stringsAsFactors = FALSE
  )
  estimates_df$lcl95 <- estimates_df$estimate - 1.96 * estimates_df$SE
  estimates_df$ucl95 <- estimates_df$estimate + 1.96 * estimates_df$SE
  
  # Map param names and add scale/truth
  estimates_df$param <- unname(map_prop_names[estimates_df$param])
  
  scale_tag <- function(p) {
    if (p %in% names(truth_raw)) "raw" else "standardized"
  }
  
  estimates_df$scale <- vapply(estimates_df$param, scale_tag, FUN.VALUE = character(1))
  estimates_df$truth <- mapply(function(p, sc) {
    if (sc == "raw") {
      if (p %in% names(truth_raw)) truth_raw[[p]] else NA_real_
    } else {
      if (p %in% names(truth_std)) truth_std[[p]] else NA_real_
    }
  }, estimates_df$param, estimates_df$scale)
  
  s1_se <- if (!is.null(vcov_mat) && "s1" %in% rownames(vcov_mat)) {
    sqrt(vcov_mat["s1", "s1"])
  } else {
    NA_real_
  }
  
  sib_df <- data.frame(
    param    = "s1",
    estimate = s1_est,
    SE       = s1_se,
    source   = "raw",
    lcl95    = s1_est - 1.96 * s1_se,
    ucl95    = s1_est + 1.96 * s1_se,
    scale    = "raw",
    truth    = s1_truth,
    stringsAsFactors = FALSE
  )
  
  # Bind sibling row to estimates_df
  estimates_df <- rbind(estimates_df, sib_df)
  rownames(estimates_df) <- NULL
  ### END SIBLING INTERACTION ADD-ON ###
  
  # Compute performance metrics
  performance_metrics <- within(estimates_df, {
    bias      <- estimate - truth
    MAE       <- abs(bias)
    MSE       <- bias^2
    covered95 <- ifelse(!is.na(truth) & (truth >= lcl95) & (truth <= ucl95), 1L, 0L)
    typeI     <- ifelse(!is.na(truth) & truth == 0 & (lcl95 > 0 | ucl95 < 0), 1L, 0L)
    typeII    <- ifelse(!is.na(truth) & truth != 0 & (lcl95 <= 0 & ucl95 >= 0), 1L, 0L)
  })
  
  # Write CSV outputs
  utils::write.csv(estimates_df,
                   paste0(output_prefix, "model_estimates_m1", ".csv"),
                   row.names = FALSE)
  
  utils::write.csv(performance_metrics,
                   paste0(output_prefix, "performance_metrics_m1.csv"),
                   row.names = FALSE)
  
  # Return results
  res_m1 <- list(estimates = estimates_df, metrics = performance_metrics)
  
  
  # ========================================================================
  # METHOD 2 – Univariate Multimodel Averaging
  # ========================================================================
  
  fits <- list(ACE = fitACE,
               AE = fitAE, CE = fitCE, E = fitE,
               ACES = fitACES, 
               AES = fitAES, ES = fitES)
  
  results_df <- process_all_models(fits)
  
  write.csv(results_df,paste0(output_prefix, "m2_models.csv"), row.names = FALSE)
  
  fitACE5_aug_m2 <- lapply(fits, add_prop_block_uni)
  
  colP_all_m2 <- c("va11",
                   "vc11",
                   "ve11",
                   "s1"
  )
  
  mma_var_m2 <- mxModelAverage(
    reference   = colP_all_m2,
    models      = fitACE5_aug_m2,
    include     = "all",
    SE          = TRUE,
    refAsBlock  = FALSE,
    type        = "AIC"
  )
  
  use_vp <- if (isTRUE(truth$flags$sex_differences)) truth$varparam_m else truth$varparam
  if (is.null(use_vp)) use_vp <- truth$varparam
  
  truth_raw <- c(
    va11 = use_vp$va11,
    vc11 = if (!is.null(use_vp$vc11)) use_vp$vc11 else 0,
    ve11 = use_vp$ve11
  )
  
  totV_true <- truth_raw["va11"] + truth_raw["vc11"] + truth_raw["ve11"]
  
  Pva11_true <- safe_div(truth_raw["va11"], totV_true)
  Pvc11_true <- safe_div(truth_raw["vc11"], totV_true)
  Pve11_true <- safe_div(truth_raw["ve11"], totV_true)
  
  s1_truth <- NA_real_
  if (!is.null(truth$realized) && !is.null(truth$realized$s1)) {
    s1_truth <- truth$realized$s1
  } else if (!is.null(truth$matrices$Bsib) && nrow(truth$matrices$Bsib) >= 1 && ncol(truth$matrices$Bsib) >= 2) {
    s1_truth <- truth$matrices$Bsib[1, 2]
  }
  
  truth_raw["s1"] <- s1_truth
  
  truth_std <- c(
    Pva11 = Pva11_true,
    Pvc11 = Pvc11_true,
    Pve11 = Pve11_true,
    totV  = totV_true
  )
  
  map_prop_names <- c(
    "Prop[1,1]"="Pva11",
    "Prop[1,2]"="Pvc11",
    "Prop[1,3]"="Pve11",
    "Prop[1,4]"="totV"
  )
  
  ma_tab <- mma_var_m2$`Model-Average Estimates`
  fit_params_est <- setNames(as.numeric(ma_tab[, "Estimate"]), rownames(ma_tab))
  fit_params_SE  <- setNames(as.numeric(ma_tab[, "SE"]),       rownames(ma_tab))
  
  s1_est <- if ("s1" %in% names(fit_params_est)) fit_params_est["s1"] else NA_real_
  s1_se  <- if ("s1" %in% names(fit_params_SE)) fit_params_SE["s1"] else NA_real_
  
  truth_raw["s1"] <- s1_truth

  coef_vals <- fit_params_est
  se_vec    <- fit_params_SE
  
  # Helper function to compute Prop vector from parameter vector
  compute_prop_from_params <- function(param_vals) {
    # Extract variance components
    va11_p <- if ("va11" %in% names(param_vals)) param_vals["va11"] else NA_real_
    vc11_p <- if ("vc11" %in% names(param_vals)) param_vals["vc11"] else 0
    ve11_p <- if ("ve11" %in% names(param_vals)) param_vals["ve11"] else NA_real_
    
    # Check for valid required variances
    if (!is.finite(va11_p) || !is.finite(ve11_p)) {
      return(rep(NA_real_, 5))
    }
    
    va11_p <- if (va11_p > 0) va11_p else 0
    vc11_p <- if (vc11_p > 0) vc11_p else 0
    ve11_p <- if (ve11_p > 0) ve11_p else 0
    
    # Compute total variance
    totV_p <- va11_p + vc11_p +  ve11_p
    
    # Compute proportions
    Pva11_p <- safe_div(va11_p, totV_p)
    Pvc11_p <- safe_div(vc11_p, totV_p)
    Pve11_p <- safe_div(ve11_p, totV_p)
    
    c(Pva11_p, Pvc11_p, Pve11_p, totV_p)
  }
  
  # Compute Prop vector
  prop_vec <- compute_prop_from_params(fit_params_est)
  
  # Compute delta-method SEs
  eps_base <- 1e-6
  
  # Build diagonal vcov (independence approximation)
  if (anyNA(se_vec) || any(!is.finite(se_vec))) {
    # keep NA variances as NA on the diagonal
    vcov_mat <- matrix(0, length(se_vec), length(se_vec),
                       dimnames=list(names(se_vec), names(se_vec)))
    diag(vcov_mat) <- se_vec^2
    diag(vcov_mat)[!is.finite(diag(vcov_mat))] <- NA_real_
  } else {
    vcov_mat <- diag(se_vec^2)
    rownames(vcov_mat) <- colnames(vcov_mat) <- names(se_vec)
  }
  
  if (is.null(coef_vals) || length(coef_vals) == 0 ||
      is.null(vcov_mat) || !all(is.finite(coef_vals)) ) {
    prop_ses <- rep(NA_real_, length(prop_vec))
  } else {
    
    diagV <- diag(vcov_mat)
    bad   <- is.na(diagV) | !is.finite(diagV)
    
    n_params <- length(coef_vals)
    n_props  <- length(prop_vec)
    gradients <- matrix(0, nrow = n_props, ncol = n_params)
    colnames(gradients) <- names(coef_vals)
    
    for (j in seq_len(n_params)) {
      pj <- coef_vals[j]
      h  <- eps_base * (abs(pj) + 1)
      
      coef_plus  <- coef_vals; coef_plus[j]  <- pj + h
      coef_minus <- coef_vals; coef_minus[j] <- pj - h
      
      prop_plus  <- tryCatch(compute_prop_from_params(coef_plus),  error=function(e) rep(NA_real_, n_props))
      prop_minus <- tryCatch(compute_prop_from_params(coef_minus), error=function(e) rep(NA_real_, n_props))
      
      ok <- is.finite(prop_plus) & is.finite(prop_minus)
      gradients[ok, j] <- (prop_plus[ok] - prop_minus[ok]) / (2 * h)
    }
    
    if (any(bad)) gradients[, bad] <- 0
    
    prop_ses <- apply(gradients, 1, function(g) {
      se_sq <- sum((g^2) * diag(vcov_mat), na.rm = TRUE)  # diagonal-only
      if (is.finite(se_sq) && se_sq >= 0) sqrt(se_sq) else NA_real_
    })
  }
  
  # Build estimates table
  prop_rownames <- names(map_prop_names)
  estimates_df <- data.frame(
    param    = prop_rownames,
    estimate = prop_vec,
    SE       = prop_ses,
    source   = "standardized",
    stringsAsFactors = FALSE
  )
  estimates_df$lcl95 <- estimates_df$estimate - 1.96 * estimates_df$SE
  estimates_df$ucl95 <- estimates_df$estimate + 1.96 * estimates_df$SE
  
  # Map param names and add scale/truth
  estimates_df$param <- unname(map_prop_names[estimates_df$param])
  
  scale_tag <- function(p) {
    if (p %in% names(truth_raw)) "raw" else "standardized"
  }
  
  estimates_df$scale <- vapply(estimates_df$param, scale_tag, FUN.VALUE = character(1))
  estimates_df$truth <- mapply(function(p, sc) {
    if (sc == "raw") {
      if (p %in% names(truth_raw)) truth_raw[[p]] else NA_real_
    } else {
      if (p %in% names(truth_std)) truth_std[[p]] else NA_real_
    }
  }, estimates_df$param, estimates_df$scale)
  
  sib_df <- data.frame(
    param    = "s1",
    estimate = s1_est,
    SE       = s1_se,
    source   = "raw",
    lcl95    = s1_est - 1.96 * s1_se,
    ucl95    = s1_est + 1.96 * s1_se,
    scale    = "raw",
    truth    = s1_truth,
    stringsAsFactors = FALSE
  )
  
  estimates_df <- rbind(estimates_df, sib_df)
  rownames(estimates_df) <- NULL

  # Compute performance metrics
  performance_metrics <- within(estimates_df, {
    bias      <- estimate - truth
    MAE       <- abs(bias)
    MSE       <- bias^2
    covered95 <- ifelse(!is.na(truth) & (truth >= lcl95) & (truth <= ucl95), 1L, 0L)
    typeI     <- ifelse(!is.na(truth) & abs(truth) < 1e-12 & (lcl95 > 0 | ucl95 < 0), 1L, 0L)
    typeII    <- ifelse(!is.na(truth) & abs(truth) >= 1e-12 & (lcl95 <= 0 & ucl95 >= 0), 1L, 0L)
  })
  
  # Write CSV outputs
  utils::write.csv(estimates_df,
                   paste0(output_prefix, "model_estimates_m2.csv"),
                   row.names = FALSE)
  utils::write.csv(performance_metrics,
                   paste0(output_prefix, "performance_metrics_m2.csv"),
                   row.names = FALSE)
  
  # Return results
  res_m2 <- list(estimates = estimates_df, metrics = performance_metrics)
  
  
  
  model_fits_m2_df <- data.frame(
    model      = mma_var_m2$'Akaike-Weights Table'$'model',
    AIC        = mma_var_m2$'Akaike-Weights Table'$'AIC',
    deltaAIC   = mma_var_m2$'Akaike-Weights Table'$'delta',
    weight     = mma_var_m2$'Akaike-Weights Table'$'AkaikeWeight',
    in_confset = mma_var_m2$'Akaike-Weights Table'$'inConfidenceSet' == "*",
    stringsAsFactors = FALSE
  )
  utils::write.csv(model_fits_m2_df,
                   paste0(output_prefix, "model_fits_m2.csv"),
                   row.names = FALSE)
  
  cat("Simulated data saved to:", paste0(output_prefix, "sim_data.csv"), "\n")
  cat("Model fits (Method 1) saved to:",
      paste0(output_prefix, "model_fits_method1_global_all.csv"), "\n")
  cat("Model fits (Method 2) saved to:",
      paste0(output_prefix, "model_fits_method2_global_all.csv"), "\n")
  cat("Model fits (Method 3) saved to:",
      paste0(output_prefix, "model_fits_method3_global_all.csv"), "\n")
  cat("Multimodel averaging estimates (Method 1) saved to:",
      paste0(output_prefix, "model_averaged_estimates_method1_global_all.csv"), "\n")
  cat("Multimodel averaging estimates (Method 2) saved to:",
      paste0(output_prefix, "model_averaged_estimates_method2_global_all.csv"), "\n")
  cat("Multimodel averaging estimates (Method 3) saved to:",
      paste0(output_prefix, "model_averaged_estimates_method3_global_all.csv"), "\n")
  cat("Performance metrics (Method 1) saved to:",
      paste0(output_prefix, "performance_metrics_method1_global_all.csv"), "\n")
  cat("Performance metrics (Method 2) saved to:",
      paste0(output_prefix, "performance_metrics_method2_global_all.csv"), "\n\n")
  cat("Performance metrics (Method 3) saved to:",
      paste0(output_prefix, "performance_metrics_method3_global_all.csv"), "\n\n")
  
  invisible(list(
    mma_var_m2        = mma_var_m2,
    fitACE5_aug_m2    = fitACE5_aug_m2,
    method1_estimates = res_m0$estimates,
    method0_metrics   = res_m0$metrics,
    method1_estimates = res_m1$estimates,
    method1_metrics   = res_m1$metrics,
    method2_estimates = res_m2$estimates,
    method2_metrics   = res_m2$metrics,
    sim_data          = sim_data,
    truth             = truth,
    model_fits_m2     = model_fits_m2_df))
}




# ------------------------------------------------------------------------------
# 10) Run the simulation and analysis
# ------------------------------------------------------------------------------

n_reps <- 1000


seed_i <- 0L

for (i in seq_len(n_reps)) {
  attempt <- 0L
  repeat {
    seed_i  <- seed_i + 1L
    attempt <- attempt + 1L
    
    message(sprintf("Running replication %d of %d (attempt %d) with seed %d",
                    i, n_reps, attempt, seed_i))
    
    res <- tryCatch(
      {
        run_simulation_and_analysis(seed = seed_i)
      },
      error = function(e) {
        message(sprintf("  -> attempt %d failed: %s", attempt, conditionMessage(e)))
        NULL
      }
    )
    
    # Break the repeat loop on success (accepted dataset)
    if (!is.null(res)) break
  }
  
  prefix <- sprintf("simulation_results_(A=0.70, C=.20, E=0.10, S=0.30)/run_%04d", i)
  
  if (!is.null(res$data)) {
    utils::write.csv(res$data, paste0(prefix, "_data.csv"), row.names = FALSE)
  }
  res$performance_metrics$seed <- seed_i
  utils::write.csv(res$model_fits_m2,            paste0(prefix, "_model_fits_m2.csv"),        row.names = FALSE)
  utils::write.csv(res$averaged_estimates,    paste0(prefix, "_model_averaged.csv"),    row.names = FALSE)
  utils::write.csv(res$performance_metrics,   paste0(prefix, "_performance.csv"),       row.names = FALSE)
}
