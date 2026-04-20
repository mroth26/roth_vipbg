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
require(mvtnorm)  
options(width = 245)
mxVersion()
mxOption(NULL, "Default optimizer", "NPSOL")

# ------------------------------------------------------------------------------
# Data simulation function for bivariate models
# ------------------------------------------------------------------------------
source('models/bivariate_sim.R')

# ------------------------------------------------------------------------------
# Single Model Selection Algorithm
# ------------------------------------------------------------------------------

single_model_select <- function(fit1, fit2, fit3, fits, fitFace, alpha = 0.05, silent = TRUE) {
  # --- helpers (kept inside function for portability) ---
  get_aic <- function(m) as.data.frame(summary(m)$informationCriteria)$par[1]
  safe_mxCompare <- function(m1, m2) tryCatch(as.data.frame(mxCompare(m1, m2)), error = function(e) NULL)
  
  # --- 0) saturated model ---
  sat0 <- mxRefModels(fitFace)$Saturated
  sat  <- mxRun(sat0, silent = TRUE)
  
  # --- 1) absolute fit screening for each tier ---
  cmp3 <- safe_mxCompare(sat, fit3)
  if (is.null(cmp3) || nrow(cmp3) < 2) stop("mxCompare(sat, fit3) failed")
  abs3_tab <- data.frame(
    tier  = "fit3",
    model = names(fit3),
    AIC   = as.numeric(sapply(fit3, get_aic)[names(fit3)]),
    p     = as.numeric(cmp3[-1, "p"]),
    df    = as.numeric(cmp3[-1, "df"]),
    diffLL= as.numeric(cmp3[-1, "diffLL"]),
    pass  = is.finite(as.numeric(cmp3[-1, "p"])) & (as.numeric(cmp3[-1, "p"]) >= alpha),
    stringsAsFactors = FALSE
  )
  passed3 <- abs3_tab$model[abs3_tab$pass]
  
  if (length(passed3) > 0) {
    tier_used <- "fit3"
    full_candidates <- passed3
    abs_used <- abs3_tab
  } else {
    cmp2 <- safe_mxCompare(sat, fit2)
    if (is.null(cmp2) || nrow(cmp2) < 2) stop("mxCompare(sat, fit2) failed")
    abs2_tab <- data.frame(
      tier  = "fit2",
      model = names(fit2),
      AIC   = as.numeric(sapply(fit2, get_aic)[names(fit2)]),
      p     = as.numeric(cmp2[-1, "p"]),
      df    = as.numeric(cmp2[-1, "df"]),
      diffLL= as.numeric(cmp2[-1, "diffLL"]),
      pass  = is.finite(as.numeric(cmp2[-1, "p"])) & (as.numeric(cmp2[-1, "p"]) >= alpha),
      stringsAsFactors = FALSE
    )
    passed2 <- abs2_tab$model[abs2_tab$pass]
    
    if (length(passed2) > 0) {
      tier_used <- "fit2"
      full_candidates <- passed2
      abs_used <- abs2_tab
    } else {
      cmp1 <- safe_mxCompare(sat, fit1)
      if (is.null(cmp1) || nrow(cmp1) < 2) stop("mxCompare(sat, fit1) failed")
      abs1_tab <- data.frame(
        tier  = "fit1",
        model = names(fit1),
        AIC   = as.numeric(sapply(fit1, get_aic)[names(fit1)]),
        p     = as.numeric(cmp1[-1, "p"]),
        df    = as.numeric(cmp1[-1, "df"]),
        diffLL= as.numeric(cmp1[-1, "diffLL"]),
        pass  = is.finite(as.numeric(cmp1[-1, "p"])) & (as.numeric(cmp1[-1, "p"]) >= alpha),
        stringsAsFactors = FALSE
      )
      passed1 <- abs1_tab$model[abs1_tab$pass]
      
      if (length(passed1) > 0) {
        tier_used <- "fit1_pass"
        aics_pass <- setNames(abs1_tab$AIC, abs1_tab$model)[passed1]
        selected_name <- names(which.min(aics_pass))
        single_fit <- fits[[selected_name]]
        if (!silent) {
          cat("Tier used:", tier_used, "\n")
          cat("Selected single model (fit1 among pass):", selected_name, "\n")
        }
        return(list(
          single_fit = single_fit,
          selected_name = selected_name,
          tier_used = tier_used,
          full_candidates = passed1,
          abs_table = abs1_tab,
          nested_trace = NULL,
          final_aics = setNames(abs1_tab$AIC, abs1_tab$model)[passed1]
        ))
      } else {
        tier_used <- "fit1_fallback_lowestAIC"
        aics_all <- setNames(abs1_tab$AIC, abs1_tab$model)
        selected_name <- names(which.min(aics_all))
        single_fit <- fits[[selected_name]]
        if (!silent) {
          cat("Tier used:", tier_used, "\n")
          cat("Selected single model (fit1 fallback):", selected_name, "\n")
        }
        return(list(
          single_fit = single_fit,
          selected_name = selected_name,
          tier_used = tier_used,
          full_candidates = names(fit1),
          abs_table = abs1_tab,
          nested_trace = NULL,
          final_aics = aics_all
        ))
      }
    }
  }
  
  if (!silent) {
    cat("Tier used:", tier_used, "\n")
    cat("Full models retained after absolute fit:", paste(full_candidates, collapse = ", "), "\n")
  }
  
  # --- 2) nested reduction (only for fit3/fit2 tiers) ---
  nesting_map <- list(
    # correlated factor nesting
    Face = c("Fac", "Fae", "Fce", "Fa", "Fc", "Fe"),
    Fac  = c("Fa", "Fc"),
    Fae  = c("Fa", "Fe"),
    Fce  = c("Fc", "Fe"),
    
    # directional nesting
    Dor  = c("Do", "Dr"),
    
    # bidirectional hybrid nesting
    Haor = c("Fa","Hao","Har","Do","Dr"),
    Hcor = c("Fc","Hco","Hcr","Do","Dr"),
    Heor = c("Fe","Heo","Her","Do","Dr"),
    
    # single-direction hybrid nesting
    Haco = c("Fa","Fc","Fac","Hao","Hco","Do"),
    Hacr = c("Fa","Fc","Fac","Har","Hcr","Dr"),
    Haeo = c("Fa","Fe","Fae","Hao","Heo","Do"),
    Haer = c("Fa","Fe","Fae","Har","Her","Dr"),
    Hceo = c("Fc","Fe","Fce","Hco","Heo","Do"),
    Hcer = c("Fc","Fe","Fce","Hcr","Her","Dr")
  )
  
  nested_trace <- data.frame(
    parent = character(),
    child  = character(),
    p = numeric(),
    df = numeric(),
    diffLL = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  
  candidate_pool <- list()
  
  for (parent_nm in full_candidates) {
    parent_model <- fits[[parent_nm]]
    child_names <- nesting_map[[parent_nm]]
    
    if (is.null(child_names) || length(child_names) == 0) {
      candidate_pool[[parent_nm]] <- parent_model
      next
    }
    
    child_names <- child_names[child_names %in% names(fits)]
    if (length(child_names) == 0) {
      candidate_pool[[parent_nm]] <- parent_model
      next
    }
    
    retained_children <- character(0)
    
    for (child_nm in child_names) {
      child_model <- fits[[child_nm]]
      cmp <- safe_mxCompare(parent_model, child_model)
      
      if (is.null(cmp) || nrow(cmp) < 2) {
        nested_trace <- rbind(nested_trace, data.frame(
          parent = parent_nm, child = child_nm,
          p = NA_real_, df = NA_real_, diffLL = NA_real_,
          decision = "compare_error",
          stringsAsFactors = FALSE
        ))
        next
      }
      
      pval <- as.numeric(cmp[2, "p"])
      dfv  <- as.numeric(cmp[2, "df"])
      dll  <- as.numeric(cmp[2, "diffLL"])
      keep_child <- is.finite(pval) && pval >= alpha
      
      nested_trace <- rbind(nested_trace, data.frame(
        parent = parent_nm, child = child_nm,
        p = pval, df = dfv, diffLL = dll,
        decision = if (keep_child) "retain_child" else "reject_child",
        stringsAsFactors = FALSE
      ))
      
      if (keep_child) retained_children <- c(retained_children, child_nm)
    }
    
    if (length(retained_children) > 0) {
      for (cn in retained_children) candidate_pool[[cn]] <- fits[[cn]]
    } else {
      candidate_pool[[parent_nm]] <- parent_model
    }
  }
  
  # unique by name
  candidate_pool <- candidate_pool[!duplicated(names(candidate_pool))]
  
  if (!silent) {
    cat("Candidates after nested reduction:", paste(names(candidate_pool), collapse = ", "), "\n")
  }
  
  # --- 3) final selection by AIC ---
  cand_aics <- sapply(candidate_pool, get_aic)
  selected_name <- names(which.min(cand_aics))
  single_fit <- candidate_pool[[selected_name]]
  
  if (!silent) cat("Selected single model:", selected_name, "\n")
  
  list(
    single_fit = single_fit,
    selected_name = selected_name,
    tier_used = tier_used,
    full_candidates = full_candidates,
    abs_table = abs_used,
    nested_trace = nested_trace,
    final_aics = cand_aics
  )
}



# ------------------------------------------------------------------------------
# Bivariate Methods Function
# ------------------------------------------------------------------------------
fit_bivariate <- function(
    data,
    topName = "bi",
    selVars = c("p1_t1","p1_t2","p2_t1","p2_t2"),
    nv = 2,
    lboundVar = NA,
    
    means = list(me1 = c(TRUE, 0), me2 = c(TRUE, 0)),
    
    A_diag = list(a11 = c(TRUE, 0.3), a22 = c(TRUE, 0.3)),
    C_diag = list(c11 = c(TRUE, 0.3), c22 = c(TRUE, 0.3)),
    E_diag = list(e11 = c(TRUE, 0.3), e22 = c(TRUE, 0.3)),
    
    R = list(
      rA12 = c(FALSE, 0),
      rC12 = c(FALSE, 0),
      rE12 = c(FALSE, 0)
    ),
    
    B = list(
      b12 = c(FALSE, 0),
      b21 = c(FALSE, 0)
    ),
    
    S = list(
      s1 = c(FALSE, 0),
      s2 = c(FALSE, 0)
    ),
    
    ci = character(0),
    tryHard = TRUE,
    extraTries = 15,
    intervals = TRUE
) {
  stopifnot(nv == 2)
  
  dataMZ <- subset(data, zyg %in% c(1))
  dataDZ <- subset(data, zyg %in% c(2))
  
  get_fs <- function(x, nm) {
    if (is.null(x) || length(x) != 2) stop("Parameter '", nm, "' must be c(free, start).")
    list(free = as.logical(x[[1]]), start = as.numeric(x[[2]]))
  }
  
  me1 <- get_fs(means$me1, "me1")
  me2 <- get_fs(means$me2, "me2")
  
  M <- mxMatrix("Full", 1, nv,
                free   = c(me1$free, me2$free),
                values = c(me1$start, me2$start),
                labels = c("me1","me2"),
                name   = "M")
  
  expMean <- mxAlgebra(
    cbind(M[1,1], M[1,1],
          M[1,2], M[1,2]),
    name="expMean"
  )
  
  for (nm in c("b12","b21")) if (is.null(B[[nm]])) B[[nm]] <- c(FALSE, 0)
  b12 <- get_fs(B$b12, "b12")
  b21 <- get_fs(B$b21, "b21")
  
  Bmat_free   <- matrix(FALSE, 2, 2)
  Bmat_values <- matrix(0,     2, 2)
  Bmat_labels <- matrix(NA_character_, 2, 2)
  
  setB <- function(i, j, spec, nm) {
    Bmat_free[i,j]   <<- spec$free
    Bmat_values[i,j] <<- spec$start
    if (spec$free) Bmat_labels[i,j] <<- nm
  }
  
  setB(1,2, b12, "b12")
  setB(2,1, b21, "b21")
  
  Bmx <- mxMatrix("Full", 2, 2,
                  free   = as.vector(t(Bmat_free)),
                  values = as.vector(t(Bmat_values)),
                  labels = as.vector(t(Bmat_labels)),
                  lbound = as.vector(t(ifelse(Bmat_free, -0.95, NA))),
                  ubound = as.vector(t(ifelse(Bmat_free,  0.95, NA))),
                  byrow  = TRUE,
                  name   = "B")
  
  I  <- mxMatrix("Iden", nv, nv, name="I")
  invIB <- mxAlgebra(solve(I - B), name="invIB")
  
  for (nm in c("s1","s2")) if (is.null(S[[nm]])) S[[nm]] <- c(FALSE, 0)
  s1 <- get_fs(S$s1, "s1")
  s2 <- get_fs(S$s2, "s2")
  
  Bsib_free   <- matrix(FALSE, 4, 4)
  Bsib_values <- matrix(0,     4, 4)
  Bsib_labels <- matrix(NA_character_, 4, 4)
  
  setSib <- function(i, j, spec, nm) {
    Bsib_free[i,j]   <<- spec$free
    Bsib_values[i,j] <<- spec$start
    if (spec$free) Bsib_labels[i,j] <<- nm
  }
  
  setSib(1, 2, s1, "s1")
  setSib(2, 1, s1, "s1")
  setSib(3, 4, s2, "s2")
  setSib(4, 3, s2, "s2")
  
  Bsib <- mxMatrix("Full", 4, 4,
                   free   = as.vector(t(Bsib_free)),
                   values = as.vector(t(Bsib_values)),
                   labels = as.vector(t(Bsib_labels)),
                   lbound = as.vector(t(ifelse(Bsib_free, -0.95, NA))),
                   ubound = as.vector(t(ifelse(Bsib_free,  0.95, NA))),
                   byrow  = TRUE,
                   name   = "Bsib")
  
  I4 <- mxMatrix("Iden", 4, 4, name="I4")
  G_sib <- mxAlgebra(solve(I4 - Bsib), name="G_sib")
  
  a11 <- get_fs(A_diag$a11, "A_diag$a11")
  a22 <- get_fs(A_diag$a22, "A_diag$a22")
  c11 <- get_fs(C_diag$c11, "C_diag$c11")
  c22 <- get_fs(C_diag$c22, "C_diag$c22")
  e11 <- get_fs(E_diag$e11, "E_diag$e11")
  e22 <- get_fs(E_diag$e22, "E_diag$e22")
  
  mkDiagVar <- function(specs, labels, name) {
    free   <- c(specs[[1]]$free, specs[[2]]$free)
    vals_v <- c(specs[[1]]$start, specs[[2]]$start)
    if (any(!is.finite(vals_v)) || any(vals_v < 0, na.rm=TRUE)) stop(name, " start variances must be finite and >= 0.")
    lbound <- ifelse(free, lboundVar, NA)
    labs   <- ifelse(free, labels, NA_character_)
    mxMatrix("Diag", 2, 2,
             free   = free,
             values = vals_v,
             lbound = lbound,
             labels = labs,
             name   = name)
  }
  
  aV <- mkDiagVar(list(a11,a22), c("va11","va22"), "aV")
  cV <- mkDiagVar(list(c11,c22), c("vc11","vc22"), "cV")
  eV <- mkDiagVar(list(e11,e22), c("ve11","ve22"), "eV")
  
  aSD <- mxAlgebra(vec2diag(sqrt(diag2vec(aV))), name="aSD")
  cSD <- mxAlgebra(vec2diag(sqrt(diag2vec(cV))), name="cSD")
  eSD <- mxAlgebra(vec2diag(sqrt(diag2vec(eV))), name="eSD")
  
  for (nm in c("rA12","rC12","rE12")) if (is.null(R[[nm]])) R[[nm]] <- c(FALSE, 0)
  rA12 <- get_fs(R$rA12, "rA12")
  rC12 <- get_fs(R$rC12, "rC12")
  rE12 <- get_fs(R$rE12, "rE12")
  
  mkStand2 <- function(r12, prefix, name) {
    freeM <- matrix(FALSE, 2, 2)
    valM  <- diag(1,2)
    labM  <- matrix(NA_character_, 2, 2)
    
    freeM[1,2] <- r12$free
    freeM[2,1] <- r12$free
    valM[1,2]  <- r12$start
    valM[2,1]  <- r12$start
    if (r12$free) {
      labM[1,2] <- paste0(prefix,"12")
      labM[2,1] <- paste0(prefix,"12")
    }
    
    mxMatrix("Stand", 2, 2,
             free   = as.vector(t(freeM)),
             values = as.vector(t(valM)),
             labels = as.vector(t(labM)),
             lbound = as.vector(t(ifelse(freeM, -1, NA))),
             ubound = as.vector(t(ifelse(freeM,  1, NA))),
             byrow  = TRUE,
             name   = name)
  }
  
  RA <- mkStand2(rA12, "rA", "RA")
  RC <- mkStand2(rC12, "rC", "RC")
  RE <- mkStand2(rE12, "rE", "RE")
  
  A <- mxAlgebra(aSD %*% RA %*% aSD, name="A")
  C <- mxAlgebra(cSD %*% RC %*% cSD, name="C")
  E <- mxAlgebra(eSD %*% RE %*% eSD, name="E")
  
  V   <- mxAlgebra(invIB %*% (A + C + E)   %*% t(invIB), name="V")
  cMZ <- mxAlgebra(invIB %*% (A + C)       %*% t(invIB), name="cMZ")
  cDZ <- mxAlgebra(invIB %*% (0.5 * A + C) %*% t(invIB), name="cDZ")
  
  expCovMZ <- mxAlgebra(rbind(
    cbind( V[1,1],   cMZ[1,1], V[1,2],   cMZ[1,2] ),
    cbind( cMZ[1,1], V[1,1],   cMZ[1,2], V[1,2]   ),
    cbind( V[2,1],   cMZ[2,1], V[2,2],   cMZ[2,2] ),
    cbind( cMZ[2,1], V[2,1],   cMZ[2,2], V[2,2]   )
  ), name="expCovMZ")
  
  expCovDZ <- mxAlgebra(rbind(
    cbind( V[1,1],   cDZ[1,1], V[1,2],   cDZ[1,2] ),
    cbind( cDZ[1,1], V[1,1],   cDZ[1,2], V[1,2]   ),
    cbind( V[2,1],   cDZ[2,1], V[2,2],   cDZ[2,2] ),
    cbind( cDZ[2,1], V[2,1],   cDZ[2,2], V[2,2]   )
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
                 Bmx, I, invIB,
                 Bsib, I4, G_sib,
                 aV, cV, eV,
                 aSD, cSD, eSD,
                 RA, RC, RE,
                 A, C, E,
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


# ------------------------------------------------------------------------------
# Attach standardized-proportion algebra for mxModelAverage
# ------------------------------------------------------------------------------

add_prop_block_biv <- function(fit) {
  add_if_missing <- function(model, name, algebra_or_matrix) {
    if (!(name %in% names(model$algebras)) && !(name %in% names(model$matrices))) {
      model <- mxModel(model, algebra_or_matrix)
    }
    model
  }
  
  # Require A, C, E at the top level
  req <- c("A","C","E")
  missing_req <- setdiff(req, c(names(fit$algebras), names(fit$matrices)))
  if (length(missing_req)) stop("downstream: missing required objects: ", paste(missing_req, collapse=", "))
  
  # Ensure B exists; if absent, add 2x2 Zero
  if (!("B" %in% names(fit$algebras) || "B" %in% names(fit$matrices))) {
    fit <- mxModel(fit, mxMatrix(type="Zero", nrow=2, ncol=2, name="B"))
  }
  
  # Prefer an existing V (phenotypic variance). If none, build DoC-aware V if invIB exists; else fallback to A+C+E
  hasV    <- "V"    %in% c(names(fit$algebras), names(fit$matrices))
  hasinvB <- "invIB"%in% names(fit$algebras)
  if (!hasV) {
    if (hasinvB) {
      fit <- mxModel(fit, mxAlgebra(invIB %*% (A + C + E) %*% t(invIB), name="V"))
    } else {
      fit <- mxModel(fit, mxAlgebra(A + C + E, name="V"))
    }
  }
  
  # Inverse SD matrix from phenotypic V
  fit <- add_if_missing(fit, "iSD", mxAlgebra(vec2diag(1 / sqrt(diag2vec(V))), name="iSD"))
  
  # Proportions, cross-trait correlations, standardized DoC
  if (!("Prop" %in% names(fit$algebras))) {
    estProp <- mxAlgebra(
      cbind(
        A[1,1]/V[1,1], C[1,1]/V[1,1], E[1,1]/V[1,1],
        A[2,2]/V[2,2], C[2,2]/V[2,2], E[2,2]/V[2,2],
        A[2,1]/sqrt(A[1,1]*A[2,2]),     # rA21
        C[2,1]/sqrt(C[1,1]*C[2,2]),     # rC21
        E[2,1]/sqrt(E[1,1]*E[2,2]),     # rE21
        iSD[2,2]*B[2,1]/iSD[1,1],       # sb21
        iSD[1,1]*B[1,2]/iSD[2,2],       # sb12
        V[2,1]/sqrt(V[1,1]*V[2,2]),     # rP21
        V[1,1] + V[2,2]
      ),
      name="Prop",
      dimnames=list("pro", c("Pva11","Pvc11","Pve11",
                             "Pva22","Pvc22","Pve22",
                             "rA21","rC21","rE21",
                             "sb21","sb12","rP21","totV"))
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
  
  if (identical(output_prefix, "")) {
    output_prefix <- sprintf("simulation_results_Face/run_%04d_", as.integer(seed))
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
    seed  = seed,
    
    rA12 = 0.30 + runif(1, -0.05, 0.05),
    rC12 = 0.20 + runif(1, -0.05, 0.05),
    rE12 = 0.10 + runif(1, -0.05, 0.05),
    
    b12 = 0 + runif(1, -0.01, 0.01),
    b21 = 0 + runif(1, -0.01, 0.01),
    
    s = c(0, 0) + runif(2, -0.01, 0.01),
    
    
    A_diag = c(0.80, 0.33) + runif(2, -0.01, 0.01),
    C_diag = c(0.10, 0.33) + runif(2, -0.01, 0.01),
    E_diag = c(0.10, 0.34) + runif(2, -0.01, 0.01),
    
    
    mean = c(0, 0),
  )
  
  sim_data <- sim_results$data
  truth    <- sim_results$truth
  
  utils::write.csv(sim_data, paste0(output_prefix, "sim_data.csv"), row.names = FALSE)
  # ========================================================================
  # Precomputation of Start Values
  # ========================================================================
  
  vars4 <- c("p1_t1","p1_t2","p2_t1","p2_t2")
  
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
  
  cross_cross <- cov_all["p1_t1","p2_t2"]
  
  cross_within <- cov_all["p1_t1","p2_t1"]
  
  sv_rA12 <- clip(cross_cross, -0.95, 0.95)
  sv_rC12 <- clip(cross_cross, -0.95, 0.95)
  sv_rE12 <- clip(cross_within - cross_cross, -0.95, 0.95)
  
  sv_b21  <- clip(cross_cross, -0.50, 0.50)
  sv_b12  <- clip(cross_cross, -0.50, 0.50)
  
  rMZ_p1 <- cov_mz["p1_t1","p1_t2"]
  rDZ_p1 <- cov_dz["p1_t1","p1_t2"]
  rMZ_p2 <- cov_mz["p2_t1","p2_t2"]
  rDZ_p2 <- cov_dz["p2_t1","p2_t2"]
  
  h2_1 <- clip(2 * (rMZ_p1 - rDZ_p1), 0.05, 0.95)
  c2_1 <- clip(rMZ_p1 - h2_1,         0.00, 0.90)
  e2_1 <- clip(1 - h2_1 - c2_1,       0.05, 0.95)
  
  h2_2 <- clip(2 * (rMZ_p2 - rDZ_p2), 0.05, 0.95)
  c2_2 <- clip(rMZ_p2 - h2_2,         0.00, 0.90)
  e2_2 <- clip(1 - h2_2 - c2_2,       0.05, 0.95)
  
  sv_va11 <- h2_1; sv_vc11 <- c2_1; sv_ve11 <- e2_1
  sv_va22 <- h2_2; sv_vc22 <- c2_2; sv_ve22 <- e2_2
  
  starts <- c(
    sv_rA12 = sv_rA12, sv_rC12 = sv_rC12, sv_rE12 = sv_rE12,
    sv_b21  = sv_b21,  sv_b12  = sv_b12,
    sv_va11 = sv_va11, sv_vc11 = sv_vc11, sv_ve11 = sv_ve11,
    sv_va22 = sv_va22, sv_vc22 = sv_vc22, sv_ve22 = sv_ve22
  )
  
  
  # ========================================================================
  # METHOD 0 – Quality Control if DGP were known
  # ========================================================================
  
  fitFa   <- fit_bivariate( data = sim_data, topName = "Fa", R = list(rA12 = c(TRUE, sv_rA12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitFc   <- fit_bivariate( data = sim_data, topName = "Fc", R = list(rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitFe   <- fit_bivariate( data = sim_data, topName = "Fe", R = list(rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitDo   <- fit_bivariate( data = sim_data, topName = "Do", B = list(b21 = c(TRUE, sv_b21)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitDr   <- fit_bivariate( data = sim_data, topName = "Dr", B = list(b12 = c(TRUE, sv_b12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitFac  <- fit_bivariate( data = sim_data, topName = "Fac", R = list(rA12 = c(TRUE, sv_rA12), rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitFae  <- fit_bivariate( data = sim_data, topName = "Fae", R = list(rA12 = c(TRUE, sv_rA12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitFce  <- fit_bivariate( data = sim_data, topName = "Fce", R = list(rC12 = c(TRUE, sv_rC12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHao  <- fit_bivariate( data = sim_data, topName = "Hao", B = list(b21 = c(TRUE, sv_b21)),R = list(rA12 = c(TRUE, sv_rA12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )  
  fitHar  <- fit_bivariate( data = sim_data, topName = "Har", B = list(b12 = c(TRUE, sv_b12)),R = list(rA12 = c(TRUE, sv_rA12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHco  <- fit_bivariate( data = sim_data, topName = "Hco", B = list(b21 = c(TRUE, sv_b21)),R = list(rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHcr  <- fit_bivariate( data = sim_data, topName = "Hcr", B = list(b12 = c(TRUE, sv_b12)),R = list(rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHeo  <- fit_bivariate( data = sim_data, topName = "Heo", B = list(b21 = c(TRUE, sv_b21)),R = list(rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHer  <- fit_bivariate( data = sim_data, topName = "Her", B = list(b12 = c(TRUE, sv_b12)),R = list(rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitDor  <- fit_bivariate( data = sim_data, topName = "Dor", B = list(b21 = c(TRUE, sv_b21), b12 = c(TRUE, sv_b12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitFace <- fit_bivariate( data = sim_data, topName = "Face", R = list(rA12 = c(TRUE, sv_rA12), rC12 = c(TRUE, sv_rC12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHaco <- fit_bivariate( data = sim_data, topName = "Haco", B = list(b21 = c(TRUE, sv_b21)),R = list(rA12 = c(TRUE, sv_rA12), rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHacr <- fit_bivariate( data = sim_data, topName = "Hacr", B = list(b12 = c(TRUE, sv_b12)),R = list(rA12 = c(TRUE, sv_rA12), rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHaeo <- fit_bivariate( data = sim_data, topName = "Haeo", B = list(b21 = c(TRUE, sv_b21)),R = list(rA12 = c(TRUE, sv_rA12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHaer <- fit_bivariate( data = sim_data, topName = "Haer", B = list(b12 = c(TRUE, sv_b12)),R = list(rA12 = c(TRUE, sv_rA12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHceo <- fit_bivariate( data = sim_data, topName = "Hceo", B = list(b21 = c(TRUE, sv_b21)),R = list(rC12 = c(TRUE, sv_rC12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHcer <- fit_bivariate( data = sim_data, topName = "Hcer", B = list(b12 = c(TRUE, sv_b12)),R = list(rC12 = c(TRUE, sv_rC12), rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHaor <- fit_bivariate( data = sim_data, topName = "Haor", B = list(b21 = c(TRUE, sv_b21), b12 = c(TRUE, sv_b12)),R = list(rA12 = c(TRUE, sv_rA12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHcor <- fit_bivariate( data = sim_data, topName = "Hcor", B = list(b21 = c(TRUE, sv_b21), b12 = c(TRUE, sv_b12)),R = list(rC12 = c(TRUE, sv_rC12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  fitHeor <- fit_bivariate( data = sim_data, topName = "Heor", B = list(b21 = c(TRUE, sv_b21), b12 = c(TRUE, sv_b12)),R = list(rE12 = c(TRUE, sv_rE12)),
                            A_diag = list(a11 = c(TRUE, sv_va11), a22 = c(TRUE, sv_va22)),
                            C_diag = list(c11 = c(TRUE, sv_vc11), c22 = c(TRUE, sv_vc22)),
                            E_diag = list(e11 = c(TRUE, sv_ve11), e22 = c(TRUE, sv_ve22)),
  )
  
  fits <- list(
    Face = fitFace
  )
  
  single_fit <- fits[[1]]
  results_df <- process_all_models(fits)
  
  write.csv(results_df,paste0(output_prefix, "models.csv"), row.names = FALSE)
  
  use_vp <- if (isTRUE(truth$flags$sex_differences)) truth$varparam_m else truth$varparam
  if (is.null(use_vp)) use_vp <- truth$varparam
  
  truth_raw <- c(
    va11 = use_vp$va11, va22 = use_vp$va22, va12 = use_vp$va12,
    vc11 = use_vp$vc11, vc22 = use_vp$vc22, vc12 = use_vp$vc12,
    ve11 = use_vp$ve11, ve22 = use_vp$ve22, ve12 = use_vp$ve12,
    sb12  = use_vp$b12,  sb21  = use_vp$b21
  )
  
  Pva11_true <- safe_div(use_vp$va11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pvc11_true <- safe_div(use_vp$vc11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pve11_true <- safe_div(use_vp$ve11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pva22_true <- safe_div(use_vp$va22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  Pvc22_true <- safe_div(use_vp$vc22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  Pve22_true <- safe_div(use_vp$ve22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  
  rA_true <- safe_div(use_vp$va12, sqrt(use_vp$va11 * use_vp$va22))
  rC_true <- safe_div(use_vp$vc12, sqrt(use_vp$vc11 * use_vp$vc22))
  rE_true <- safe_div(use_vp$ve12, sqrt(use_vp$ve11 * use_vp$ve22))
  
  A_true <- mat2(use_vp$va11, use_vp$va22, use_vp$va12)
  C_true <- mat2(use_vp$vc11, use_vp$vc22, use_vp$vc12)
  E_true <- mat2(use_vp$ve11, use_vp$ve22, use_vp$ve12)
  B_true <- matrix(c(0, use_vp$b12, use_vp$b21, 0), nrow=2, byrow=TRUE)
  
  build_truth_V <- function(A, C, E, B) {
    I2 <- diag(2)
    IB <- I2 - B
    if (any(!is.finite(IB)) || abs(det(IB)) < .Machine$double.eps) return(NULL)
    invIB <- solve(IB)
    invIB %*% (A + C + E) %*% t(invIB)
  }
  
  std_betas_from_V <- function(B, V) {
    if (is.null(V)) return(c(sb12 = NA_real_, sb21 = NA_real_))
    iSD <- diag(1 / sqrt(diag(V)))
    c(sb21 = iSD[2,2] * B[2,1] / iSD[1,1],
      sb12 = iSD[1,1] * B[1,2] / iSD[2,2])
  }
  
  V_truth <- NULL
  if (isTRUE(truth$flags$sex_differences)) {
    if (!is.null(truth$matrices$V_m)) V_truth <- truth$matrices$V_m
  } else {
    if (!is.null(truth$matrices$V))   V_truth <- truth$matrices$V
    if (is.null(V_truth) && !is.null(truth$matrices$V_m)) V_truth <- truth$matrices$V_m
  }
  if (is.null(V_truth)) V_truth <- build_truth_V(A_true, C_true, E_true, B_true)
  
  rP_true   <- if (!is.null(V_truth)) safe_div(V_truth[2,1], sqrt(V_truth[1,1]*V_truth[2,2])) else NA_real_
  sb_true   <- if (!is.null(V_truth)) std_betas_from_V(B_true, V_truth) else c(sb12=NA_real_, sb21=NA_real_)
  totV_true <- if (!is.null(V_truth)) sum(diag(V_truth)) else NA_real_
  
  if (!is.null(truth$realized)) {
    rA_true <- if (!is.null(truth$realized$rA_m)) truth$realized$rA_m else rA_true
    rC_true <- if (!is.null(truth$realized$rC_m)) truth$realized$rC_m else rC_true
    rE_true <- if (!is.null(truth$realized$rE_m)) truth$realized$rE_m else rE_true
  }
  
  truth_std <- c(
    Pva11 = Pva11_true, Pvc11 = Pvc11_true, Pve11 = Pve11_true,
    Pva22 = Pva22_true, Pvc22 = Pvc22_true, Pve22 = Pve22_true,
    rA21  = rA_true,    rC21  = rC_true,    rE21  = rE_true,
    sb12  = sb_true["sb12"], sb21 = sb_true["sb21"],
    rP21  = rP_true,
    totV  = totV_true
  )
  
  map_prop_names <- c(
    "Prop[1,1]"="Pva11","Prop[1,2]"="Pvc11","Prop[1,3]"="Pve11",
    "Prop[1,4]"="Pva22","Prop[1,5]"="Pvc22","Prop[1,6]"="Pve22",
    "Prop[1,7]"="rA21", "Prop[1,8]"="rC21","Prop[1,9]"="rE21",
    "Prop[1,10]"="sb21","Prop[1,11]"="sb12","Prop[1,12]"="rP21","Prop[1,13]"="totV"
  )
  
  # Extract fitted parameters
  fit_params <- single_fit$output$estimate
  
  # Helper function to compute Prop vector from parameter vector
  compute_prop_from_params <- function(param_vals) {
    # Extract variance components (diagonal elements)
    va11_p <- if ("va11" %in% names(param_vals)) param_vals["va11"] else NA_real_
    va22_p <- if ("va22" %in% names(param_vals)) param_vals["va22"] else NA_real_
    vc11_p <- if ("vc11" %in% names(param_vals)) param_vals["vc11"] else NA_real_
    vc22_p <- if ("vc22" %in% names(param_vals)) param_vals["vc22"] else NA_real_
    ve11_p <- if ("ve11" %in% names(param_vals)) param_vals["ve11"] else NA_real_
    ve22_p <- if ("ve22" %in% names(param_vals)) param_vals["ve22"] else NA_real_
    
    # Extract correlations
    rA12_p <- if ("rA12" %in% names(param_vals)) param_vals["rA12"] else 0
    rC12_p <- if ("rC12" %in% names(param_vals)) param_vals["rC12"] else 0
    rE12_p <- if ("rE12" %in% names(param_vals)) param_vals["rE12"] else 0
    
    # Extract B matrix elements
    b12_p <- if ("b12" %in% names(param_vals)) param_vals["b12"] else 0
    b21_p <- if ("b21" %in% names(param_vals)) param_vals["b21"] else 0
    
    # Reconstruct A, C, E matrices
    A_p <- matrix(c(va11_p, rA12_p * sqrt(va11_p * va22_p),
                    rA12_p * sqrt(va11_p * va22_p), va22_p), nrow=2, byrow=TRUE)
    C_p <- matrix(c(vc11_p, rC12_p * sqrt(vc11_p * vc22_p),
                    rC12_p * sqrt(vc11_p * vc22_p), vc22_p), nrow=2, byrow=TRUE)
    E_p <- matrix(c(ve11_p, rE12_p * sqrt(ve11_p * ve22_p),
                    rE12_p * sqrt(ve11_p * ve22_p), ve22_p), nrow=2, byrow=TRUE)
    
    # Reconstruct B matrix
    B_p <- matrix(c(0, b12_p, b21_p, 0), nrow=2, byrow=TRUE)
    
    # Compute (I - B)^{-1}
    I2 <- diag(2)
    IB_p <- I2 - B_p
    if (any(!is.finite(IB_p)) || abs(det(IB_p)) < .Machine$double.eps) {
      return(rep(NA_real_, 13))
    }
    invIB_p <- solve(IB_p)
    
    # Compute V = (I-B)^{-1} (A + C + E) (I-B)^{-T}
    V_p <- invIB_p %*% (A_p + C_p + E_p) %*% t(invIB_p)
    
    # Compute proportions
    Pva11_p <- safe_div(va11_p, va11_p + vc11_p + ve11_p)
    Pvc11_p <- safe_div(vc11_p, va11_p + vc11_p + ve11_p)
    Pve11_p <- safe_div(ve11_p, va11_p + vc11_p + ve11_p)
    Pva22_p <- safe_div(va22_p, va22_p + vc22_p + ve22_p)
    Pvc22_p <- safe_div(vc22_p, va22_p + vc22_p + ve22_p)
    Pve22_p <- safe_div(ve22_p, va22_p + vc22_p + ve22_p)
    
    # Compute correlations
    rA21_p <- safe_div(A_p[2,1], sqrt(A_p[1,1] * A_p[2,2]))
    rC21_p <- safe_div(C_p[2,1], sqrt(C_p[1,1] * C_p[2,2]))
    rE21_p <- safe_div(E_p[2,1], sqrt(E_p[1,1] * E_p[2,2]))
    
    # Compute standardized betas
    iSD_p <- diag(1 / sqrt(diag(V_p)))
    sb21_p <- iSD_p[2,2] * B_p[2,1] / iSD_p[1,1]
    sb12_p <- iSD_p[1,1] * B_p[1,2] / iSD_p[2,2]
    
    # Compute phenotypic correlation
    rP21_p <- safe_div(V_p[2,1], sqrt(V_p[1,1] * V_p[2,2]))
    
    # Compute total variance
    totV_p <- sum(diag(V_p))
    
    c(
      Pva11_p, Pvc11_p, Pve11_p,
      Pva22_p, Pvc22_p, Pve22_p,
      rA21_p, rC21_p, rE21_p,
      sb21_p, sb12_p, rP21_p, totV_p
    )
  }
  
  # Compute Prop vector
  prop_vec <- compute_prop_from_params(fit_params)
  
  # Compute delta-method SEs
  eps_base <- 1e-6
  coef_vals <- coef(single_fit)
  vcov_mat <- vcov(single_fit)
  
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
  
  # Compute performance metrics
  performance_metrics <- within(estimates_df, {
    bias      <- estimate - truth
    MAE       <- abs(bias)
    MSE      <- bias^2
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
  
  fit3 <- list(
    Face = fitFace, Haco = fitHaco, Hacr = fitHacr, Haeo = fitHaeo, Haer = fitHaer,
    Haor = fitHaor, Hceo = fitHceo, Hcer = fitHcer, Hcor = fitHcor, Heor = fitHeor
  )
  fit2 <- list(
    Fac = fitFac, Fae = fitFae, Fce = fitFce, Hcr = fitHcr, Hco = fitHco,
    Har = fitHar, Hao = fitHao, Her = fitHer, Heo = fitHeo, Dor = fitDor
  )
  fit1 <- list(
    Fa = fitFa, Fc = fitFc, Fe = fitFe, Do = fitDo, Dr = fitDr
  )
  fits <- list(
    Fa = fitFa, Fc = fitFc, Fe = fitFe, Do = fitDo, Dr = fitDr,
    Fac = fitFac, Fae = fitFae, Fce = fitFce, Hcr = fitHcr, Hco = fitHco,
    Har = fitHar, Hao = fitHao, Her = fitHer, Heo = fitHeo, Dor = fitDor,
    Face = fitFace, Haco = fitHaco, Hacr = fitHacr, Haeo = fitHaeo, Haer = fitHaer,
    Haor = fitHaor, Hceo = fitHceo, Hcer = fitHcer, Hcor = fitHcor, Heor = fitHeor
  )
  
  
  fits <- list(fit = single_model_select(fit1, fit2, fit3, fits, fitFace)$single_fit)
  
  single_fit <- fits[[1]]
  results_df <- process_all_models(fits)
  
  write.csv(results_df,paste0(output_prefix, "m1_models.csv"), row.names = FALSE)
  
  use_vp <- if (isTRUE(truth$flags$sex_differences)) truth$varparam_m else truth$varparam
  if (is.null(use_vp)) use_vp <- truth$varparam
  
  truth_raw <- c(
    va11 = use_vp$va11, va22 = use_vp$va22, va12 = use_vp$va12,
    vc11 = use_vp$vc11, vc22 = use_vp$vc22, vc12 = use_vp$vc12,
    ve11 = use_vp$ve11, ve22 = use_vp$ve22, ve12 = use_vp$ve12,
    sb12  = use_vp$b12,  sb21  = use_vp$b21
  )
  
  Pva11_true <- safe_div(use_vp$va11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pvc11_true <- safe_div(use_vp$vc11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pve11_true <- safe_div(use_vp$ve11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pva22_true <- safe_div(use_vp$va22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  Pvc22_true <- safe_div(use_vp$vc22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  Pve22_true <- safe_div(use_vp$ve22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  
  rA_true <- safe_div(use_vp$va12, sqrt(use_vp$va11 * use_vp$va22))
  rC_true <- safe_div(use_vp$vc12, sqrt(use_vp$vc11 * use_vp$vc22))
  rE_true <- safe_div(use_vp$ve12, sqrt(use_vp$ve11 * use_vp$ve22))
  
  A_true <- mat2(use_vp$va11, use_vp$va22, use_vp$va12)
  C_true <- mat2(use_vp$vc11, use_vp$vc22, use_vp$vc12)
  E_true <- mat2(use_vp$ve11, use_vp$ve22, use_vp$ve12)
  B_true <- matrix(c(0, use_vp$b12, use_vp$b21, 0), nrow=2, byrow=TRUE)
  
  build_truth_V <- function(A, C, E, B) {
    I2 <- diag(2)
    IB <- I2 - B
    if (any(!is.finite(IB)) || abs(det(IB)) < .Machine$double.eps) return(NULL)
    invIB <- solve(IB)
    invIB %*% (A + C + E) %*% t(invIB)
  }
  
  std_betas_from_V <- function(B, V) {
    if (is.null(V)) return(c(sb12 = NA_real_, sb21 = NA_real_))
    iSD <- diag(1 / sqrt(diag(V)))
    c(sb21 = iSD[2,2] * B[2,1] / iSD[1,1],
      sb12 = iSD[1,1] * B[1,2] / iSD[2,2])
  }
  
  V_truth <- NULL
  if (isTRUE(truth$flags$sex_differences)) {
    if (!is.null(truth$matrices$V_m)) V_truth <- truth$matrices$V_m
  } else {
    if (!is.null(truth$matrices$V))   V_truth <- truth$matrices$V
    if (is.null(V_truth) && !is.null(truth$matrices$V_m)) V_truth <- truth$matrices$V_m
  }
  if (is.null(V_truth)) V_truth <- build_truth_V(A_true, C_true, E_true, B_true)
  
  rP_true   <- if (!is.null(V_truth)) safe_div(V_truth[2,1], sqrt(V_truth[1,1]*V_truth[2,2])) else NA_real_
  sb_true   <- if (!is.null(V_truth)) std_betas_from_V(B_true, V_truth) else c(sb12=NA_real_, sb21=NA_real_)
  totV_true <- if (!is.null(V_truth)) sum(diag(V_truth)) else NA_real_
  
  if (!is.null(truth$realized)) {
    rA_true <- if (!is.null(truth$realized$rA_m)) truth$realized$rA_m else rA_true
    rC_true <- if (!is.null(truth$realized$rC_m)) truth$realized$rC_m else rC_true
    rE_true <- if (!is.null(truth$realized$rE_m)) truth$realized$rE_m else rE_true
  }
  
  truth_std <- c(
    Pva11 = Pva11_true, Pvc11 = Pvc11_true, Pve11 = Pve11_true,
    Pva22 = Pva22_true, Pvc22 = Pvc22_true, Pve22 = Pve22_true,
    rA21  = rA_true,    rC21  = rC_true,    rE21  = rE_true,
    sb12  = sb_true["sb12"], sb21 = sb_true["sb21"],
    rP21  = rP_true,
    totV  = totV_true
  )
  
  map_prop_names <- c(
    "Prop[1,1]"="Pva11","Prop[1,2]"="Pvc11","Prop[1,3]"="Pve11",
    "Prop[1,4]"="Pva22","Prop[1,5]"="Pvc22","Prop[1,6]"="Pve22",
    "Prop[1,7]"="rA21", "Prop[1,8]"="rC21","Prop[1,9]"="rE21",
    "Prop[1,10]"="sb21","Prop[1,11]"="sb12","Prop[1,12]"="rP21","Prop[1,13]"="totV"
  )
  
  # Extract fitted parameters
  fit_params <- coef(single_fit)
  
  # Helper function to compute Prop vector from parameter vector
  compute_prop_from_params <- function(param_vals) {
    # Extract variance components (diagonal elements)
    va11_p <- if ("va11" %in% names(param_vals)) param_vals["va11"] else NA_real_
    va22_p <- if ("va22" %in% names(param_vals)) param_vals["va22"] else NA_real_
    vc11_p <- if ("vc11" %in% names(param_vals)) param_vals["vc11"] else NA_real_
    vc22_p <- if ("vc22" %in% names(param_vals)) param_vals["vc22"] else NA_real_
    ve11_p <- if ("ve11" %in% names(param_vals)) param_vals["ve11"] else NA_real_
    ve22_p <- if ("ve22" %in% names(param_vals)) param_vals["ve22"] else NA_real_
    
    # Extract correlations
    rA12_p <- if ("rA12" %in% names(param_vals)) param_vals["rA12"] else 0
    rC12_p <- if ("rC12" %in% names(param_vals)) param_vals["rC12"] else 0
    rE12_p <- if ("rE12" %in% names(param_vals)) param_vals["rE12"] else 0
    
    # Extract B matrix elements
    b12_p <- if ("b12" %in% names(param_vals)) param_vals["b12"] else 0
    b21_p <- if ("b21" %in% names(param_vals)) param_vals["b21"] else 0
    
    # Reconstruct A, C, E matrices with guarded square roots
    A_p <- matrix(c(va11_p, rA12_p * sqrt(pmax(va11_p * va22_p, 0)),
                    rA12_p * sqrt(pmax(va11_p * va22_p, 0)), va22_p), nrow=2, byrow=TRUE)
    C_p <- matrix(c(vc11_p, rC12_p * sqrt(pmax(vc11_p * vc22_p, 0)),
                    rC12_p * sqrt(pmax(vc11_p * vc22_p, 0)), vc22_p), nrow=2, byrow=TRUE)
    E_p <- matrix(c(ve11_p, rE12_p * sqrt(pmax(ve11_p * ve22_p, 0)),
                    rE12_p * sqrt(pmax(ve11_p * ve22_p, 0)), ve22_p), nrow=2, byrow=TRUE)
    
    # Reconstruct B matrix
    B_p <- matrix(c(0, b12_p, b21_p, 0), nrow=2, byrow=TRUE)
    
    # Compute (I - B)^{-1}
    I2 <- diag(2)
    IB_p <- I2 - B_p
    if (any(!is.finite(IB_p)) || rcond(IB_p) < 1e-12) {
      return(rep(NA_real_, 13))
    }
    invIB_p <- solve(IB_p)
    
    # Compute V = (I-B)^{-1} (A + C + E) (I-B)^{-T}
    V_p <- invIB_p %*% (A_p + C_p + E_p) %*% t(invIB_p)
    
    # Compute proportions
    Pva11_p <- safe_div(va11_p, va11_p + vc11_p + ve11_p)
    Pvc11_p <- safe_div(vc11_p, va11_p + vc11_p + ve11_p)
    Pve11_p <- safe_div(ve11_p, va11_p + vc11_p + ve11_p)
    Pva22_p <- safe_div(va22_p, va22_p + vc22_p + ve22_p)
    Pvc22_p <- safe_div(vc22_p, va22_p + vc22_p + ve22_p)
    Pve22_p <- safe_div(ve22_p, va22_p + vc22_p + ve22_p)
    
    # Compute correlations
    rA21_p <- safe_div(A_p[2,1], sqrt(pmax(A_p[1,1] * A_p[2,2], 0)))
    rC21_p <- safe_div(C_p[2,1], sqrt(pmax(C_p[1,1] * C_p[2,2], 0)))
    rE21_p <- safe_div(E_p[2,1], sqrt(pmax(E_p[1,1] * E_p[2,2], 0)))
    
    # Compute standardized betas with guarded diagonal extraction
    dV <- diag(V_p)
    if (any(!is.finite(dV)) || any(dV <= 0)) return(rep(NA_real_, 13))
    iSD_p <- diag(1 / sqrt(dV))
    sb21_p <- iSD_p[2,2] * B_p[2,1] / iSD_p[1,1]
    sb12_p <- iSD_p[1,1] * B_p[1,2] / iSD_p[2,2]
    
    # Compute phenotypic correlation
    rP21_p <- safe_div(V_p[2,1], sqrt(pmax(V_p[1,1] * V_p[2,2], 0)))
    
    # Compute total variance
    totV_p <- sum(diag(V_p))
    
    c(
      Pva11_p, Pvc11_p, Pve11_p,
      Pva22_p, Pvc22_p, Pve22_p,
      rA21_p, rC21_p, rE21_p,
      sb21_p, sb12_p, rP21_p, totV_p
    )
  }
  
  # Compute delta-method SEs
  eps_base <- 1e-5
  coef_vals <- coef(single_fit)
  vcov_mat  <- vcov(single_fit)
  prop_vec  <- compute_prop_from_params(coef_vals)
  
  if (is.null(coef_vals) || length(coef_vals) == 0 || 
      is.null(vcov_mat) || !all(is.finite(vcov_mat))) {
    prop_ses <- rep(NA_real_, length(prop_vec))
  } else {
    # Ensure vcov_mat matches the order of coef_vals
    if (!is.null(colnames(vcov_mat)) && !is.null(names(coef_vals))) {
      vcov_mat <- vcov_mat[names(coef_vals), names(coef_vals), drop = FALSE]
    }
    
    n_params <- length(coef_vals)
    n_props <- length(prop_vec)
    gradients <- matrix(NA_real_, nrow = n_props, ncol = n_params)
    
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
        # Central difference
        if (is.finite(prop_plus[k]) && is.finite(prop_minus[k])) {
          gradients[k, j] <- (prop_plus[k] - prop_minus[k]) / (2 * h)
          
        } else {
          # Forward difference
          if (is.finite(prop_plus[k]) && is.finite(prop_vec[k])) {
            gradients[k, j] <- (prop_plus[k] - prop_vec[k]) / h
            
            # Backward difference
          } else if (is.finite(prop_minus[k]) && is.finite(prop_vec[k])) {
            gradients[k, j] <- (prop_vec[k] - prop_minus[k]) / h
            
            # Derivative cannot be computed reliably
          } else {
            gradients[k, j] <- NA_real_
          }
        }
      }
    }
    
    # Compute delta-method SEs: sqrt(g' * vcov * g)
    prop_ses <- apply(gradients, 1, function(g) {
      if (anyNA(g)) return(NA_real_)
      se_sq <- as.numeric(g %*% vcov_mat %*% g)
      sqrt(pmax(se_sq, 0))
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
  
  # Compute performance metrics
  tol0 <- 1e-12
  performance_metrics <- within(estimates_df, {
    bias      <- estimate - truth
    MAE       <- abs(bias)
    MSE       <- bias^2
    covered95 <- ifelse(!is.na(truth) & (truth >= lcl95) & (truth <= ucl95), 1L, 0L)
    typeI     <- ifelse(!is.na(truth) & abs(truth) < tol0 & (lcl95 > 0 | ucl95 < 0), 1L, 0L)
    typeII    <- ifelse(!is.na(truth) & abs(truth) >= tol0 & (lcl95 <= 0 & ucl95 >= 0), 1L, 0L)
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
  # METHOD 2 – Bivariate Multimodel Averaging
  # ========================================================================
  
  fits <- list(
    Fa = fitFa, Fc = fitFc, Fe = fitFe, Do = fitDo, Dr = fitDr,
    Fac = fitFac, Fae = fitFae, Fce = fitFce, Hcr = fitHcr, Hco = fitHco,
    Har = fitHar, Hao = fitHao, Her = fitHer, Heo = fitHeo, Dor = fitDor,
    Face = fitFace, Haco = fitHaco, Hacr = fitHacr, Haeo = fitHaeo, Haer = fitHaer,
    Haor = fitHaor, Hceo = fitHceo, Hcer = fitHcer, Hcor = fitHcor, Heor = fitHeor
  )
  
  results_df <- process_all_models(fits)
  
  write.csv(results_df,paste0(output_prefix, "m2_models.csv"), row.names = FALSE)
  
  fitACE5_aug_m2 <- lapply(fits, add_prop_block_biv)
  
  colP_all_m2 <- c("va11","va22","rA12",
                   "vc11","vc22","rC12", 
                   "ve11","ve22","rE12",
                   "b21","b12")
  
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
    va11 = use_vp$va11, va22 = use_vp$va22, va12 = use_vp$va12,
    vc11 = use_vp$vc11, vc22 = use_vp$vc22, vc12 = use_vp$vc12,
    ve11 = use_vp$ve11, ve22 = use_vp$ve22, ve12 = use_vp$ve12,
    sb12  = use_vp$b12,  sb21  = use_vp$b21
  )
  
  Pva11_true <- safe_div(use_vp$va11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pvc11_true <- safe_div(use_vp$vc11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pve11_true <- safe_div(use_vp$ve11, use_vp$va11 + use_vp$vc11 + use_vp$ve11)
  Pva22_true <- safe_div(use_vp$va22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  Pvc22_true <- safe_div(use_vp$vc22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  Pve22_true <- safe_div(use_vp$ve22, use_vp$va22 + use_vp$vc22 + use_vp$ve22)
  
  rA_true <- safe_div(use_vp$va12, sqrt(use_vp$va11 * use_vp$va22))
  rC_true <- safe_div(use_vp$vc12, sqrt(use_vp$vc11 * use_vp$vc22))
  rE_true <- safe_div(use_vp$ve12, sqrt(use_vp$ve11 * use_vp$ve22))
  
  A_true <- mat2(use_vp$va11, use_vp$va22, use_vp$va12)
  C_true <- mat2(use_vp$vc11, use_vp$vc22, use_vp$vc12)
  E_true <- mat2(use_vp$ve11, use_vp$ve22, use_vp$ve12)
  B_true <- matrix(c(0, use_vp$b12, use_vp$b21, 0), nrow=2, byrow=TRUE)
  
  build_truth_V <- function(A, C, E, B) {
    I2 <- diag(2)
    IB <- I2 - B
    if (any(!is.finite(IB)) || abs(det(IB)) < .Machine$double.eps) return(NULL)
    invIB <- solve(IB)
    invIB %*% (A + C + E) %*% t(invIB)
  }
  
  std_betas_from_V <- function(B, V) {
    if (is.null(V)) return(c(sb12 = NA_real_, sb21 = NA_real_))
    diag_V <- diag(V)
    if (any(diag_V <= 0)) return(c(sb12 = NA_real_, sb21 = NA_real_))
    iSD <- diag(1 / sqrt(diag_V))
    c(sb21 = iSD[2,2] * B[2,1] / iSD[1,1],
      sb12 = iSD[1,1] * B[1,2] / iSD[2,2])
  }
  
  V_truth <- NULL
  if (isTRUE(truth$flags$sex_differences)) {
    if (!is.null(truth$matrices$V_m)) V_truth <- truth$matrices$V_m
  } else {
    if (!is.null(truth$matrices$V))   V_truth <- truth$matrices$V
    if (is.null(V_truth) && !is.null(truth$matrices$V_m)) V_truth <- truth$matrices$V_m
  }
  if (is.null(V_truth)) V_truth <- build_truth_V(A_true, C_true, E_true, B_true)
  
  rP_true   <- if (!is.null(V_truth)) safe_div(V_truth[2,1], sqrt(V_truth[1,1]*V_truth[2,2])) else NA_real_
  sb_true   <- if (!is.null(V_truth)) std_betas_from_V(B_true, V_truth) else c(sb12=NA_real_, sb21=NA_real_)
  totV_true <- if (!is.null(V_truth)) sum(diag(V_truth)) else NA_real_
  
  if (!is.null(truth$realized)) {
    rA_true <- if (!is.null(truth$realized$rA_m)) truth$realized$rA_m else rA_true
    rC_true <- if (!is.null(truth$realized$rC_m)) truth$realized$rC_m else rC_true
    rE_true <- if (!is.null(truth$realized$rE_m)) truth$realized$rE_m else rE_true
  }
  
  truth_std <- c(
    Pva11 = Pva11_true, Pvc11 = Pvc11_true, Pve11 = Pve11_true,
    Pva22 = Pva22_true, Pvc22 = Pvc22_true, Pve22 = Pve22_true,
    rA21  = rA_true,    rC21  = rC_true,    rE21  = rE_true,
    sb12  = sb_true["sb12"], sb21 = sb_true["sb21"],
    rP21  = rP_true,
    totV  = totV_true
  )
  
  map_prop_names <- c(
    "Prop[1,1]"="Pva11","Prop[1,2]"="Pvc11","Prop[1,3]"="Pve11",
    "Prop[1,4]"="Pva22","Prop[1,5]"="Pvc22","Prop[1,6]"="Pve22",
    "Prop[1,7]"="rA21", "Prop[1,8]"="rC21","Prop[1,9]"="rE21",
    "Prop[1,10]"="sb21","Prop[1,11]"="sb12","Prop[1,12]"="rP21","Prop[1,13]"="totV"
  )
  
  ma_tab <- mma_var_m2$`Model-Average Estimates`
  ma_tab <- mma_var_m2$`Model-Average Estimates`
  fit_params_est <- setNames(as.numeric(ma_tab[, "Estimate"]), rownames(ma_tab))
  fit_params_SE  <- setNames(as.numeric(ma_tab[, "SE"]),       rownames(ma_tab))
  
  coef_vals <- fit_params_est
  se_vec    <- fit_params_SE
  
  # Helper function to compute Prop vector from parameter vector
  compute_prop_from_params <- function(param_vals) {
    # Extract variance components (diagonal elements)
    va11_p <- if ("va11" %in% names(param_vals)) param_vals["va11"] else NA_real_
    va22_p <- if ("va22" %in% names(param_vals)) param_vals["va22"] else NA_real_
    vc11_p <- if ("vc11" %in% names(param_vals)) param_vals["vc11"] else NA_real_
    vc22_p <- if ("vc22" %in% names(param_vals)) param_vals["vc22"] else NA_real_
    ve11_p <- if ("ve11" %in% names(param_vals)) param_vals["ve11"] else NA_real_
    ve22_p <- if ("ve22" %in% names(param_vals)) param_vals["ve22"] else NA_real_
    
    # Extract correlations (treat missing as 0)
    rA12_p <- if ("rA12" %in% names(param_vals)) param_vals["rA12"] else 0
    rC12_p <- if ("rC12" %in% names(param_vals)) param_vals["rC12"] else 0
    rE12_p <- if ("rE12" %in% names(param_vals)) param_vals["rE12"] else 0
    
    # Extract B matrix elements (treat missing as 0)
    b12_p <- if ("b12" %in% names(param_vals)) param_vals["b12"] else 0
    b21_p <- if ("b21" %in% names(param_vals)) param_vals["b21"] else 0
    
    # Check for valid variances
    if (!all(is.finite(c(va11_p, va22_p, vc11_p, vc22_p, ve11_p, ve22_p)))) {
      return(rep(NA_real_, 13))
    }
    if (any(c(va11_p, va22_p, vc11_p, vc22_p, ve11_p, ve22_p) < 0)) {
      return(rep(NA_real_, 13))
    }
    
    # Reconstruct A, C, E matrices using rsqrt(v11*v22) for off-diagonals
    va12_p <- if (va11_p > 0 && va22_p > 0) rA12_p * sqrt(va11_p * va22_p) else 0
    vc12_p <- if (vc11_p > 0 && vc22_p > 0) rC12_p * sqrt(vc11_p * vc22_p) else 0
    ve12_p <- if (ve11_p > 0 && ve22_p > 0) rE12_p * sqrt(ve11_p * ve22_p) else 0
    
    A_p <- matrix(c(va11_p, va12_p, va12_p, va22_p), nrow=2, byrow=TRUE)
    C_p <- matrix(c(vc11_p, vc12_p, vc12_p, vc22_p), nrow=2, byrow=TRUE)
    E_p <- matrix(c(ve11_p, ve12_p, ve12_p, ve22_p), nrow=2, byrow=TRUE)
    
    # Reconstruct B matrix
    B_p <- matrix(c(0, b12_p, b21_p, 0), nrow=2, byrow=TRUE)
    
    # Compute (I - B)^{-1}
    I2 <- diag(2)
    IB_p <- I2 - B_p
    if (any(!is.finite(IB_p)) || abs(det(IB_p)) < .Machine$double.eps) {
      return(rep(NA_real_, 13))
    }
    invIB_p <- solve(IB_p)
    
    # Compute V = (I-B)^{-1} (A + C + E) (I-B)^{-T}
    V_p <- invIB_p %*% (A_p + C_p + E_p) %*% t(invIB_p)
    
    if (!all(is.finite(V_p))) {
      return(rep(NA_real_, 13))
    }
    
    # Compute proportions
    Pva11_p <- safe_div(va11_p, va11_p + vc11_p + ve11_p)
    Pvc11_p <- safe_div(vc11_p, va11_p + vc11_p + ve11_p)
    Pve11_p <- safe_div(ve11_p, va11_p + vc11_p + ve11_p)
    Pva22_p <- safe_div(va22_p, va22_p + vc22_p + ve22_p)
    Pvc22_p <- safe_div(vc22_p, va22_p + vc22_p + ve22_p)
    Pve22_p <- safe_div(ve22_p, va22_p + vc22_p + ve22_p)
    
    # Compute correlations from A, C, E matrices
    rA21_p <- safe_div(A_p[2,1], sqrt(A_p[1,1] * A_p[2,2]))
    rC21_p <- safe_div(C_p[2,1], sqrt(C_p[1,1] * C_p[2,2]))
    rE21_p <- safe_div(E_p[2,1], sqrt(E_p[1,1] * E_p[2,2]))
    
    # Compute standardized betas from V
    diag_V_p <- diag(V_p)
    if (any(diag_V_p <= 0)) {
      sb21_p <- NA_real_
      sb12_p <- NA_real_
    } else {
      iSD_p <- diag(1 / sqrt(diag_V_p))
      sb21_p <- iSD_p[2,2] * B_p[2,1] / iSD_p[1,1]
      sb12_p <- iSD_p[1,1] * B_p[1,2] / iSD_p[2,2]
    }
    
    # Compute phenotypic correlation
    
    W <- mma_var_m2$`Akaike-Weights Table`
    w <- setNames(W$AkaikeWeight, W$model)
    
    mw <- mma_var_m2$`Model-wise Estimates`  # params x models
    
    calc_rP_from_col <- function(parcol) {
      # parcol is a named vector like c(va11=..., va22=..., rA12=..., ...)
      get <- function(nm, default=0) if (nm %in% names(parcol) && is.finite(parcol[nm])) parcol[nm] else default
      
      va11 <- get("va11", NA); va22 <- get("va22", NA)
      vc11 <- get("vc11", NA); vc22 <- get("vc22", NA)
      ve11 <- get("ve11", NA); ve22 <- get("ve22", NA)
      if (any(!is.finite(c(va11,va22,vc11,vc22,ve11,ve22)))) return(NA_real_)
      
      rA12 <- get("rA12", 0); rC12 <- get("rC12", 0); rE12 <- get("rE12", 0)
      b21  <- get("b21",  0); b12  <- get("b12",  0)
      
      va12 <- rA12 * sqrt(va11*va22)
      vc12 <- rC12 * sqrt(vc11*vc22)
      ve12 <- rE12 * sqrt(ve11*ve22)
      
      S <- matrix(c(va11+vc11+ve11, va12+vc12+ve12,
                    va12+vc12+ve12, va22+vc22+ve22), 2,2, byrow=TRUE)
      
      IB <- diag(2) - matrix(c(0,b12,b21,0),2,2,byrow=TRUE)
      if (abs(det(IB)) < .Machine$double.eps) return(NA_real_)
      invIB <- solve(IB)
      
      V <- invIB %*% S %*% t(invIB)
      V[2,1] / sqrt(V[1,1]*V[2,2])
    }
    
    rP_by_model <- sapply(colnames(mw), function(m) calc_rP_from_col(mw[, m]))
    rP21_p <- sum(w[names(rP_by_model)] * rP_by_model, na.rm=TRUE)    
    
    # Compute total variance
    totV_p <- sum(diag(V_p))
    
    c(
      Pva11_p, Pvc11_p, Pve11_p,
      Pva22_p, Pvc22_p, Pve22_p,
      rA21_p, rC21_p, rE21_p,
      sb21_p, sb12_p, rP21_p, totV_p
    )
  }
  
  # Compute Prop vector
  prop_vec <- compute_prop_from_params(fit_params_est)
  
  # Compute delta-method SEs
  coef_vals <- fit_params_est
  se_vec <- fit_params_SE  # you need to keep the SE vector; see below
  vcov_mat <- diag(se_vec^2)
  rownames(vcov_mat) <- colnames(vcov_mat) <- names(se_vec)
  
  # Compute delta-method SEs
  eps_base <- 1e-6
  
  ma_tab <- mma_var_m2$`Model-Average Estimates`
  coef_vals <- setNames(as.numeric(ma_tab[, "Estimate"]), rownames(ma_tab))
  se_vec    <- setNames(as.numeric(ma_tab[, "SE"]),       rownames(ma_tab))
  
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
    
    # If vcov diagonal has NAs, the corresponding gradients contribute NA variance
    # We'll handle that by zeroing gradient entries where variance is NA (conservative)
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
    
    # zero out gradient columns for parameters with missing variance
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
# 10) Run the simulation and analysis with infinite retry per replication
#     (no max attempts; never duplicates seeds)
# ------------------------------------------------------------------------------

n_reps <- 1000

# Optional: Create a directory for all runs
dir.create("simulation_results_Face", showWarnings = FALSE)

# Use a single global seed counter that increments on every attempt
seed_i <- 0L

for (i in seq_len(n_reps)) {
  attempt <- 0L
  repeat {
    # Advance seed for every single attempt (success OR failure)
    seed_i  <- seed_i + 1L
    attempt <- attempt + 1L
    
    message(sprintf("Running replication %d of %d (attempt %d) with seed %d",
                    i, n_reps, attempt, seed_i))
    
    # Try one full simulation + analysis
    res <- tryCatch(
      {
        run_simulation_and_analysis(seed = seed_i)
      },
      error = function(e) {
        # Log and retry immediately with next seed
        message(sprintf("  -> attempt %d failed: %s", attempt, conditionMessage(e)))
        NULL
      }
    )
    
    # Break the repeat loop on success (accepted dataset)
    if (!is.null(res)) break
  }
  
  # Build unique file prefix for saving this replication
  prefix <- sprintf("simulation_results_Face/run_%04d", i)
  
  # Save the returned data and results with unique names
  # NOTE: ensure run_simulation_and_analysis returns these components.
  # If not, adapt the saving accordingly.
  if (!is.null(res$data)) {
    utils::write.csv(res$data, paste0(prefix, "_data.csv"), row.names = FALSE)
  }
  res$performance_metrics$seed <- seed_i
  utils::write.csv(res$model_fits_m2,            paste0(prefix, "_model_fits_m2.csv"),        row.names = FALSE)
  utils::write.csv(res$averaged_estimates,    paste0(prefix, "_model_averaged.csv"),    row.names = FALSE)
  utils::write.csv(res$performance_metrics,   paste0(prefix, "_performance.csv"),       row.names = FALSE)
  saveRDS(res, paste0(prefix, "_full.RDS"))
}
