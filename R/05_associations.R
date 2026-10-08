assoc_setting <- function(cfg, name, default) {
  if (is.null(cfg[[name]])) default else cfg[[name]]
}

assoc_rcs_spec <- function(x, probs = c(0.05, 0.35, 0.65, 0.95)) {
  mu <- mean(x)
  sig <- stats::sd(x)
  if (!is.finite(sig) || sig <= 0) stop("Association exposure has zero variance.")
  z <- (x - mu) / sig
  kk <- as.numeric(stats::quantile(z, probs, names = FALSE, type = 7))
  if (length(kk) != 4L || any(diff(kk) <= 0)) {
    stop("Association exposure requires four distinct quantile knots; inspect ties.")
  }
  list(center = mu, scale = sig, knots_z = kk,
       knots_raw = mu + sig * kk, knot_probs = probs,
       terms = c("assoc_x1", "assoc_xnl1", "assoc_xnl2"))
}

assoc_rcs_matrix <- function(x, spec) {
  z <- (x - spec$center) / spec$scale
  k <- spec$knots_z
  den <- (k[4L] - k[1L])^2
  H <- vapply(seq_len(2L), function(j) {
    (pmax(z - k[j], 0)^3 -
       pmax(z - k[3L], 0)^3 * (k[4L] - k[j]) / (k[4L] - k[3L]) +
       pmax(z - k[4L], 0)^3 * (k[3L] - k[j]) / (k[4L] - k[3L])) / den
  }, numeric(length(z)))
  ans <- cbind(z, matrix(H, nrow = length(z), ncol = 2L))
  colnames(ans) <- spec$terms
  ans
}

assoc_prepare <- function(d, exposure, outcome, adjust_bmi, cfg) {
  needed <- unique(c("id", "side", "sex", "age", "FL", exposure, outcome,
                     if (adjust_bmi) "BMI"))
  missing <- setdiff(needed, names(d))
  if (length(missing)) stop("Association input missing: ", paste(missing, collapse = ", "))
  keep <- stats::complete.cases(d[, needed, drop = FALSE])
  numeric_vars <- needed[vapply(d[, needed, drop = FALSE], is.numeric, logical(1))]
  for (v in numeric_vars) keep <- keep & is.finite(d[[v]])
  z <- d[keep, , drop = FALSE]
  z <- z[order(z$id, z$side), , drop = FALSE]
  if (!nrow(z)) stop("No complete cases for ", exposure, " -> ", outcome)
  if (length(unique(z$sex)) != 2L || length(unique(z$side)) != 2L) {
    stop("Association models require both sex and recorded-side levels.")
  }
  if (any(z$FL <= 0)) stop("Nonpositive FL in association complete cases.")
  if (length(unique(z$id)) < assoc_setting(cfg, "min_clusters", 30L)) {
    stop("Too few independent children for association GEE.")
  }
  z$FL_10 <- (z$FL - 200) / 10
  if (adjust_bmi) z$BMI_5 <- (z$BMI - 18) / 5
  ak <- assoc_setting(cfg, "age_knots", c(6, 9, 12, 15))
  ab <- assoc_setting(cfg, "age_boundary", c(3, 19))
  if (length(ab) != 2L || any(ak <= ab[1L] | ak >= ab[2L])) {
    stop("Invalid fixed age-spline knots or boundaries in association config.")
  }
  if (any(z$age < ab[1L] | z$age > ab[2L])) {
    stop("Association ages extend outside configured age-spline boundaries.")
  }
  A <- splines::ns(z$age, knots = ak, Boundary.knots = ab, intercept = FALSE)
  age_names <- paste0("assoc_age", seq_len(ncol(A)))
  colnames(A) <- age_names
  z[age_names] <- as.data.frame(A)
  spec <- assoc_rcs_spec(z[[exposure]],
    assoc_setting(cfg, "association_exposure_knot_probs", c(.05, .35, .65, .95)))
  z[spec$terms] <- as.data.frame(assoc_rcs_matrix(z[[exposure]], spec))
  if (!is.finite(stats::sd(z[[outcome]])) || stats::sd(z[[outcome]]) <= 0) {
    stop("Association outcome has zero variance: ", outcome)
  }
  list(data = z, exposure_spec = spec, age_names = age_names,
       age_basis = A, n_input_feet = nrow(d), n_excluded_feet = sum(!keep),
       n_input_children = length(unique(d$id)),
       n_single_foot_children = sum(table(z$id) == 1L))
}

assoc_formula <- function(outcome, exposure_terms, age_names, adjust_bmi,
                          interaction = NULL) {
  rhs <- paste(c(paste(exposure_terms, collapse = " + "),
    paste0("(", paste(age_names, collapse = " + "), ") * sex"),
    "side", "FL_10", if (adjust_bmi) "BMI_5"), collapse = " + ")
  if (identical(interaction, "sex")) {
    rhs <- paste0(rhs, " + (", paste(exposure_terms, collapse = " + "), "):sex")
  } else if (identical(interaction, "age")) {
    rhs <- paste0(rhs, " + (", paste(exposure_terms, collapse = " + "),
                  "):(", paste(age_names, collapse = " + "), ")")
  }
  stats::as.formula(paste(outcome, "~", rhs))
}

assoc_wald_row <- function(obj, indices, test, model_id) {
  w <- wald_block(obj, indices)
  data.frame(model_id = model_id, test = test,
             statistic = as.numeric(w[["statistic"]]), df = as.numeric(w[["df"]]),
             p = as.numeric(w[["p"]]), status = "ok", note = "",
             n_children = obj$n_children, n_feet = obj$n_feet,
             stringsAsFactors = FALSE)
}

assoc_linear_coefficients <- function(obj, prep, exposure, outcome, model_id,
                                      analysis_set) {
  j <- match("assoc_x1", names(obj$beta))
  if (is.na(j)) stop("Missing association exposure coefficient.")
  sx <- prep$exposure_spec$scale
  sy <- stats::sd(prep$data[[outcome]])
  unit <- switch(exposure, BMI = 5, hallux = 5, AI = .01, AS_ratio = .01,
                 AV = 1000, AV_FL3 = .001, 1)
  mult <- c(raw_per_1 = 1 / sx, raw_per_increment = unit / sx,
            outcome_SD_per_increment = unit / (sx * sy),
            fully_standardised = 1 / sy)
  b <- unname(obj$beta[j])
  se <- sqrt(obj$V[j, j])
  ans <- data.frame(model_id = model_id, analysis_set = analysis_set,
    exposure = exposure, outcome = outcome, scale = names(mult),
    exposure_increment_raw = c(1, unit, unit, sx),
    estimate = b * mult, se = se * mult,
    lower = (b - stats::qnorm(.975) * se) * mult,
    upper = (b + stats::qnorm(.975) * se) * mult,
    p = 2 * stats::pnorm(-abs(b / se)),
    exposure_mean = prep$exposure_spec$center,
    exposure_sd = sx, outcome_mean = mean(prep$data[[outcome]]), outcome_sd = sy,
    n_children = obj$n_children, n_feet = obj$n_feet,
    stringsAsFactors = FALSE, row.names = NULL)
  ans
}

assoc_mean_design <- function(obj, data, spec, x) {
  X <- design(obj, data)
  L <- matrix(colMeans(X), nrow = length(x), ncol = ncol(X), byrow = TRUE,
              dimnames = list(NULL, colnames(X)))
  pos <- match(spec$terms, colnames(L))
  if (anyNA(pos)) stop("Exposure columns missing from association curve design.")
  L[, pos] <- assoc_rcs_matrix(x, spec)
  L
}

assoc_curve <- function(obj, prep, exposure, outcome, model_id, analysis_set,
                        cfg, seed) {
  support <- assoc_setting(cfg, "association_curve_support", c(.025, .975))
  if (length(support) != 2L || any(support < 0 | support > 1) || diff(support) <= 0) {
    stop("association_curve_support must be two increasing probabilities.")
  }
  xx <- prep$data[[exposure]]
  lim <- as.numeric(stats::quantile(xx, support, names = FALSE, type = 7))
  gx <- seq(lim[1L], lim[2L], length.out = assoc_setting(cfg, "curve_n", 101L))
  ref <- stats::median(xx)

  LL <- assoc_mean_design(obj, prep$data, prep$exposure_spec, c(gx, ref))
  L <- LL[seq_along(gx), , drop = FALSE]
  refL <- LL[nrow(LL), ]
  B <- as.data.frame(linear_bands(obj$beta, obj$V, L,
     n_sim = assoc_setting(cfg, "n_sim", 5000L), seed = seed))
  D <- sweep(L, 2L, refL, FUN = "-")
  DB <- as.data.frame(linear_bands(obj$beta, obj$V, D,
     n_sim = assoc_setting(cfg, "n_sim", 5000L), seed = seed + 1L))
  mu <- mean(prep$data[[outcome]])
  sy <- stats::sd(prep$data[[outcome]])
  ans <- data.frame(model_id = model_id, analysis_set = analysis_set,
    exposure = exposure, outcome = outcome, exposure_raw = gx,
    exposure_SD = (gx - prep$exposure_spec$center) / prep$exposure_spec$scale,
    estimate_raw = B$estimate, se_raw = B$se,
    lower_raw = B$lower, upper_raw = B$upper,
    simultaneous_lower_raw = B$sim_lower, simultaneous_upper_raw = B$sim_upper,
    estimate_outcome_SD = (B$estimate - mu) / sy,
    se_outcome_SD = B$se / sy,
    lower_outcome_SD = (B$lower - mu) / sy,
    upper_outcome_SD = (B$upper - mu) / sy,
    simultaneous_lower_outcome_SD = (B$sim_lower - mu) / sy,
    simultaneous_upper_outcome_SD = (B$sim_upper - mu) / sy,
    median_exposure_reference = ref,
    difference_from_median_raw = DB$estimate,
    difference_lower_raw = DB$lower, difference_upper_raw = DB$upper,
    difference_simultaneous_lower_raw = DB$sim_lower,
    difference_simultaneous_upper_raw = DB$sim_upper,
    difference_from_median_outcome_SD = DB$estimate / sy,
    difference_lower_outcome_SD = DB$lower / sy,
    difference_upper_outcome_SD = DB$upper / sy,
    exposure_mean = prep$exposure_spec$center,
    exposure_sd = prep$exposure_spec$scale, outcome_mean = mu, outcome_sd = sy,
    support_lower_probability = support[1L], support_upper_probability = support[2L],
    support_lower_raw = lim[1L], support_upper_raw = lim[2L],
    n_children = obj$n_children, n_feet = obj$n_feet,
    standardisation_population = "model_complete_case_feet_equal_weight",
    band_scope = "one_curve_finite_grid",
    stringsAsFactors = FALSE)
  ans
}

assoc_compact_fit <- function(obj) {
  list(beta = obj$beta, V = obj$V,
       formula = stats::formula(obj$fit),
       n_children = obj$n_children, n_feet = obj$n_feet,
       working_correlation = obj$fit$geese$alpha,
       convergence_code = obj$fit$geese$error)
}

assoc_fit_one <- function(d, cfg, exposure, outcome, adjust_bmi,
                         analysis_set, model_id, seed) {
  prep <- assoc_prepare(d, exposure, outcome, adjust_bmi, cfg)
  z <- prep$data
  spec <- prep$exposure_spec
  f_lin <- assoc_formula(outcome, "assoc_x1", prep$age_names, adjust_bmi)
  f_nl <- assoc_formula(outcome, spec$terms, prep$age_names, adjust_bmi)
  lin <- fit_gee(z, f_lin, cfg)
  nl <- fit_gee(z, f_nl, cfg)
  if (lin$n_feet != nrow(z) || nl$n_feet != nrow(z)) {
    stop("Association helper unexpectedly changed the complete-case sample.")
  }
  coeff <- assoc_linear_coefficients(lin, prep, exposure, outcome, model_id,
                                    analysis_set)
  curves <- assoc_curve(nl, prep, exposure, outcome, model_id, analysis_set, cfg, seed)
  tests <- list(
    assoc_wald_row(nl, match(spec$terms, names(nl$beta)), "overall_exposure", model_id),
    assoc_wald_row(nl, match(spec$terms[-1L], names(nl$beta)), "nonlinearity", model_id))

  for (kind in c("sex", "age")) {
    tests[[length(tests) + 1L]] <- tryCatch({
      fi <- assoc_formula(outcome, spec$terms, prep$age_names, adjust_bmi, kind)
      oi <- fit_gee(z, fi, cfg)
      bn <- names(oi$beta)
      is_x <- vapply(strsplit(bn, ":", fixed = TRUE),
        function(parts) any(parts %in% spec$terms), logical(1))
      is_i <- grepl(":", bn, fixed = TRUE)
      target <- if (kind == "sex") grepl("sex", bn, fixed = TRUE) else
        grepl("assoc_age", bn, fixed = TRUE)
      assoc_wald_row(oi, which(is_x & is_i & target),
                     paste0(kind, "_by_exposure"), model_id)
    }, error = function(e) {
      data.frame(model_id = model_id, test = paste0(kind, "_by_exposure"),
        statistic = NA_real_, df = NA_real_, p = NA_real_,
        status = "failed", note = conditionMessage(e),
        n_children = nl$n_children, n_feet = nl$n_feet,
        stringsAsFactors = FALSE)
    })
  }
  tests <- bind_rows_base(tests)
  tests$analysis_set <- analysis_set
  tests$exposure <- exposure
  tests$outcome <- outcome
  info <- data.frame(model_id = model_id, analysis_set = analysis_set,
    exposure = exposure, outcome = outcome,
    linear_formula = paste(deparse(f_lin, width.cutoff = 500), collapse = " "),
    nonlinear_formula = paste(deparse(f_nl, width.cutoff = 500), collapse = " "),
    n_input_children = prep$n_input_children, n_input_feet = prep$n_input_feet,
    n_children = nl$n_children, n_feet = nl$n_feet,
    n_excluded_feet = prep$n_excluded_feet,
    n_single_foot_children = prep$n_single_foot_children,
    n_linear_coefficients = length(lin$beta), n_nonlinear_coefficients = length(nl$beta),
    exposure_mean = spec$center, exposure_sd = spec$scale,
    outcome_mean = mean(z[[outcome]]), outcome_sd = stats::sd(z[[outcome]]),
    knot1_raw = spec$knots_raw[1L], knot2_raw = spec$knots_raw[2L],
    knot3_raw = spec$knots_raw[3L], knot4_raw = spec$knots_raw[4L],
    working_correlation_linear = as.numeric(lin$fit$geese$alpha)[1L],
    working_correlation_nonlinear = as.numeric(nl$fit$geese$alpha)[1L],
    AV_definition = if (outcome == "AV" || exposure == "AV")
      "raw_arch_volume_mm3_FL_covariate_adjusted" else if (outcome == "AV_FL3")
      "dimensionless_AV_divided_by_FL_cubed_sensitivity" else "not_applicable",
    stringsAsFactors = FALSE)
  model <- list(linear = assoc_compact_fit(lin), nonlinear = assoc_compact_fit(nl),
    exposure_spec = spec, age_knots = assoc_setting(cfg, "age_knots", c(6, 9, 12, 15)),
    age_boundary = assoc_setting(cfg, "age_boundary", c(3, 19)),
    complete_case_row_id = if ("row_id" %in% names(z)) z$row_id else NULL,
    standardisation = list(exposure_mean = spec$center, exposure_sd = spec$scale,
      outcome_mean = mean(z[[outcome]]), outcome_sd = stats::sd(z[[outcome]])))
  list(coefficients = coeff, curves = curves, tests = tests, info = info, model = model)
}

run_associations <- function(d, cfg, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  plan <- expand.grid(exposure = c("BMI", "hallux"),
    outcome = c("AI", "AS_ratio", "AV"), stringsAsFactors = FALSE)
  plan$adjust_bmi <- plan$exposure == "hallux"
  plan$analysis_set <- "primary"

  if (isTRUE(assoc_setting(cfg, "association_av_ratio_sensitivity", TRUE))) {
    d$AV_FL3 <- ifelse(is.finite(d$FL) & d$FL > 0, d$AV / d$FL^3, NA_real_)
    plan <- rbind(plan, data.frame(exposure = c("BMI", "hallux"),
      outcome = "AV_FL3", adjust_bmi = c(FALSE, TRUE),
      analysis_set = "AV_FL3_sensitivity", stringsAsFactors = FALSE))
  }

  if (isTRUE(assoc_setting(cfg, "association_reverse_supplementary", FALSE))) {
    plan <- rbind(plan, data.frame(exposure = c("AI", "AS_ratio", "AV"),
      outcome = "hallux", adjust_bmi = TRUE,
      analysis_set = "reverse_direction_supplementary", stringsAsFactors = FALSE))
  }
  plan$model_id <- paste(plan$analysis_set, plan$exposure, "to", plan$outcome, sep = "_")
  answers <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    p <- plan[i, ]
    if (isTRUE(cfg$verbose)) message("  Association: ", p$model_id)
    answers[[i]] <- assoc_fit_one(d, cfg, p$exposure, p$outcome, p$adjust_bmi,
      p$analysis_set, p$model_id,
      seed = as.integer(assoc_setting(cfg, "seed", 20260928L)) + 100L * i)
  }
  coeff <- bind_rows_base(lapply(answers, `[[`, "coefficients"))
  curves <- bind_rows_base(lapply(answers, `[[`, "curves"))
  tests <- bind_rows_base(lapply(answers, `[[`, "tests"))
  info <- bind_rows_base(lapply(answers, `[[`, "info"))

  coeff$q_BH <- NA_real_
  for (g in unique(coeff$analysis_set)) for (s in unique(coeff$scale)) {
    k <- coeff$analysis_set == g & coeff$scale == s
    coeff$q_BH[k] <- stats::p.adjust(coeff$p[k], method = "BH")
  }
  tests$q_BH <- NA_real_
  for (g in unique(tests$analysis_set)) for (t in unique(tests$test)) {
    k <- tests$analysis_set == g & tests$test == t
    tests$q_BH[k] <- stats::p.adjust(tests$p[k], method = "BH")
  }
  write_table(coeff, file.path(out, "association_linear_coefficients.csv"))
  write_table(curves, file.path(out, "association_adjusted_curves.csv"))
  write_table(tests, file.path(out, "association_robust_wald_tests.csv"))
  write_table(info, file.path(out, "association_model_audit.csv"))
  write_table(curves[curves$analysis_set == "primary", , drop = FALSE],
              file.path(out, "figure4_primary_six_curves.csv"))
  models <- lapply(answers, `[[`, "model")
  names(models) <- plan$model_id
  saveRDS(models, file.path(out, "association_models_compact.rds"), compress = "gzip")
  notes <- c(
    "Primary direction follows 09_six_curves_no_points_v2.R: BMI/hallux -> morphology.",
    "BMI -> AI/AS_ratio/raw AV: fixed age natural spline x sex + recorded side + FL.",
    "Hallux -> AI/AS_ratio/raw AV: the same covariates plus BMI.",
    "Raw AV is the measured arch volume in mm^3. It is never replaced by an age/FL residual or by AV/FL^3 in the primary models.",
    "Exposure centring/scaling is numerical parameterisation only. Raw units are restored in every raw-scale output.",
    "Outcome SD curves use (adjusted fitted raw outcome - complete-case observed outcome mean)/observed outcome SD.",
    "Fully standardised linear coefficients are raw slopes x observed exposure SD / observed outcome SD.",
    "SDs and marginal-standardisation weights are foot-level and specific to each model's complete cases.",
    "Curves use empirical covariate standardisation over the complete-case feet, preserving observed age-by-sex design moments.",
    "Curve intervals describe adjusted conditional means, not individual reference or prediction intervals.",
    "Simultaneous curve bands apply to one finite curve grid, not to all six panels jointly.",
    "The nonlinearity test jointly tests the two explicitly nonlinear restricted-cubic-spline terms with robust sandwich covariance.",
    "Age-by-exposure and sex-by-exposure are separate joint interaction tests in augmented nonlinear models.",
    "An interaction fit failure is recorded in association_robust_wald_tests.csv; it is never interpreted as nonsignificance.",
    "BH adjustments are separate for each test type and analysis set; raw p values remain available.",
    "Optional AV/FL^3 sensitivities retain FL as a covariate and have a distinct dimensionless estimand.",
    "Cross-sectional association directions describe regression roles and do not establish temporal or causal direction."
  )
  writeLines(notes, file.path(out, "association_analysis_notes.txt"), useBytes = TRUE)
  invisible(list(coefficients = coeff, curves = curves, tests = tests,
                 audit = info, models = models))
}
