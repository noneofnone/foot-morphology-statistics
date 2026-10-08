gamma_bind_rows <- function(parts) {
  parts <- Filter(function(x) is.data.frame(x) && nrow(x), parts)
  if (!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, parts); rownames(out) <- NULL; out
}

gamma_replace_av_rows <- function(old, new, metric_col = "metric") {
  if (!metric_col %in% names(old)) stop("Missing metric selector: ", metric_col)
  gamma_bind_rows(list(old[old[[metric_col]] != "AV_sizeadjusted", , drop = FALSE], new))
}

gamma_partial_gradient <- function(obj, age, FL, h = .01, beta = obj$beta) {
  if (length(age) != length(FL) || any(!is.finite(FL) | FL <= 0)) stop("Invalid reference FL.")
  lo <- pmax(min(age), age - h); hi <- pmin(max(age), age + h)
  dt <- hi - lo
  if (any(dt <= 0)) stop("Age derivative needs distinct supported ages.")
  estimate <- numeric(length(age))
  J <- matrix(0, nrow = length(age), ncol = length(beta), dimnames = list(NULL, names(beta)))
  side_gradients <- list()
  for (side in c("Left", "Right")) {

    XL <- gamma_design(obj, lo, FL, side)
    XH <- gamma_design(obj, hi, FL, side)
    muL <- exp(drop(XL %*% beta)); muH <- exp(drop(XH %*% beta))
    side_gradients[[side]] <- (muH - muL) / dt
    estimate <- estimate + side_gradients[[side]] / 2
    J <- J + (sweep(XH, 1, muH, "*") - sweep(XL, 1, muL, "*")) / dt / 2
  }
  list(estimate = estimate, J = J, left = side_gradients$Left, right = side_gradients$Right)
}

gamma_cumulative <- function(age, value) {
  z <- as.matrix(value)
  if (nrow(z) != length(age)) stop("Integration grid mismatch.")
  increments <- (z[-1, , drop = FALSE] + z[-nrow(z), , drop = FALSE]) / 2 * diff(age)
  out <- rbind(rep(0, ncol(z)), apply(increments, 2, cumsum))
  if (is.null(dim(value))) drop(out) else out
}

gamma_sustained_onset <- function(z, min_tail = 1) {
  z <- z[order(z$age), ]
  ok <- is.finite(z$sim_lower) & is.finite(z$sim_upper) & z$sim_lower <= 0 & z$sim_upper >= 0
  sustained <- rev(cumprod(rev(as.integer(ok)))) == 1L
  i <- which(sustained & z$age <= max(z$age) - min_tail)
  data.frame(metric = "AV_sizeadjusted", sex = z$sex[1],
    onset = if (length(i)) z$age[i[1]] else NA_real_,
    status = if (length(i)) "sustained_nondetectable_change" else "not_identified_in_observed_range",
    tail_years = if (length(i)) max(z$age) - z$age[i[1]] else NA_real_,
    equivalence_margin = NA_real_, equivalence_onset = NA_real_, equivalence_status = "margin_not_prespecified")
}

gamma_prepare_data <- function(d, cfg) {
  needed <- c("row_id", "id", "sex", "age", "side_code", "FL", "AV")
  if (!is.data.frame(d) || length(setdiff(needed, names(d)))) stop("Missing Gamma AV input columns")
  if (!all(vapply(d[c("age", "FL", "AV")], is.numeric, logical(1)))) stop("Non-numeric Gamma AV measurements")
  if (anyNA(d$id) || any(!nzchar(trimws(as.character(d$id)))) || any(!is.finite(d$age)) ||
      any(!as.character(d$sex) %in% c("Boys", "Girls")) ||
      any(!as.character(d$side_code) %in% c("1", "2"))) stop("Invalid Gamma AV participant, age, sex or side")
  d$id <- as.character(d$id)
  d$sex <- factor(d$sex, levels = c("Boys", "Girls"))
  d$side <- factor(c("Left", "Right")[match(as.character(d$side_code), c("1", "2"))],
                   levels = c("Left", "Right"))
  if (anyDuplicated(d[c("id", "side")]) || anyNA(d$row_id) || anyDuplicated(d$row_id))
    stop("Duplicate or missing Gamma AV participant/side or row identifiers")
  if (any(d$age < cfg$age_boundary[1] | d$age > cfg$age_boundary[2])) stop("Age outside Gamma AV spline boundaries")
  coherent <- vapply(split(seq_len(nrow(d)), d$id), function(ii) {
    length(unique(d$sex[ii])) == 1L && diff(range(d$age[ii])) < 1e-8
  }, logical(1))
  if (!all(coherent)) stop("Within-participant age or sex disagreement")
  d <- d[order(d$id, d$side), , drop = FALSE]
  rownames(d) <- NULL
  d
}

gamma_fit <- function(d, sex, cfg) {
  z <- d[d$sex == sex, , drop = FALSE]
  n_input <- length(unique(z$id))
  valid <- is.finite(z$AV) & z$AV > 0 & is.finite(z$FL) & z$FL > 0
  counts <- tapply(valid, z$id, sum)
  z <- z[z$id %in% names(counts)[counts == 2L], , drop = FALSE]
  n_pairs <- length(unique(z$id))
  if (n_pairs < cfg$min_clusters) stop("Too few complete AV pairs: ", sex)
  B <- splines::ns(z$age, knots = cfg$age_knots, Boundary.knots = cfg$age_boundary)
  z[paste0("a", seq_len(ncol(B)))] <- B
  center <- mean(log(z$FL))
  z$lc <- log(z$FL) - center
  z$y <- z$AV
  basis <- paste0("a", seq_len(ncol(B)), collapse = " + ")
  f <- stats::as.formula(paste0("y ~ side * (", basis, ") + lc + lc: (", basis, ")"))
  fit <- stats::glm(f, data = z, family = stats::Gamma(link = "log"),
    contrasts = list(side = stats::contr.treatment(c("Left", "Right"), base = 1L)),
    control = stats::glm.control(maxit = 100), x = TRUE, y = TRUE)
  beta <- stats::coef(fit)
  if (!isTRUE(fit$converged) || fit$rank != length(beta) || any(!is.finite(beta)))
    stop("Gamma AV model convergence or rank failure: ", sex)
  V <- sandwich::vcovCL(fit, cluster = z$id, type = "HC1", cadjust = TRUE)
  V <- V[names(beta), names(beta), drop = FALSE]
  V <- (V + t(V)) / 2
  if (any(!is.finite(V))) stop("Non-finite Gamma AV covariance")
  eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < -1e-8 * max(abs(eig))) stop("Indefinite Gamma AV covariance")
  list(fit = fit, beta = beta, V = V, metric = "AV", sex = sex, center = center,
    age_range = range(z$age), n_pairs = n_pairs, n_input = n_input, paired_data = z,
    age_knots = cfg$age_knots, age_boundary = cfg$age_boundary)
}

gamma_design <- function(obj, age, FL, side) {
  nd <- data.frame(age = age, side = factor(rep(side, length(age)), levels = c("Left", "Right")))
  nd[paste0("a", seq_len(length(obj$age_knots) + 1L))] <- splines::ns(age,
    knots = obj$age_knots, Boundary.knots = obj$age_boundary)
  nd$lc <- log(FL) - obj$center
  X <- stats::model.matrix(stats::delete.response(stats::terms(obj$fit)), nd,
                           contrasts.arg = obj$fit$contrasts, xlev = obj$fit$xlevels)
  X[, names(obj$beta), drop = FALSE]
}

gamma_reference_fl <- function(reference, sex, age) {
  ref <- reference$gradients
  support <- reference$FL_support
  if (!all(c("sex", "age", "local_FL50") %in% names(ref)) ||
      !all(c("sex", "age", "FL_low", "FL_high") %in% names(support))) stop("Missing AV reference FL or support")
  if (anyDuplicated(ref[c("sex", "age")]) || anyDuplicated(support[c("sex", "age")]))
    stop("Duplicate AV reference ages")
  rr <- ref[ref$sex == sex, , drop = FALSE]
  ss <- support[support$sex == sex, , drop = FALSE]
  if (nrow(rr) < 2L || nrow(ss) < 2L) stop("Insufficient AV reference FL support: ", sex)
  value <- stats::approx(rr$age, rr$local_FL50, xout = age, rule = 1)$y
  low <- stats::approx(ss$age, ss$FL_low, xout = age, rule = 1)$y
  high <- stats::approx(ss$age, ss$FL_high, xout = age, rule = 1)$y
  ok <- is.finite(value) & value > 0 & is.finite(low) & is.finite(high) & value >= low & value <= high
  if (!all(ok)) stop("Gamma AV reference FL is outside age/FL support")
  data.frame(sex = sex, age = age, reference_FL = value, FL_low = low, FL_high = high,
             in_support = ok, stringsAsFactors = FALSE)
}

gamma_pair_predictions <- function(obj, reference, cfg) {
  age <- cfg$age_grid[cfg$age_grid >= obj$age_range[1] & cfg$age_grid <= obj$age_range[2]]
  if (length(age) < 3L) stop("Insufficient supported Gamma AV ages")
  rr <- gamma_reference_fl(reference, obj$sex, age)
  XL <- gamma_design(obj, age, rr$reference_FL, "Left")
  XR <- gamma_design(obj, age, rr$reference_FL, "Right")
  muL <- exp(drop(XL %*% obj$beta))
  muR <- exp(drop(XR %*% obj$beta))
  J <- sweep(XR, 1, muR, "*") - sweep(XL, 1, muL, "*")
  variance <- rowSums((J %*% obj$V) * J)
  if (any(variance < -1e-10 * max(1, max(abs(variance))))) stop("Negative bilateral AV contrast variance")
  se <- sqrt(pmax(variance, 0))
  difference <- muR - muL
  if (any(!is.finite(c(muL, muR, difference, se)))) stop("Non-finite Gamma AV predictions")
  list(reference = rr,
    trajectories = rbind(
      data.frame(metric = "AV", sex = obj$sex, side = "Left", age = age, mean = muL, reference_FL = rr$reference_FL),
      data.frame(metric = "AV", sex = obj$sex, side = "Right", age = age, mean = muR, reference_FL = rr$reference_FL)),
    contrasts = data.frame(metric = "AV", sex = obj$sex, age = age, estimate = difference, se = se,
      lower = difference - qnorm(.975) * se, upper = difference + qnorm(.975) * se,
      reference_FL = rr$reference_FL, n_pairs = obj$n_pairs, interval = "pointwise 95% CI"))
}

gamma_shared_models <- function(d, reference, cfg) {
  d <- gamma_prepare_data(d, cfg)
  models <- predictions <- audits <- coefficients <- list()
  for (sex in c("Boys", "Girls")) {
    key <- paste("AV", sex, sep = "|")
    obj <- gamma_fit(d, sex, cfg)
    models[[key]] <- obj
    predictions[[key]] <- gamma_pair_predictions(obj, reference, cfg)
    audits[[key]] <- data.frame(metric = "AV", sex = sex, n_input_children = obj$n_input,
      n_complete_pairs = obj$n_pairs, n_feet = 2L * obj$n_pairs,
      family = obj$fit$family$family, link = obj$fit$family$link,
      rank = obj$fit$rank, n_coefficients = length(obj$beta), converged = obj$fit$converged,
      age_min = obj$age_range[1], age_max = obj$age_range[2], log_FL_center = obj$center,
      covariance = "participant-clustered HC1 with cluster correction",
      formula = paste(deparse(stats::formula(obj$fit)), collapse = " "), stringsAsFactors = FALSE)
    coefficients[[key]] <- data.frame(metric = "AV", sex = sex, term = names(obj$beta),
      estimate = unname(obj$beta), robust_se = sqrt(diag(obj$V)), stringsAsFactors = FALSE)
  }
  list(models = models, support = reference$FL_support,
    reference_path = gamma_bind_rows(lapply(predictions, `[[`, "reference")),
    trajectories = gamma_bind_rows(lapply(predictions, `[[`, "trajectories")),
    contrasts = gamma_bind_rows(lapply(predictions, `[[`, "contrasts")),
    model_audit = gamma_bind_rows(audits), coefficient_table = gamma_bind_rows(coefficients))
}

run_av_gamma <- function(d, reference, cfg, out = NULL) {
  n_sim <- cfg$av_age_n_sim %||% 5000L
  seed <- cfg$av_age_seed %||% 20261001L
  h <- cfg$av_age_derivative_h %||% .01
  min_tail <- cfg$min_tail_years %||% 1
  if (length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 1000 || n_sim != as.integer(n_sim))
    stop("av_age_n_sim must be an integer >= 1000")
  if (length(seed) != 1L || !is.finite(seed) || length(h) != 1L || !is.finite(h) || h <= 0 || h > .05 ||
      length(min_tail) != 1L || !is.finite(min_tail) || min_tail < 1) stop("Invalid Gamma AV inference settings")
  if (length(cfg$age_grid) < 3L || any(!is.finite(cfg$age_grid)) || any(diff(cfg$age_grid) <= 0))
    stop("Gamma AV age grid must be finite and strictly increasing")
  shared <- gamma_shared_models(d, reference, cfg)
  have_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (have_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({ if (have_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
  set.seed(seed)
  calc <- list(); max_by_sex <- matrix(0, n_sim, 2)
  sex_levels <- c("Boys", "Girls")
  for (sx in sex_levels) {
    obj <- shared$models[[paste("AV", sx, sep = "|")]]
    rr <- shared$reference_path[shared$reference_path$sex == sx, ]

    age <- round(rr$age, 8)
    g <- gamma_partial_gradient(obj, age, rr$reference_FL, h)
    variance <- rowSums((g$J %*% obj$V) * g$J)
    if (any(variance < -1e-8)) stop("Negative gradient variance.")
    se <- sqrt(pmax(variance, 0))
    if (any(!is.finite(se) | se <= 0)) stop("Nonpositive gradient SE.")
    eig <- eigen(obj$V, symmetric = TRUE)
    if (min(eig$values) < -1e-8 * max(abs(eig$values))) stop("Invalid covariance eigenvalues.")
    perturb <- eig$vectors %*% (sqrt(pmax(eig$values, 0)) * matrix(rnorm(length(obj$beta) * n_sim), nrow = length(obj$beta)))
    deviations <- g$J %*% perturb
    max_by_sex[, match(sx, sex_levels)] <- apply(abs(deviations / se), 2, max)

    gradient_draws <- deviations + g$estimate
    integral <- gamma_cumulative(age, g$estimate)
    integral_draws <- gamma_cumulative(age, gradient_draws)
    span <- tail(integral, 1)
    spans <- integral_draws[nrow(integral_draws), ]
    span_ci <- as.numeric(stats::quantile(spans, c(.025, .975)))
    tol <- 1e-6 * max(abs(integral), 1e-6)
    eligible <- is.finite(spans) & abs(spans) > tol & sign(spans) == sign(span)
    reliable <- is.finite(span) && abs(span) > tol && prod(span_ci) > 0 && mean(eligible) >= .95
    p <- lo <- hi <- rep(NA_real_, length(age))
    if (reliable) {
      p <- 100 * integral / span
      draws <- 100 * sweep(integral_draws[, eligible, drop = FALSE], 2, spans[eligible], "/")
      lo <- apply(draws, 1, stats::quantile, .025)
      hi <- apply(draws, 1, stats::quantile, .975)
    }
    calc[[sx]] <- list(obj = obj, rr = rr, gradient = g, se = se,
      gradient_table = data.frame(metric = "AV_sizeadjusted", sex = sx, age = age,
        estimate = g$estimate, se = se, lower = g$estimate - stats::qnorm(.975) * se,
        upper = g$estimate + stats::qnorm(.975) * se, units = "mm^3/year", unit = "mm^3/year",
        local_FL50 = rr$reference_FL, local_FL_in_support = rr$in_support,
        band_family = "all supported ages and both sexes within AV", estimator = "Shared S5 Gamma mean model; cluster-robust delta method"),
      progress = data.frame(metric = "AV_sizeadjusted", sex = sx, age = age, estimate = p,
        lower = lo, upper = hi, endpoint_span = span, span_lower = span_ci[1], span_upper = span_ci[2],
        cumulative_response_change = integral, progress_defined = reliable, progress_reliable = reliable,
        display_eligible = reliable, nonmonotone = if (reliable) any(diff(p) < -1e-5) else NA,
        valid_draw_fraction = mean(eligible), interpretation = "Endpoint-normalised integral of partial mean-AV age gradient; not change along the moving-FL path"),
      heat = data.frame(metric = "AV_sizeadjusted", sex = sx, age = age,
        normalized_gradient = if (reliable) 100 * g$estimate / span else NA_real_,
        definition = "percentage of integrated partial mean-AV age change per year", raw_gradient = g$estimate))
  }
  critical <- unname(stats::quantile(apply(max_by_sex, 1, max), .95))
  gradients <- progress <- heat <- stability <- histogram <- list()
  for (sx in sex_levels) {
    z <- calc[[sx]]
    gt <- z$gradient_table
    gt$sim_lower <- gt$estimate - critical * gt$se
    gt$sim_upper <- gt$estimate + critical * gt$se
    gt$sim_critical <- critical
    gt$p <- 2 * stats::pnorm(-abs(gt$estimate / gt$se))
    gradients[[sx]] <- gt; progress[[sx]] <- z$progress
    heat[[sx]] <- transform(z$heat, raw_sim_lower = gt$sim_lower, raw_sim_upper = gt$sim_upper)
    stability[[sx]] <- gamma_sustained_onset(gt, min_tail)
    onset <- stability[[sx]]$onset
    if (is.finite(onset)) {
      dd <- z$obj$paired_data
      dd <- dd[dd$age >= onset & dd$age <= max(gt$age), ]
      if (!nrow(dd)) next
      ref_FL <- stats::median(dd$FL)
      sp <- shared$support[shared$support$sex == sx, ]
      low <- stats::approx(sp$age, sp$FL_low, xout = dd$age, rule = 1)$y
      high <- stats::approx(sp$age, sp$FL_high, xout = dd$age, rule = 1)$y
      keep <- is.finite(low) & is.finite(high) & ref_FL >= low & ref_FL <= high
      dd <- dd[keep, , drop = FALSE]
      if (!nrow(dd)) next
      adjustment <- numeric(nrow(dd))
      for (side in c("Left", "Right")) {
        i <- which(as.character(dd$side) == side)
        Xo <- gamma_design(z$obj, dd$age[i], dd$FL[i], side)
        Xr <- gamma_design(z$obj, dd$age[i], rep(ref_FL, length(i)), side)
        adjustment[i] <- exp(drop((Xr - Xo) %*% z$obj$beta))
      }
      histogram[[sx]] <- data.frame(row_id = dd$row_id, id = dd$id, sex = sx, age = dd$age,
        selection_metric = "AV_sizeadjusted", distribution_variable = "AV_at_stable_window_reference_FL",
        value = dd$AV * adjustment, onset = onset, reference_FL = ref_FL,
        interpretation = "Observed AV times shared Gamma model mean ratio at stable-window median FL versus observed FL; mm^3")
    }
  }
  result <- list(gradients = gamma_bind_rows(gradients), progress = gamma_bind_rows(progress),
    heatmap = gamma_bind_rows(heat), stability = gamma_bind_rows(stability),
    histogram_data = gamma_bind_rows(histogram),
    model_audit = shared$model_audit[shared$model_audit$metric == "AV", ],
    coefficients = shared$coefficient_table[shared$coefficient_table$metric == "AV", ],
    reference_path = shared$reference_path,
    inference = data.frame(n_sim = n_sim, seed = seed, derivative_h = h, min_tail_years = min_tail,
      simultaneous_critical = critical, mean_side_weights = "Left=0.5; Right=0.5",
      inference = "Delta method; child-cluster HC1 covariance; joint Gaussian perturbations; conditional on reference FL"),
    models = shared$models, trajectories = shared$trajectories, contrasts = shared$contrasts)
  result$progress_qc <- gamma_bind_rows(lapply(split(result$progress, result$progress$sex), function(z) {
    data.frame(sex = z$sex[1], endpoint_span = z$endpoint_span[1],
      endpoint_lower = z$span_lower[1], endpoint_upper = z$span_upper[1],
      endpoint_CI_contains_zero = z$span_lower[1] <= 0 & z$span_upper[1] >= 0,
      valid_draw_fraction = z$valid_draw_fraction[1], progress_reliable = z$progress_reliable[1],
      nonmonotone = z$nonmonotone[1], primary_gradient_unit = "mm^3/year")
  }))
  if (!is.null(out)) {
    directory <- file.path(out, "results", "av_gamma")
    for (nm in c("gradients", "progress", "heatmap", "stability", "model_audit", "coefficients",
                 "reference_path", "inference", "trajectories", "contrasts", "progress_qc"))
      write_table(result[[nm]], file.path(directory, paste0("av_gamma_", nm, ".csv")))
    atomic_rds(result, file.path(out, "models", "av_gamma_analysis.rds"))
  }
  result
}

combine_av_models <- function(reference, gamma, cfg, out) {
  result <- reference
  result$log_scale <- reference[c("gradients", "progress", "progress_qc", "estimand")]
  result$gamma <- gamma
  result$gradients <- gamma$gradients
  result$progress <- gamma$progress
  result$progress_qc <- gamma$progress_qc
  result$heatmap <- gamma$heatmap
  result$estimand <- "Fixed-FL partial age derivative of the equally weighted bilateral Gamma arithmetic mean; mm^3/year"
  for (nm in c("gradients", "progress", "progress_qc")) {
    av_write(reference[[nm]], out, paste0("av_log_scale_", nm, ".csv"))
    av_write(result[[nm]], out, paste0("av_", nm, ".csv"))
  }
  atomic_rds(reference, file.path(out, "models", "av_reference_analysis.rds"))
  atomic_rds(result, file.path(out, "models", "av_analysis.rds"))
  result
}
