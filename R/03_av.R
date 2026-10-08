av_default <- function(x, fallback) if (is.null(x)) fallback else x
av_bind <- function(x) {
  x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  rownames_out <- do.call(rbind, x); rownames(rownames_out) <- NULL; rownames_out
}
av_write <- function(x, out, name) {
  dir.create(file.path(out, "results", "av"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file.path(out, "results", "av", name), row.names = FALSE, na = "")
}
av_hash <- function(x) {
  path <- tempfile("av_hash_"); on.exit(unlink(path)); saveRDS(x, path, compress = FALSE)
  unname(tools::md5sum(path))
}
av_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp"); saveRDS(x, tmp, compress = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Cannot write AV checkpoint: ", path)
}
av_interp <- function(x, y, target) {
  if (length(unique(x)) < 2L) return(rep(y[1L], length(target)))
  stats::approx(x, y, xout = target, rule = 2, ties = mean)$y
}
av_resample <- function(d, seed) {
  set.seed(seed)
  child <- d[!duplicated(d$id), c("id", "sex", "age_group"), drop = FALSE]
  strata <- split(seq_len(nrow(child)), interaction(child$sex, child$age_group, drop = TRUE))
  selected <- unlist(lapply(strata, function(ii) ii[sample.int(length(ii), length(ii), replace = TRUE)]), use.names = FALSE)
  rows <- split(seq_len(nrow(d)), as.character(d$id))
  parts <- lapply(seq_along(selected), function(j) {
    z <- d[rows[[as.character(child$id[selected[j]])]], , drop = FALSE]
    z$id <- paste0("av_boot_", j); z
  })
  av_bind(parts)
}
av_folds <- function(d, k, seed) {
  set.seed(seed)
  child <- d[!duplicated(d$id), c("id", "sex", "age_group"), drop = FALSE]
  child$fold <- NA_integer_
  strata <- split(seq_len(nrow(child)), interaction(child$sex, child$age_group, drop = TRUE))
  for (ii in strata) child$fold[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  child$fold[match(d$id, child$id)]
}

av_fit <- function(z, cfg) {
  if (length(unique(z$id)) < av_default(cfg$min_clusters, 30L)) stop("Insufficient AV participant clusters")
  if (length(unique(z$age)) < 6L) stop("AV requires at least six distinct ages")
  z <- droplevels(z)
  center <- mean(log(z$FL))
  z$logAV <- log(z$AV); z$lc <- log(z$FL) - center
  sides <- levels(z$side)
  if (is.null(sides)) { z$side <- factor(z$side); sides <- levels(z$side) }
  ka <- min(av_default(cfg$av_age_k, 6L), length(unique(z$age)) - 1L)
  kb <- min(av_default(cfg$av_beta_k, 5L), length(unique(z$age)) - 1L)

  s <- mgcv::s
  f <- stats::as.formula(paste0("logAV ~ s(age,k=", ka, ",bs='cr') + lc + s(age,by=lc,k=", kb,
                               ",bs='cr')", if (length(sides) > 1L) " + side" else ""))
  fit <- mgcv::gam(f, data = z, method = "REML", select = TRUE)
  if (any(!is.finite(stats::coef(fit))) || !isTRUE(fit$converged)) stop("AV GAM did not converge")
  if (!is.null(fit$outer.info$conv) && !identical(fit$outer.info$conv, "full convergence"))
    stop("AV smoothing-parameter optimisation failed: ", fit$outer.info$conv)
  obj <- list(fit = fit, center = center, sides = sides, sex = as.character(z$sex[1L]))

  pooled_mu <- av_predict(obj, z$age, z$FL)
  residual <- z$logAV - pooled_mu
  probs <- sort(unique(c(.001, .005, seq(.01, .99, .01), .995, .999)))
  groups <- split(seq_len(nrow(z)), z$age_group)
  exact <- lapply(groups, function(ii) {
    data.frame(age = mean(z$age[ii]), age_group = as.character(z$age_group[ii[1L]]),
               n_children = length(unique(z$id[ii])), n_feet = length(ii),
               FL_low = as.numeric(stats::quantile(z$FL[ii], .025)),
               FL50 = stats::median(z$FL[ii]),
               FL_high = as.numeric(stats::quantile(z$FL[ii], .975)))
  })
  support <- av_bind(exact)
  ord <- order(support$age); support <- support[ord, , drop = FALSE]; groups <- groups[ord]
  qmat <- t(vapply(groups, function(ii) as.numeric(stats::quantile(residual[ii], probs, type = 8)), numeric(length(probs))))
  if (nrow(support) < 4L) stop("AV residual distribution requires at least four age groups")
  dfq <- min(av_default(cfg$av_quantile_df, 4), nrow(support) - 1)
  qmodels <- lapply(seq_along(probs), function(j) {
    stats::smooth.spline(support$age, qmat[, j], w = support$n_children, df = dfq)
  })

  flmodel <- stats::smooth.spline(support$age, log(support$FL50), w = support$n_children,
                                df = min(av_default(cfg$av_fl_path_df, 6), nrow(support) - 1))
  obj$probs <- probs; obj$qmodels <- qmodels; obj$support <- support
  obj$q_exact <- qmat; obj$flmodel <- flmodel
  obj$n_children <- length(unique(z$id)); obj$n_feet <- nrow(z)
  obj$age_range <- range(z$age)
  obj$residuals <- residual
  obj
}

av_predict <- function(obj, age, FL) {
  n <- max(length(age), length(FL)); age <- rep_len(age, n); FL <- rep_len(FL, n)
  vals <- lapply(obj$sides, function(ss) {
    nd <- data.frame(age = age, lc = log(FL) - obj$center,
                     side = factor(rep(ss, n), levels = obj$sides))
    as.numeric(stats::predict(obj$fit, nd, type = "response"))
  })
  Reduce(`+`, vals) / length(vals)
}
av_quantiles <- function(obj, age) {
  raw <- vapply(obj$qmodels, function(f) as.numeric(stats::predict(f, x = age)$y), numeric(length(age)))
  if (is.null(dim(raw))) raw <- matrix(raw, nrow = length(age))
  crossed <- if (ncol(raw) > 1L) rowSums(raw[, -1, drop = FALSE] < raw[, -ncol(raw), drop = FALSE]) > 0L else rep(FALSE, nrow(raw))
  ordered <- t(apply(raw, 1, sort))
  if (length(age) == 1L) ordered <- matrix(ordered, nrow = 1L)
  list(q = ordered, crossed = crossed, max_correction = apply(abs(ordered - raw), 1, max))
}
av_support_at <- function(obj, age) {
  s <- obj$support
  data.frame(FL_low = av_interp(s$age, s$FL_low, age),
             FL_high = av_interp(s$age, s$FL_high, age),
             age_supported = age >= obj$age_range[1] & age <= obj$age_range[2])
}
av_score <- function(z, obj, source, fold = NA_integer_) {
  mu <- av_predict(obj, z$age, z$FL)

  au <- sort(unique(z$age)); qr <- av_quantiles(obj, au); idx <- match(z$age, au)
  q <- qr$q[idx, , drop = FALSE]; e <- log(z$AV) - mu
  p <- vapply(seq_len(nrow(z)), function(i) {
    stats::approx(q[i, ], obj$probs, xout = e[i], rule = 2, ties = mean)$y
  }, numeric(1))
  support <- av_support_at(obj, z$age)
  q10 <- q[, which.min(abs(obj$probs - .1))]; q50 <- q[, which.min(abs(obj$probs - .5))]
  q90 <- q[, which.min(abs(obj$probs - .9))]
  ans <- data.frame(row_id = z$row_id, id = as.character(z$id), metric = "Conditional_AV",
             sex = as.character(z$sex), age = z$age, age_group = z$age_group,
             side = as.character(z$side), percentile = p, percentile_percent = 100 * p, source = source, fold = fold,
             AV_raw = z$AV, FL = z$FL, log_location = mu, geometric_location_AV = exp(mu),
             AV_conditional_P10 = exp(mu + q10), AV_conditional_P50 = exp(mu + q50),
             AV_conditional_P90 = exp(mu + q90), AV_log_residual = e,
             AV_relative_deviation_pct = 100 * expm1(e),
             in_support = support$age_supported & z$FL >= support$FL_low & z$FL <= support$FL_high,
             FL_support_low = support$FL_low, FL_support_high = support$FL_high,
             tail_flag = ifelse(e < q[, 1], "below_probability_grid", ifelse(e > q[, ncol(q)], "above_probability_grid", "inside_probability_grid")),
             quantiles_reordered = qr$crossed[idx], quantile_max_log_correction = qr$max_correction[idx],
             reference_target = "Age-sex-FL reference; location standardised to equal recorded sides",
             stringsAsFactors = FALSE)
  for (pp in c(.03, .1, .25, .5, .75, .9, .97)) {
    jj <- which.min(abs(obj$probs - pp))
    ans[[paste0("AV_P", round(100 * pp))]] <- exp(mu + q[, jj])
  }

  position <- (z$FL - support$FL_low) / pmax(support$FL_high - support$FL_low, 1e-8)
  ans$FL_band <- as.character(cut(position, c(-Inf, 1/3, 2/3, Inf), labels = c("smaller_FL", "middle_FL", "larger_FL")))
  ans$percentile_extrapolated <- ans$percentile
  ans$percentile[!ans$in_support] <- NA_real_
  ans$percentile_percent[!ans$in_support] <- NA_real_
  ans
}

av_gradient <- function(obj, age, cfg) {
  h <- av_default(cfg$av_derivative_h, .01)
  lo <- pmax(min(age), age - h); hi <- pmin(max(age), age + h)
  FL50 <- exp(stats::predict(obj$flmodel, age)$y)

  p_lo <- av_predict(obj, lo, FL50); p_hi <- av_predict(obj, hi, FL50)
  g <- (p_hi - p_lo) / (hi - lo)

  fl_lo <- exp(stats::predict(obj$flmodel, lo)$y); fl_hi <- exp(stats::predict(obj$flmodel, hi)$y)
  total <- (av_predict(obj, hi, fl_hi) - av_predict(obj, lo, fl_lo)) / (hi - lo)
  logh <- 1e-4
  beta <- (av_predict(obj, age, FL50 * exp(logh)) - av_predict(obj, age, FL50 * exp(-logh))) / (2 * logh)
  support <- av_support_at(obj, age)
  data.frame(metric = "AV_sizeadjusted", sex = obj$sex, age = age, estimate = g,
             unit = "log(AV)/year", instantaneous_percent_per_year = 100 * g,
             local_FL50 = FL50, allometric_beta = beta, total_age_gradient_logAV = total,
             foot_length_growth_component = total - g,
             local_FL_in_support = support$age_supported & FL50 >= support$FL_low & FL50 <= support$FL_high,
             log_location_at_local_FL50 = av_predict(obj, age, FL50), stringsAsFactors = FALSE)
}
av_integral <- function(age, g) c(0, cumsum(diff(age) * (head(g, -1) + tail(g, -1)) / 2))
av_boot_one <- function(task, d, cfg) {
  tryCatch({
    z <- av_resample(d, cfg$seed + 50000L + task)
    ans <- lapply(c("Boys", "Girls"), function(ss) {
      ob <- av_fit(z[as.character(z$sex) == ss, , drop = FALSE], cfg)
      sex_grid <- av_default(cfg$av_grids[[ss]], cfg$age_grid)
      gr <- av_gradient(ob, sex_grid, cfg)

      qr <- av_quantiles(ob, sex_grid)
      list(gradient = gr$estimate, integral = av_integral(gr$age, gr$estimate),
           crossing_fraction = mean(qr$crossed), beta = gr$allometric_beta)
    })
    names(ans) <- c("Boys", "Girls")
    list(ok = TRUE, values = ans, error = "")
  }, error = function(e) list(ok = FALSE, values = NULL, error = conditionMessage(e)))
}
av_bands <- function(point, draws, level = .95) {
  alpha <- (1 - level) / 2
  qs <- t(apply(draws, 1, stats::quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE))
  se <- apply(draws, 1, stats::sd, na.rm = TRUE)
  good <- is.finite(se) & se > 1e-12
  if (any(good)) {
    z <- sweep(draws[good, , drop = FALSE], 1, point[good], "-")
    z <- sweep(z, 1, se[good], "/")
    crit <- as.numeric(stats::quantile(apply(abs(z), 2, max, na.rm = TRUE), level, na.rm = TRUE))
  } else crit <- 0
  data.frame(lower = qs[, 1], upper = qs[, 2], sim_lower = point - crit * se,
             sim_upper = point + crit * se, bootstrap_se = se, simultaneous_critical = crit)
}

av_curves <- function(obj, cfg) {
  ages <- sort(unique(obj$support$age)); nodes <- av_default(cfg$av_chart_n_FL, 80L)
  display_probs <- c(.03, .1, .25, .5, .75, .9, .97)
  ans <- lapply(ages, function(a) {
    sp <- av_support_at(obj, a); fl <- seq(sp$FL_low, sp$FL_high, length.out = nodes)
    mu <- av_predict(obj, rep(a, nodes), fl); qr <- av_quantiles(obj, a)
    av_bind(lapply(display_probs, function(p) {
      j <- which.min(abs(obj$probs - p))
      data.frame(sex = obj$sex, age = a, FL = fl, probability = p,
                 percentile = paste0("P", round(100 * p)), AV_centile = exp(mu + qr$q[1, j]),
                 FL_low = sp$FL_low, FL_high = sp$FL_high, unit = "mm^3")
    }))
  })
  av_bind(ans)
}
av_coverage <- function(scores) {
  out <- list(); j <- 0L
  for (source in unique(scores$source)) for (sex in unique(scores$sex)) {
    z <- scores[scores$source == source & scores$sex == sex, , drop = FALSE]
    for (group in c("Overall", as.character(sort(unique(z$age_group))))) {
      zz <- if (group == "Overall") z else z[as.character(z$age_group) == group, , drop = FALSE]
      for (scope in c("all", "FL_support")) {
        zz2 <- if (scope == "all") zz else zz[zz$in_support, , drop = FALSE]
        if (!nrow(zz2)) next
        for (p in c(.03, .1, .25, .5, .75, .9, .97)) {
          j <- j + 1L

          out[[j]] <- data.frame(source = source, sex = sex, age_group = group, scope = scope,
             probability = p, observed_below = mean(zz2$AV_raw < zz2[[paste0("AV_P", round(100 * p))]]),
             observed_at_or_below = mean(zz2$AV_raw <= zz2[[paste0("AV_P", round(100 * p))]]),
             n_children = length(unique(zz2$id)), n_feet = nrow(zz2))
        }
      }
    }
  }
  av_bind(out)
}

av_FL_calibration <- function(scores) {
  groups <- split(seq_len(nrow(scores)), interaction(scores$source, scores$sex, scores$FL_band, drop = TRUE))
  av_bind(lapply(groups, function(ii) {
    z <- scores[ii, , drop = FALSE]
    av_bind(lapply(c(.03, .1, .25, .5, .75, .9, .97), function(p) {
      q <- z[[paste0("AV_P", round(100 * p))]]
      data.frame(source = z$source[1], sex = z$sex[1], FL_band = z$FL_band[1], probability = p,
        observed_below = mean(z$AV_raw < q), observed_at_or_below = mean(z$AV_raw <= q),
        mean_log_residual = mean(z$AV_log_residual), n_children = length(unique(z$id)), n_feet = nrow(z),
        interpretation = "FL bands are thirds of the age-sex empirical P2.5-P97.5 range, with tails included")
    }))
  }))
}

run_av <- function(d, cfg, out) {
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Install mgcv before running AV analysis")
  cfg$seed <- av_default(cfg$seed, 20260928L)
  cfg$age_grid <- av_default(cfg$age_grid, seq(3, 18, .05))
  cfg$av_n_boot <- as.integer(av_default(cfg$av_n_boot, 500L))
  cfg$workers <- as.integer(av_default(cfg$workers, 6L))
  cfg$folds <- as.integer(av_default(cfg$folds, 5L))
  cfg$min_boot_success <- av_default(cfg$min_boot_success, .95)
  cfg$av_progress_min_span <- av_default(cfg$av_progress_min_span, 1e-8)
  if (cfg$av_n_boot < 2L) stop("AV uncertainty needs at least two bootstrap resamples")
  needed <- c("row_id", "id", "sex", "age", "age_group", "side", "FL", "AV")
  if (length(setdiff(needed, names(d)))) stop("Missing AV columns: ", paste(setdiff(needed, names(d)), collapse = ", "))
  keep <- complete.cases(d[, needed]) & is.finite(d$age) & is.finite(d$FL) & d$FL > 0 & is.finite(d$AV) & d$AV > 0
  excluded <- d[!keep, needed, drop = FALSE]; z <- d[keep, , drop = FALSE]
  z$id <- as.character(z$id); z$side <- factor(z$side)
  if (!all(c("Boys", "Girls") %in% as.character(z$sex))) stop("AV module requires both sex groups")
  if (anyDuplicated(z$row_id)) stop("AV row_id must be unique")
  cfg$age_grid <- sort(unique(cfg$age_grid))
  if (length(cfg$age_grid) < 3L) stop("AV age grid needs at least three points")
  if (isTRUE(cfg$verbose)) message("AV: fitting independent varying-coefficient models and conditional references")
  models <- setNames(lapply(c("Boys", "Girls"), function(ss) av_fit(z[as.character(z$sex) == ss, , drop = FALSE], cfg)), c("Boys", "Girls"))
  cfg$av_grids <- lapply(models, function(m) {
    grid <- cfg$age_grid[cfg$age_grid >= m$age_range[1] & cfg$age_grid <= m$age_range[2]]
    if (length(grid) < 3L) stop("Insufficient supported AV ages in configured grid for ", m$sex)
    grid
  })
  score_parts <- lapply(names(models), function(ss) av_score(z[as.character(z$sex) == ss, , drop = FALSE], models[[ss]], "apparent"))
  folds <- av_folds(z, cfg$folds, cfg$seed + 741L)
  cv_ledger <- list(); q <- 0L
  if (isTRUE(cfg$verbose)) message("AV: participant-stratified full-pipeline out-of-fold calibration")
  for (f in seq_len(cfg$folds)) for (ss in names(models)) {
    tr <- z[folds != f & as.character(z$sex) == ss, , drop = FALSE]
    te <- z[folds == f & as.character(z$sex) == ss, , drop = FALSE]
    if (!nrow(te)) next
    fit <- av_fit(tr, cfg)
    sc <- av_score(te, fit, "OOF", f)
    score_parts[[length(score_parts) + 1L]] <- sc
    q <- q + 1L
    cv_ledger[[q]] <- data.frame(sex = ss, fold = f, n_train_children = length(unique(tr$id)),
      n_test_children = length(unique(te$id)), participant_overlap = length(intersect(unique(tr$id), unique(te$id))),
      log_RMSE = sqrt(mean(sc$AV_log_residual^2)), log_MAE = mean(abs(sc$AV_log_residual)),
      reordered_fraction = mean(sc$quantiles_reordered), outside_support_fraction = mean(!sc$in_support))
  }
  scores <- av_bind(score_parts)
  apparent <- scores[scores$source == "apparent", , drop = FALSE]
  oof <- scores[scores$source == "OOF", , drop = FALSE]
  if (!setequal(oof$row_id, z$row_id) || anyDuplicated(oof$row_id)) stop("AV OOF predictions are incomplete or duplicated")
  gradient_parts <- lapply(names(models), function(ss) av_gradient(models[[ss]], cfg$av_grids[[ss]], cfg))
  names(gradient_parts) <- names(models)
  if (isTRUE(cfg$verbose)) message("AV: ", cfg$av_n_boot, " stratified participant bootstraps; ", cfg$workers, " worker(s)")
  workers <- max(1L, min(cfg$workers, cfg$av_n_boot))

  env <- environment(run_av)
  exports <- ls(env, all.names = TRUE); exports <- exports[grepl("^av_", exports)]
  exports <- exports[vapply(exports, function(nm) is.function(get(nm, envir = env)), logical(1))]
  code <- lapply(exports, function(nm) {
    fun <- get(nm, envir = env)
    if (is.function(fun)) paste(deparse(body(fun)), collapse = "\n") else NULL
  }); names(code) <- exports
  code$run_av <- paste(deparse(body(run_av)), collapse = "\n")
  checksum_cfg <- cfg; checksum_cfg$workers <- NULL; checksum_cfg$av_boot_batch <- NULL
  checksum_cfg$av_grids <- cfg$av_grids
  checkpoint_key <- av_hash(list(data = z[, needed], config = checksum_cfg, code = code,
    R = R.version.string, mgcv = as.character(utils::packageVersion("mgcv"))))
  checkpoint_dir <- file.path(out, "checkpoints", "av", checkpoint_key)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  if (workers > 1L) {
    pool <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(pool), add = TRUE)
    parallel::clusterExport(pool, exports, envir = env)
    parallel::clusterEvalQ(pool, { Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1); NULL })
  }
  batch_size <- max(1L, as.integer(av_default(cfg$av_boot_batch, 16L)))
  batches <- split(seq_len(cfg$av_n_boot), ceiling(seq_len(cfg$av_n_boot) / batch_size))
  boot <- vector("list", cfg$av_n_boot)
  for (tasks in batches) {
    checkpoint_path <- file.path(checkpoint_dir, sprintf("batch_%06d_%06d.rds", min(tasks), max(tasks)))
    cached <- if (!identical(cfg$resume, FALSE) && file.exists(checkpoint_path)) tryCatch(readRDS(checkpoint_path), error = function(e) NULL) else NULL
    valid_cache <- is.list(cached) && identical(cached$key, checkpoint_key) && identical(cached$tasks, tasks) && length(cached$results) == length(tasks)
    if (valid_cache) {
      b <- cached$results
    } else {
      b <- if (workers > 1L) parallel::parLapplyLB(pool, tasks, av_boot_one, d = z, cfg = cfg) else
        lapply(tasks, av_boot_one, d = z, cfg = cfg)
      av_atomic_rds(list(key = checkpoint_key, tasks = tasks, results = b), checkpoint_path)
    }
    boot[tasks] <- b
    if (isTRUE(cfg$verbose)) message("AV bootstrap ", max(tasks), "/", cfg$av_n_boot,
            if (valid_cache) " (checkpoint reused)" else " (checkpoint saved)")
  }
  success <- vapply(boot, `[[`, logical(1), "ok")
  boot_ledger <- data.frame(replicate = seq_along(boot), success = success,
                            error = vapply(boot, `[[`, character(1), "error"))
  av_write(boot_ledger, out, "av_bootstrap_ledger.csv")
  if (mean(success) < cfg$min_boot_success || sum(success) < 2L) {
    stop("AV bootstrap success ", sum(success), "/", length(success), " is below the configured requirement; inspect av_bootstrap_ledger.csv")
  }
  boot <- boot[success]; progress_parts <- list(); progress_qc <- list()
  for (ss in names(models)) {
    g <- gradient_parts[[ss]]
    draws <- vapply(boot, function(b) b$values[[ss]]$gradient, numeric(nrow(g)))
    g <- cbind(g, av_bands(g$estimate, draws)); g$n_boot_success <- ncol(draws)
    gradient_parts[[ss]] <- g
    integral <- av_integral(g$age, g$estimate)
    iboot <- vapply(boot, function(b) b$values[[ss]]$integral, numeric(nrow(g)))
    denom <- tail(integral, 1); den_boot <- iboot[nrow(iboot), ]
    valid <- is.finite(den_boot) & abs(den_boot) > cfg$av_progress_min_span
    iband <- av_bands(integral, iboot)
    reliable <- abs(denom) > cfg$av_progress_min_span
    progress <- if (reliable) 100 * integral / denom else rep(NA_real_, length(integral))
    if (reliable && sum(valid) >= 2L) {
      pdraw <- sweep(iboot[, valid, drop = FALSE], 2, den_boot[valid], "/") * 100
      pband <- av_bands(progress, pdraw)
    } else pband <- data.frame(lower = NA_real_, upper = NA_real_, sim_lower = NA_real_, sim_upper = NA_real_,
                               bootstrap_se = NA_real_, simultaneous_critical = NA_real_)[rep(1, length(integral)), ]
    dend <- as.numeric(stats::quantile(den_boot, c(.025, .975)))
    interval_zero <- dend[1] <= 0 && dend[2] >= 0
    progress_parts[[ss]] <- data.frame(metric = "AV_sizeadjusted", sex = ss, age = g$age,
       estimate = progress, pband, cumulative_log_change = integral,
       cumulative_log_lower = iband$lower, cumulative_log_upper = iband$upper,
       accumulated_index = exp(integral), endpoint_log_span = denom,
       endpoint_span_CI_contains_zero = interval_zero,
       progress_reliable = reliable && !interval_zero,
       interpretation = "Endpoint-relative cumulative local fixed-FL age-gradient index; 0 at minimum age and 100 at maximum age",
       stringsAsFactors = FALSE)
    pnames <- c("estimate", "lower", "upper", "sim_lower", "sim_upper")
    for (pn in pnames) progress_parts[[ss]][[paste0(pn, "_diagnostic")]] <- progress_parts[[ss]][[pn]]
    if (!reliable || interval_zero) for (pn in pnames) progress_parts[[ss]][[pn]] <- NA_real_
    progress_qc[[ss]] <- data.frame(sex = ss, endpoint_log_span = denom,
       endpoint_lower = dend[1], endpoint_upper = dend[2], endpoint_CI_contains_zero = interval_zero,
       bootstrap_near_zero_denominator_fraction = mean(!valid),
       bootstrap_opposite_endpoint_sign_fraction = mean(sign(den_boot) != sign(denom)),
       nonmonotone = if (reliable) any(diff(progress) < -1e-8) else NA,
       outside_0_100 = if (reliable) any(progress < -1e-8 | progress > 100 + 1e-8) else NA,
       primary_gradient_unit = "log(AV)/year")
  }
  gradients <- av_bind(gradient_parts); progress <- av_bind(progress_parts)
  curves <- av_bind(lapply(models, av_curves, cfg = cfg))
  support <- av_bind(lapply(models, function(m) cbind(sex = m$sex, m$support)))
  quantile_ledger <- av_bind(lapply(models, function(m) {
    aa <- cfg$av_grids[[m$sex]]; qq <- av_quantiles(m, aa)
    data.frame(sex = m$sex, age = aa, crossing_before_rearrangement = qq$crossed,
      maximum_log_quantile_correction = qq$max_correction,
      correction_method = "Increasing rearrangement separately at each age")
  }))
  hist <- av_bind(lapply(c("AV_raw", "Conditional_AV", "AV_relative_deviation_pct"), function(m) {
    value <- switch(m, AV_raw = scores$AV_raw, Conditional_AV = scores$percentile_percent,
                    AV_relative_deviation_pct = scores$AV_relative_deviation_pct)
    data.frame(row_id = scores$row_id, id = scores$id, sex = scores$sex, age = scores$age,
      age_group = scores$age_group, side = scores$side, metric = m, value = value,
      source = scores$source, in_support = scores$in_support,
      unit = switch(m, AV_raw = "mm^3", Conditional_AV = "percentile", AV_relative_deviation_pct = "% deviation from geometric location"))
  }))
  diagnostics <- av_bind(lapply(models, function(m) {
    gr <- gradients[gradients$sex == m$sex, ]
    cc <- tryCatch(max(mgcv::concurvity(m$fit, full = TRUE)["estimate", ], na.rm = TRUE), error = function(e) NA_real_)
    data.frame(sex = m$sex, n_children = m$n_children, n_feet = m$n_feet,
      min_allometric_beta = min(gr$allometric_beta), max_allometric_beta = max(gr$allometric_beta),
      negative_beta_fraction = mean(gr$allometric_beta < 0), max_concurvity_estimate = cc,
      local_path_outside_support_fraction = mean(!gr$local_FL_in_support),
      bootstrap_success_fraction = mean(success),
      logFL_center = m$center,
      specification = paste(paste(deparse(stats::formula(m$fit)), collapse = " "), "; REML; select=TRUE"))
  }))
  results <- list(scores = scores, scores_apparent = apparent, scores_oof = oof,
       gradients = gradients, progress = progress, histogram_data = hist,
       models = models, conditional_centiles = curves, FL_support = support, quantile_ledger = quantile_ledger,
       calibration = av_coverage(scores), calibration_by_FL = av_FL_calibration(scores), cv_ledger = av_bind(cv_ledger),
       bootstrap_ledger = boot_ledger, diagnostics = diagnostics, progress_qc = av_bind(progress_qc),
       excluded = excluded, config = cfg, checkpoint_key = checkpoint_key,
       estimand = "Partial derivative of the sex-specific fitted log-volume location with respect to age, evaluated at local FL50, holding FL fixed in each difference")
  tables <- c("scores", "gradients", "progress", "histogram_data", "conditional_centiles", "FL_support",
              "calibration", "calibration_by_FL", "cv_ledger", "quantile_ledger", "diagnostics", "progress_qc", "excluded")
  for (nm in tables) av_write(results[[nm]], out, paste0("av_", nm, ".csv"))
  dir.create(file.path(out, "models"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(out, "models", "av_analysis.rds"), compress = "gzip")
  results
}

av_adjust_to_FL <- function(d, models, reference_FL) {
  av_bind(lapply(names(models), function(ss) {
    z <- d[as.character(d$sex) == ss & is.finite(d$AV) & d$AV > 0 & is.finite(d$FL) & d$FL > 0, , drop = FALSE]
    if (!nrow(z)) return(data.frame())
    ref <- if (length(reference_FL) == 1L) reference_FL else reference_FL[[ss]]
    if (!is.finite(ref) || ref <= 0) stop("reference_FL must be positive")
    mu_obs <- av_predict(models[[ss]], z$age, z$FL)
    mu_ref <- av_predict(models[[ss]], z$age, rep(ref, nrow(z)))
    sp <- av_support_at(models[[ss]], z$age)
    data.frame(row_id = z$row_id, id = z$id, sex = ss, age = z$age, reference_FL = ref,
       AV_at_reference_FL = z$AV * exp(mu_ref - mu_obs), unit = "mm^3",
       reference_FL_in_support = sp$age_supported & ref >= sp$FL_low & ref <= sp$FL_high)
  }))
}
