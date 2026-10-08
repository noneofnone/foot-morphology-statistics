ref_opt <- function(cfg, name, default) {
  if (is.null(cfg[[name]])) default else cfg[[name]]
}

ref_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(z) {
    for (nm in setdiff(cols, names(z))) z[[nm]] <- NA
    z[, cols, drop = FALSE]
  })
  rownames_out <- do.call(rbind, x)
  rownames(rownames_out) <- NULL
  rownames_out
}

ref_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (exists('write_table', mode = 'function')) write_table(x, path) else
    utils::write.csv(x, path, row.names = FALSE, na = '')
}

ref_formula <- function(df, response = FALSE) {
  e <- new.env(parent = asNamespace('gamlss'))
  e$cs <- gamlss::cs
  rhs <- if (df == 0) '1' else sprintf('cs(age, df = %.8g)', df)
  stats::as.formula(paste(if (response) 'y ~' else '~', rhs), env = e)
}

ref_predict <- function(ob, ages) {
  if (!length(ages) || any(!is.finite(ages))) stop('Invalid prediction ages.')
  ua <- sort(unique(as.numeric(ages)))
  warns <- character()
  pr <- withCallingHandlers(
    gamlss::predictAll(ob$fit, newdata = data.frame(age = ua),
                      type = 'response', data = ob$data),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart('muffleWarning') }
  )
  if (!all(c('mu', 'sigma', 'nu') %in% names(pr))) stop('Missing BCCG parameters.')
  vals <- lapply(pr[c('mu', 'sigma', 'nu')], function(v) {
    if (!length(v) %in% c(1L, length(ua))) stop('BCCG prediction length mismatch.')
    rep_len(as.numeric(v), length(ua))
  })
  ans <- data.frame(age = ua, M = vals[[1]], S = vals[[2]], L = vals[[3]])
  if (any(!is.finite(as.matrix(ans))) || any(ans$M <= 0 | ans$S <= 0))
    stop('Invalid BCCG parameter prediction.')
  ans <- ans[match(ages, ua), , drop = FALSE]
  rownames(ans) <- NULL
  attr(ans, 'warnings') <- unique(warns)
  ans
}

ref_quantiles <- function(pa, probs) {
  if (any(!is.finite(probs) | probs <= 0 | probs >= 1)) stop('Invalid centile probabilities.')
  ans <- vapply(probs, function(p)
    gamlss.dist::qBCCG(p, mu = pa$M, sigma = pa$S, nu = pa$L), numeric(nrow(pa)))
  ans <- matrix(ans, nrow = nrow(pa), ncol = length(probs))
  if (any(!is.finite(ans) | ans <= 0)) stop('Nonfinite/nonpositive BCCG centile.')
  if (ncol(ans) > 1L && any(ans[, -1, drop = FALSE] < ans[, -ncol(ans), drop = FALSE]))
    stop('BCCG centiles are not ordered.')
  ans
}

ref_cdf <- function(y, pa) {
  if (any(!is.finite(y) | y <= 0)) stop('BCCG scoring requires positive finite observations.')
  p <- gamlss.dist::pBCCG(y, mu = pa$M, sigma = pa$S, nu = pa$L)
  if (any(!is.finite(p) | p < 0 | p > 1)) stop('Invalid normalized BCCG percentile.')
  p
}

ref_audit <- function(ob, cfg) {

  fv <- lapply(c('mu', 'sigma', 'nu'), function(p) as.numeric(ob$fit[[paste0(p, '.fv')]]))
  if (any(lengths(fv) != nrow(ob$data))) stop('Missing fitted parameters for prediction audit.')

  ii <- !duplicated(ob$data$age)
  p0_unique <- data.frame(age = ob$data$age[ii],
                          M = fv[[1]][ii], S = fv[[2]][ii], L = fv[[3]][ii])
  p1_unique <- ref_predict(ob, p0_unique$age)
  q0 <- ref_quantiles(p0_unique, cfg$probs)
  q1 <- ref_quantiles(p1_unique, cfg$probs)

  rel_q <- max(abs(q1 - q0) / q0)

  p0_row <- data.frame(age = ob$data$age,
                       M = fv[[1]], S = fv[[2]], L = fv[[3]])
  jj <- match(ob$data$age, p1_unique$age)
  if (anyNA(jj)) stop('Prediction audit could not map re-predicted parameters to fitted rows.')
  p1_row <- p1_unique[jj, , drop = FALSE]
  pct0 <- ref_cdf(ob$data$y, p0_row)
  pct1 <- ref_cdf(ob$data$y, p1_row)
  pct_abs <- abs(pct1 - pct0)
  max_abs_pct <- max(pct_abs)
  mean_abs_pct <- mean(pct_abs)

  rel_m <- max(abs(p1_unique$M - p0_unique$M) / p0_unique$M)
  rel_s <- max(abs(p1_unique$S - p0_unique$S) / p0_unique$S)
  abs_l <- max(abs(p1_unique$L - p0_unique$L))

  centile_tol <- ref_opt(cfg, 'lms_prediction_tolerance', 1e-3)
  percentile_tol <- ref_opt(cfg, 'lms_percentile_tolerance', 1e-3)
  parameter_tol <- ref_opt(cfg, 'lms_parameter_diagnostic_tolerance', centile_tol)

  output_passed <- all(is.finite(c(rel_q, max_abs_pct, mean_abs_pct))) &&
    rel_q <= centile_tol && max_abs_pct <= percentile_tol
  parameter_diagnostic_passed <- all(is.finite(c(rel_m, rel_s, abs_l))) &&
    max(rel_m, rel_s) <= parameter_tol

  warns <- unique(attr(p1_unique, 'warnings'))

  data.frame(
    max_relative_M = rel_m,
    max_relative_S = rel_s,
    max_absolute_L = abs_l,
    max_relative_centile = rel_q,
    max_absolute_percentile = max_abs_pct,
    max_percentile_point_difference = 100 * max_abs_pct,
    mean_absolute_percentile = mean_abs_pct,
    centile_tolerance = centile_tol,
    percentile_tolerance = percentile_tol,
    parameter_diagnostic_tolerance = parameter_tol,
    output_passed = output_passed,
    parameter_diagnostic_passed = parameter_diagnostic_passed,
    diagnostic_review_required = !parameter_diagnostic_passed,
    tolerance = centile_tol,
    passed = output_passed,
    audit_rule = 'hard: centile relative error + percentile absolute error; diagnostic: M/S/L re-prediction',
    warnings = paste(warns, collapse = ' | '),
    stringsAsFactors = FALSE
  )
}

ref_fit <- function(d, spec, cfg) {
  d <- d[order(d$age, d$id, d$row_id), , drop = FALSE]
  if (any(!is.finite(d$y) | d$y <= 0)) stop('BCCG requires finite positive observations.')
  if (length(unique(d$age)) < 8L) stop('BCCG requires at least eight distinct ages.')
  formula_mu <- ref_formula(spec$df_mu, TRUE)
  formula_sigma <- ref_formula(spec$df_sigma)
  formula_nu <- ref_formula(spec$df_nu)
  budgets <- c(ref_opt(cfg, 'lms_n_cyc', 100L), ref_opt(cfg, 'lms_retry_n_cyc', 250L))
  history <- list()
  for (attempt in seq_along(budgets)) {
    warns <- character(); err <- ''; started <- proc.time()[['elapsed']]
    fit <- withCallingHandlers(tryCatch(
      gamlss::gamlss(formula = formula_mu, sigma.formula = formula_sigma,
        nu.formula = formula_nu,
        family = gamlss.dist::BCCG(mu.link = 'log', sigma.link = 'log', nu.link = 'identity'),
        data = d, control = gamlss::gamlss.control(n.cyc = budgets[attempt],
          c.crit = ref_opt(cfg, 'lms_c_crit', 0.001), trace = FALSE)),
      error = function(e) { err <<- conditionMessage(e); NULL }),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart('muffleWarning') })
    nonconv <- any(grepl('algorithm.*not.*converged|not yet converged', warns, ignore.case = TRUE))
    ok <- !is.null(fit) && isTRUE(fit$converged) && is.finite(fit$G.deviance) && !nonconv
    if (ok) {
      ob <- list(fit = fit, data = d, spec = spec)
      tst <- tryCatch(ref_predict(ob, sort(unique(d$age))), error = function(e) {
        err <<- conditionMessage(e); NULL
      })
      ok <- !is.null(tst)
      if (ok) warns <- c(warns, attr(tst, 'warnings'))
    }
    if (!ok && !nzchar(err)) err <- 'RS algorithm did not converge cleanly or deviance was nonfinite.'
    history[[attempt]] <- data.frame(candidate = as.character(spec$candidate),
      attempt = attempt, n_cyc = budgets[attempt], converged = ok,
      elapsed_sec = proc.time()[['elapsed']] - started,
      warnings = paste(unique(warns), collapse = ' | '), error = err, stringsAsFactors = FALSE)
    if (ok) return(list(ok = TRUE, model = ob, history = ref_bind(history), error = ''))
  }
  list(ok = FALSE, model = NULL, history = ref_bind(history), error = err)
}

ref_select <- function(d, cfg) {
  sp <- cfg$lms_candidates
  fits <- vector('list', nrow(sp)); rows <- vector('list', nrow(sp)); history <- list()
  for (i in seq_len(nrow(sp))) {
    st <- proc.time()[['elapsed']]
    z <- tryCatch(ref_fit(d, sp[i, , drop = FALSE], cfg), error = function(e)
      list(ok = FALSE, model = NULL, history = data.frame(), error = conditionMessage(e)))
    fits[i] <- list(z$model); history[[i]] <- z$history
    df <- if (z$ok) as.numeric(z$model$fit$df.fit) else NA_real_
    dev <- if (z$ok) as.numeric(z$model$fit$G.deviance) else NA_real_
    rows[[i]] <- data.frame(candidate = as.character(sp$candidate[i]),
      df_mu = sp$df_mu[i], df_sigma = sp$df_sigma[i], df_nu = sp$df_nu[i],
      formula_mu = paste(deparse(ref_formula(sp$df_mu[i], TRUE)), collapse = ' '),
      formula_sigma = paste(deparse(ref_formula(sp$df_sigma[i])), collapse = ' '),
      formula_nu = paste(deparse(ref_formula(sp$df_nu[i])), collapse = ' '),
      complexity = sp$df_mu[i] + sp$df_sigma[i] + sp$df_nu[i],
      fitted_df = df, global_deviance = dev, GAIC3 = dev + 3 * df,
      AIC = dev + 2 * df, converged = z$ok,
      elapsed_sec = proc.time()[['elapsed']] - st, error = z$error, stringsAsFactors = FALSE)
  }
  rows <- ref_bind(rows)
  valid <- which(rows$converged & is.finite(rows$GAIC3))
  if (!length(valid)) return(list(ok = FALSE, rows = rows, history = ref_bind(history),
    model = NULL, error = 'All prespecified BCCG candidates failed.'))
  rows$delta_GAIC3 <- rows$GAIC3 - min(rows$GAIC3[valid])
  eligible <- valid[rows$delta_GAIC3[valid] <= 2 + 1e-10]
  pick <- eligible[order(rows$complexity[eligible], rows$GAIC3[eligible], rows$candidate[eligible])][1]
  rows$selected <- seq_len(nrow(rows)) == pick
  list(ok = TRUE, rows = rows, history = ref_bind(history), model = fits[[pick]], error = '')
}

ref_folds <- function(d, k, seed) {
  u <- d[!duplicated(d$id), c('id', 'age_group'), drop = FALSE]
  if (nrow(u) < k * 2L) stop('Too few participant IDs for requested OOF folds.')
  set.seed(seed); u$fold <- NA_integer_; offset <- 0L
  for (ix in split(seq_len(nrow(u)), u$age_group)) {
    labs <- (seq_along(ix) - 1L + offset) %% k + 1L
    u$fold[ix] <- sample(labs, length(labs), replace = FALSE)
    offset <- offset + length(ix)
  }
  u$fold[match(d$id, u$id)]
}

ref_fold_task <- function(fold, d, cfg) {
  train <- d[d$fold != fold, , drop = FALSE]
  test <- d[d$fold == fold, , drop = FALSE]
  st <- proc.time()[['elapsed']]; select <- NULL; audit <- data.frame(); warns <- character()
  err <- ''; score <- data.frame()
  supported <- test$age >= min(train$age) & test$age <= max(train$age)
  ok <- tryCatch({
    if (length(intersect(unique(train$id), unique(test$id)))) stop('Participant leakage across OOF train/test split.')
    select <- ref_select(train, cfg)
    if (!select$ok) stop(select$error)
    audit <- ref_audit(select$model, cfg)
    if (!isTRUE(audit$passed)) stop('OOF selected model failed output prediction audit.')
    score <- data.frame(row_id = test$row_id, id = test$id, age = test$age,
      age_group = test$age_group, side = as.character(test$side), value = test$y,
      percentile = NA_real_, source = 'OOF', fold = fold, negative_log_density = NA_real_,
      M = NA_real_, S = NA_real_, L = NA_real_, candidate = as.character(select$model$spec$candidate),
      score_available = supported, score_exclusion = ifelse(supported, '', 'outside_training_age_support'),
      stringsAsFactors = FALSE)
    if (any(supported)) {
      pa <- ref_predict(select$model, test$age[supported]); warns <- attr(pa, 'warnings')
      ll <- gamlss.dist::dBCCG(test$y[supported], mu = pa$M, sigma = pa$S, nu = pa$L, log = TRUE)
      if (any(!is.finite(ll))) stop('Nonfinite OOF log-density.')
      score$percentile[supported] <- ref_cdf(test$y[supported], pa)
      score$negative_log_density[supported] <- -ll
      score[supported, c('M','S','L')] <- pa[, c('M','S','L')]
    }
    TRUE
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  list(ok = ok, scores = score, candidates = if (is.null(select)) data.frame() else select$rows,
    history = if (is.null(select)) data.frame() else select$history, audit = audit,
    ledger = data.frame(fold = fold, n_train_ids = length(unique(train$id)),
      n_test_ids = length(unique(test$id)), n_test_feet = nrow(test), success = ok,
      train_age_min = min(train$age), train_age_max = max(train$age),
      n_test_outside_age_support = sum(!supported),
      n_test_scored = if (nrow(score)) sum(is.finite(score$percentile)) else 0L,
      participant_overlap = length(intersect(unique(train$id), unique(test$id))),
      elapsed_sec = proc.time()[['elapsed']] - st, warnings = paste(warns, collapse = ' | '),
      error = err, stringsAsFactors = FALSE))
}

ref_resample <- function(d, seed) {

  u <- d[!duplicated(d$id), c('id', 'age_group'), drop = FALSE]
  byid <- split(seq_len(nrow(d)), as.character(d$id))
  set.seed(seed)
  draw <- unlist(lapply(split(as.character(u$id), u$age_group), function(ids)
    ids[sample.int(length(ids), length(ids), replace = TRUE)]), use.names = FALSE)
  ix <- unlist(byid[draw], use.names = FALSE)
  out <- d[ix, , drop = FALSE]
  out$id <- rep(paste0('bootstrap_', seq_along(draw)), lengths(byid[draw]))
  rownames(out) <- NULL
  out
}

ref_boot_task <- function(b, d, spec, cfg, seed) {
  st <- proc.time()[['elapsed']]; fit <- NULL; audit <- data.frame(); q <- NULL
  err <- ''; warns <- character()
  age_range <- c(NA_real_, NA_real_); edge_n <- NA_integer_; edge_distance <- NA_real_
  ok <- tryCatch({
    if (min(cfg$age_grid) < min(d$age) || max(cfg$age_grid) > max(d$age))
      stop('Bootstrap grid must lie inside original observed age support.')
    boot <- ref_resample(d, seed + b)
    age_range <- range(boot$age)
    edge_n <- sum(cfg$age_grid < age_range[1] | cfg$age_grid > age_range[2])
    edge_distance <- max(0, age_range[1] - min(cfg$age_grid), max(cfg$age_grid) - age_range[2])
    fit <- ref_fit(boot, spec, cfg)
    if (!fit$ok) stop(fit$error)
    audit <- ref_audit(fit$model, cfg)
    if (!isTRUE(audit$passed)) stop('Bootstrap model failed output prediction audit.')
    pa <- ref_predict(fit$model, cfg$age_grid); warns <- attr(pa, 'warnings')
    q <- ref_quantiles(pa, cfg$probs)
    TRUE
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  list(ok = ok, q = if (ok) q else NULL, audit = audit,
    history = if (is.null(fit)) data.frame() else fit$history,
    ledger = data.frame(replicate = b, success = ok,
      resampled_age_min = age_range[1], resampled_age_max = age_range[2],
      n_grid_points_outside_resampled_age_support = edge_n,
      max_distance_outside_resampled_age_support = edge_distance,
      edge_policy = 'fixed_original_support_grid_spline_edge_prediction_permitted',
      elapsed_sec = proc.time()[['elapsed']] - st,
      warnings = paste(warns, collapse = ' | '), error = err, stringsAsFactors = FALSE))
}

ref_pool <- function(cfg) {
  n <- as.integer(ref_opt(cfg, 'workers', 1L))
  if (n <= 1L) return(NULL)
  cl <- parallel::makePSOCKcluster(n)
  tryCatch({
    parallel::clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')
      if (!requireNamespace('gamlss', quietly = TRUE) || !requireNamespace('gamlss.dist', quietly = TRUE))
        stop('gamlss/gamlss.dist missing on a reference worker.')
      NULL
    })
    env <- environment(ref_pool)
    funs <- ls(envir = env, pattern = '^ref_')
    funs <- funs[vapply(funs, function(nm) is.function(get(nm, envir = env)), logical(1))]
    parallel::clusterExport(cl, funs, envir = env)
    cl
  }, error = function(e) { parallel::stopCluster(cl); stop(e) })
}

ref_map <- function(cl, X, FUN, ...) {
  if (is.null(cl)) lapply(X, FUN, ...) else parallel::parLapplyLB(cl, X, FUN, ...)
}

ref_hash <- function(object) {
  path <- tempfile(fileext = '.rds'); on.exit(unlink(path))
  saveRDS(object, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

ref_bootstrap <- function(d, model, cfg, cl, key, out, seed) {
  B <- as.integer(cfg$n_boot)
  path <- file.path(out, 'checkpoints', paste0('reference_bootstrap_', key, '.rds'))
  code_signature <- lapply(list(ref_fit, ref_resample, ref_predict, ref_quantiles,
    ref_audit, ref_boot_task), function(f) paste(deparse(body(f)), collapse = '\n'))
  signature <- ref_hash(list(version = 'reference_v1.2_output_audit', data = d,
    spec = model$spec, code = code_signature,
    cfg = cfg[setdiff(names(cfg), c('workers', 'output_dir'))], seed = seed))
  res <- vector('list', B)
  if (isTRUE(ref_opt(cfg, 'resume', TRUE)) && file.exists(path)) {
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(old) && identical(old$signature, signature) && length(old$results) == B)
      res <- old$results
  }
  remaining <- which(vapply(res, is.null, logical(1)))
  if (length(remaining)) for (batch in split(remaining, ceiling(seq_along(remaining) / 16L))) {
    rr <- ref_map(cl, batch, ref_boot_task, d = d, spec = model$spec, cfg = cfg, seed = seed)
    res[batch] <- rr
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(path, '.tmp')
    saveRDS(list(signature = signature, results = res), tmp)
    if (!file.rename(tmp, path)) { file.copy(tmp, path, overwrite = TRUE); unlink(tmp) }
    if (isTRUE(cfg$verbose)) message(sprintf('Reference %s bootstrap: %d/%d completed.', key, sum(!vapply(res, is.null, logical(1))), B))
  }
  res
}

ref_calibration <- function(scores, cfg) {
  if (!nrow(scores)) return(data.frame())
  groups <- split(seq_len(nrow(scores)), interaction(scores$source, scores$age_group, drop = TRUE))
  groups <- c(list(overall_OOF = which(scores$source == 'OOF'),
                  overall_apparent = which(scores$source == 'apparent')), groups)
  ans <- lapply(groups, function(ix) {
    if (!length(ix)) return(data.frame())
    all_z <- scores[ix, , drop = FALSE]
    z <- all_z[is.finite(all_z$percentile), , drop = FALSE]
    data.frame(source = all_z$source[1], age_group = if (length(unique(all_z$age_group)) == 1L)
      all_z$age_group[1] else NA_integer_, scope = if (length(unique(all_z$age_group)) == 1L) 'age' else 'overall',
      probability = cfg$probs, observed_fraction = vapply(cfg$probs, function(p)
        if (nrow(z)) mean(z$percentile <= p) else NA_real_, numeric(1)),
      n_feet = nrow(z), n_ids = length(unique(z$id)), n_feet_total = nrow(all_z),
      n_feet_unscored = nrow(all_z) - nrow(z),
      mean_negative_log_density = if (nrow(z)) mean(z$negative_log_density) else NA_real_, stringsAsFactors = FALSE)
  })
  ref_bind(ans)
}

run_reference <- function(d, cfg, out) {
  for (pkg in c('gamlss', 'gamlss.dist')) if (!requireNamespace(pkg, quietly = TRUE))
    stop('Reference analysis requires package: ', pkg)
  needed <- c('row_id', 'id', 'sex', 'age', 'age_group', 'side', 'FL', 'AI', 'AS_ratio')
  if (!all(needed %in% names(d))) stop('Missing reference columns: ', paste(setdiff(needed, names(d)), collapse = ', '))
  if (anyDuplicated(d$row_id)) stop('row_id must be unique before reference analysis.')
  if (any(cfg$probs <= 0 | cfg$probs >= 1) || is.unsorted(cfg$probs, strictly = TRUE))
    stop('cfg$probs must be strictly increasing probabilities inside (0, 1).')
  if (length(cfg$age_grid) < 2L || any(!is.finite(cfg$age_grid)) || is.unsorted(cfg$age_grid, strictly = TRUE))
    stop('Invalid age grid.')
  specs <- cfg$lms_candidates
  if (!all(c('candidate', 'df_mu', 'df_sigma', 'df_nu') %in% names(specs)) || !nrow(specs) ||
      anyDuplicated(specs$candidate) || any(!is.finite(as.matrix(specs[, c('df_mu','df_sigma','df_nu')]))))
    stop('Invalid prespecified LMS candidate table.')
  for (p in c('results/reference', 'models', 'qc', 'checkpoints'))
    dir.create(file.path(out, p), recursive = TRUE, showWarnings = FALSE)
  cl <- ref_pool(cfg)
  if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
  store <- list(scores = list(), centiles = list(), parameters = list(), calibration = list(),
    candidates = list(), history = list(), oof = list(), audits = list(), bootstrap = list(),
    exclusions = list(), status = list(), threshold_map = list(), grid_support = list(), distribution_support = list())
  models <- list(); counter <- 0L
  flush <- function() {
    for (nm in names(store)) ref_write(ref_bind(store[[nm]]),
      file.path(out, if (nm %in% c('exclusions','status','history','audits','bootstrap','oof','grid_support','distribution_support')) 'qc' else 'results/reference',
        paste0('reference_', nm, '.csv')))
  }
  for (metric in c('FL', 'AI', 'AS_ratio')) for (sex in c('Boys', 'Girls')) {
    counter <- counter + 1L; key <- paste(metric, sex, sep = '_')
    st <- proc.time()[['elapsed']]; err <- ''; model_ok <- FALSE; oof_ok <- FALSE
    boot_success <- 0L; ci_ok <- FALSE; selected <- NA_character_
    ai_tail_max <- NA_real_; ai_support_flag <- NA; n_oof_unscored <- NA_integer_
    group_cfg <- cfg
    z <- d[as.character(d$sex) == sex, needed, drop = FALSE]
    z$y <- z[[metric]]
    reason <- ifelse(!is.finite(z$y), 'missing_or_nonfinite', ifelse(z$y <= 0, 'nonpositive_BCCG_support', ''))
    invalid <- nzchar(reason)
    if (any(invalid)) store$exclusions[[key]] <- data.frame(metric = metric, sex = sex,
      row_id = z$row_id[invalid], id = z$id[invalid], age = z$age[invalid],
      value = z$y[invalid], reason = reason[invalid], stringsAsFactors = FALSE)
    z <- z[!invalid, , drop = FALSE]
    z <- z[, c('row_id', 'id', 'age', 'age_group', 'side', 'y'), drop = FALSE]
    annotate <- function(x, stage = NULL, fold = NULL, replicate = NULL) {
      if (!is.data.frame(x) || !nrow(x)) return(data.frame())
      x$metric <- metric; x$sex <- sex
      if (!is.null(stage)) x$stage <- stage
      if (!is.null(fold)) x$fold <- fold
      if (!is.null(replicate)) x$replicate <- replicate
      x
    }
    if (isTRUE(cfg$verbose)) message('Reference model: ', key)
    tryCatch({
      if (sum(reason == 'nonpositive_BCCG_support') > 0 &&
          identical(ref_opt(cfg, 'lms_nonpositive_policy', 'exclude'), 'fail'))
        stop('Nonpositive values present; configured BCCG policy is fail.')
      if (length(unique(z$id)) < ref_opt(cfg, 'min_clusters', 30L)) stop('Too few usable participant clusters.')
      in_support <- cfg$age_grid >= min(z$age) & cfg$age_grid <= max(z$age)
      store$grid_support[[key]] <- annotate(data.frame(age = cfg$age_grid, retained = in_support,
        observed_age_min = min(z$age), observed_age_max = max(z$age),
        reason = ifelse(in_support, '', 'outside_original_observed_age_support')))
      group_cfg$age_grid <- cfg$age_grid[in_support]
      if (length(group_cfg$age_grid) < 2L) stop('Fewer than two requested grid ages lie inside observed support.')
      sel <- ref_select(z, cfg)
      store$candidates[[key]] <- annotate(sel$rows, 'full_data')
      store$history[[key]] <- annotate(sel$history, 'full_data')
      if (!sel$ok) stop(sel$error)
      ob <- sel$model; selected <- as.character(ob$spec$candidate)
      audit <- ref_audit(ob, cfg); store$audits[[key]] <- annotate(audit, 'full_data')
      if (!isTRUE(audit$passed)) stop('Full-data selected model failed numerical prediction audit.')
      pa <- ref_predict(ob, group_cfg$age_grid); qq <- ref_quantiles(pa, cfg$probs)
      model_ok <- TRUE; models[[key]] <- ob
      saveRDS(ob, file.path(out, 'models', paste0('reference_', key, '.rds')))
      pa$P50 <- gamlss.dist::qBCCG(0.5, mu = pa$M, sigma = pa$S, nu = pa$L)
      store$parameters[[key]] <- annotate(pa)
      if (metric == 'AI') {
        tail_mass <- gamlss.dist::pBCCG(1, mu = pa$M, sigma = pa$S, nu = pa$L, lower.tail = FALSE)
        ai_tail_max <- max(tail_mass)
        threshold <- ref_opt(cfg, 'ai_tail_mass_threshold', 0.001)
        ai_support_flag <- any(tail_mass > threshold)
        store$distribution_support[[key]] <- annotate(data.frame(age = pa$age,
          probability_above_one = tail_mass, threshold = threshold, review_required = tail_mass > threshold,
          highest_requested_centile_above_one = qq[, ncol(qq)] > 1,
          action = 'diagnostic_only_no_clipping_or_renormalization'))
      }
      if (metric == 'AI') store$threshold_map[[key]] <- annotate(ref_bind(lapply(c(0.21, 0.26), function(th) {
        data.frame(age = pa$age, threshold = th, reference_percentile = gamlss.dist::pBCCG(th,
          mu = pa$M, sigma = pa$S, nu = pa$L), source = 'apparent', stringsAsFactors = FALSE)
      })))
      obs <- ref_predict(ob, z$age)
      apparent <- data.frame(row_id = z$row_id, id = z$id, age = z$age, age_group = z$age_group,
        side = as.character(z$side), value = z$y, percentile = ref_cdf(z$y, obs), source = 'apparent',
        fold = NA_integer_, negative_log_density = -gamlss.dist::dBCCG(z$y, mu = obs$M,
          sigma = obs$S, nu = obs$L, log = TRUE), M = obs$M, S = obs$S, L = obs$L,
        candidate = selected, score_available = TRUE, score_exclusion = '', stringsAsFactors = FALSE)
      z$fold <- ref_folds(z, cfg$folds, cfg$seed + counter * 1000L)
      fo <- ref_map(cl, seq_len(cfg$folds), ref_fold_task, d = z, cfg = cfg)
      oo <- ref_bind(lapply(fo, `[[`, 'scores'))
      n_oof_unscored <- nrow(z) - sum(is.finite(oo$percentile))
      oof_ok <- all(vapply(fo, `[[`, logical(1), 'ok')) && nrow(oo) == nrow(z) &&
        !anyDuplicated(oo$row_id) && n_oof_unscored == 0L
      for (f in seq_along(fo)) {
        tag <- paste(key, 'fold', f, sep = '_')
        store$candidates[[tag]] <- annotate(fo[[f]]$candidates, 'OOF_training', f)
        store$history[[tag]] <- annotate(fo[[f]]$history, 'OOF_training', f)
        store$audits[[tag]] <- annotate(fo[[f]]$audit, 'OOF_training', f)
        store$oof[[tag]] <- annotate(fo[[f]]$ledger)
      }
      scores <- ref_bind(list(apparent, oo))
      store$scores[[key]] <- annotate(scores)
      cal <- ref_calibration(scores, cfg)
      if (nrow(cal)) cal$OOF_complete <- oof_ok
      store$calibration[[key]] <- annotate(cal)
      if (!oof_ok) err <- 'OOF incomplete: review fold ledger; no apparent replacement was made.'
      boot <- ref_bootstrap(z, ob, group_cfg, cl, key, out, cfg$seed + counter * 100000L)
      good <- which(vapply(boot, `[[`, logical(1), 'ok')); boot_success <- length(good)
      ci_ok <- cfg$n_boot > 0L && boot_success / cfg$n_boot >= cfg$min_boot_success && boot_success >= 2L
      lower <- upper <- matrix(NA_real_, nrow(qq), ncol(qq))
      if (ci_ok) {
        mat <- vapply(boot[good], function(v) as.vector(v$q), numeric(length(qq)))
        lower[] <- apply(mat, 1L, stats::quantile, probs = 0.025, names = FALSE, type = 7)
        upper[] <- apply(mat, 1L, stats::quantile, probs = 0.975, names = FALSE, type = 7)
      }
      for (b in seq_along(boot)) {
        tag <- paste(key, 'boot', b, sep = '_')
        store$bootstrap[[tag]] <- annotate(boot[[b]]$ledger)
        store$history[[tag]] <- annotate(boot[[b]]$history, 'bootstrap', replicate = b)
        store$audits[[tag]] <- annotate(boot[[b]]$audit, 'bootstrap', replicate = b)
      }
      tab <- expand.grid(age = group_cfg$age_grid, probability = cfg$probs, KEEP.OUT.ATTRS = FALSE)
      tab$estimate <- as.vector(qq); tab$lower <- as.vector(lower); tab$upper <- as.vector(upper)
      tab$natural_support_violation <- if (metric == 'AI') tab$estimate > 1 else FALSE
      tab$AI_probability_above_one <- if (metric == 'AI') rep(tail_mass, length(cfg$probs)) else NA_real_
      tab$AI_support_review_required <- if (metric == 'AI')
        tab$AI_probability_above_one > ref_opt(cfg, 'ai_tail_mass_threshold', 0.001) else NA
      tab$n_boot_requested <- cfg$n_boot; tab$n_boot_success <- boot_success
      tab$ci_available <- ci_ok; tab$ci_type <- 'pointwise_percentile_95_percent'
      tab$inference <- 'age_stratified_participant_bootstrap_conditional_on_selected_smoothing'
      tab$bootstrap_edge_policy <- 'fixed_original_support_grid_spline_edge_prediction_permitted'
      tab$candidate <- selected; tab$population <- 'observed_positive_feet'
      store$centiles[[key]] <- annotate(tab)
      if (!ci_ok) err <- paste(c(err[nzchar(err)], 'Bootstrap success requirement unmet; CI withheld.'), collapse = ' | ')
    }, error = function(e) { err <<- paste(c(err[nzchar(err)], conditionMessage(e)), collapse = ' | ') })
    store$status[[key]] <- data.frame(metric = metric, sex = sex, n_feet = nrow(z),
      n_ids = length(unique(z$id)), n_excluded = sum(invalid), selected_candidate = selected,
      model_success = model_ok, OOF_complete = oof_ok, bootstrap_requested = cfg$n_boot,
      n_OOF_unscored = n_oof_unscored, n_grid_requested = length(cfg$age_grid),
      n_grid_retained = length(group_cfg$age_grid),
      AI_max_probability_above_one = ai_tail_max, AI_support_review_required = ai_support_flag,
      bootstrap_success = boot_success, bootstrap_success_fraction = if (cfg$n_boot > 0) boot_success/cfg$n_boot else NA_real_,
      CI_available = ci_ok, elapsed_sec = proc.time()[['elapsed']] - st, error = err,
      stringsAsFactors = FALSE)
    flush()
  }
  ans <- lapply(store, ref_bind); ans$models <- models
  ans
}
