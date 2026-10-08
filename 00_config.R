foot_config <- function(root = ".") {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  nc <- parallel::detectCores(logical = FALSE)
  if (is.na(nc)) nc <- 2L
  list(
    version = "foot-stats-1.1.0",
    data_file = file.path(root, "data", "observations.csv"),
    output_dir = file.path(root, "results"),
    stage = "full",
    resume = TRUE, verbose = FALSE,
    expected_children = 4929L, expected_feet = 9858L,
    require_expected_n = TRUE,
    encoding = "UTF-8-BOM",
    side_labels = c("1" = "Left", "2" = "Right"),

    age_boundary = c(3, 19), age_knots = c(6, 9, 12, 15),
    age_grid = seq(3, 18, by = .05), probs = c(.03,.10,.25,.50,.75,.90,.97),
    seed = 20260928L, folds = 5L, min_clusters = 30L,
    workers = max(1L, min(6L, nc - 1L)),
    n_boot = 500L, av_n_boot = 500L, prop_n_boot = 2000L,
    n_sim = 5000L, curve_n = 101L, min_boot_success = .95,
    av_age_n_sim = 5000L, av_age_seed = 20261001L, av_age_derivative_h = .01,
    lms_candidates = data.frame(
      candidate=c("M3_S1_L0","M4_S1_L0","M4_S2_L0","M3_S1_L1","M4_S1_L1","M5_S2_L1"),
      df_mu=c(3,4,4,3,4,5), df_sigma=c(1,1,2,1,1,2), df_nu=c(0,0,0,1,1,1)),
    bmi_tolerance = .05, ratio_tolerance = 1e-4,
    allow_derived_discrepancies = FALSE,
    min_tail_years = 1,

    equivalence_margin = c(FL=NA_real_, AI=NA_real_, AS_ratio=NA_real_, AV_sizeadjusted=NA_real_),
    export_individual = FALSE,
    run_associations = TRUE, run_international = FALSE,
    international_dir = NULL
  )
}
