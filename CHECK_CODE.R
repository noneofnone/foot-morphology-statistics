local({
  files <- Filter(function(x) !is.null(x), lapply(sys.frames(), function(x) x$ofile))
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  self <- if (length(files)) tail(files, 1)[[1]] else if (length(cli)) sub("^--file=", "", cli[1]) else "CHECK_CODE.R"
  root <- dirname(normalizePath(self, mustWork = TRUE))
  scripts <- c(list.files(root, pattern = "\\.R$", full.names = TRUE),
               list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE),
               list.files(file.path(root, "tests"), pattern = "\\.R$", full.names = TRUE))
  for (path in scripts) parse(path, encoding = "UTF-8")
  logs <- file.path(root, "validation")
  dir.create(logs, showWarnings = FALSE)
  exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  tests <- list.files(file.path(root, "tests"), pattern = "^test_.*\\.R$", full.names = TRUE)
  status <- lapply(tests, function(path) {
    log <- file.path(logs, paste0(basename(path), ".log"))
    started <- proc.time()[["elapsed"]]
    code <- system2(exe, c(shQuote(path), shQuote(root)), stdout = log, stderr = log)
    if (code != 0L) stop("Test failed: ", basename(path), "; see ", log)
    data.frame(test = basename(path), exit_code = code,
               elapsed_seconds = proc.time()[["elapsed"]] - started)
  })
  write.csv(do.call(rbind, status), file.path(logs, "test_summary.csv"), row.names = FALSE)
  invisible(status)
})
