pkgs <- c("mgcv","geepack","gamlss","gamlss.dist","sandwich")
missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)))
