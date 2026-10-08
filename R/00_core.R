`%||%` <- function(x, y) if (is.null(x)) y else x
bind_rows_base <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0, x)
  if (!length(x)) return(data.frame())
  nm <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(z) { for (n in setdiff(nm, names(z))) z[[n]] <- NA; z[,nm,drop=FALSE] })
  ans <- do.call(rbind, x); rownames(ans) <- NULL; ans
}
write_table <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  if (!is.data.frame(x)) x <- as.data.frame(x)
  utils::write.csv(x, path, row.names=FALSE, na="", fileEncoding="UTF-8")
  invisible(path)
}
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(path,".tmp"); saveRDS(x,tmp,compress=FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp,path)) stop("Cannot write checkpoint: ",path)
}
object_hash <- function(x) {
  p <- tempfile(); on.exit(unlink(p)); saveRDS(x,p,compress=FALSE)
  unname(tools::md5sum(p))
}
assign_folds <- function(d, k, seed) {
  set.seed(seed)
  u <- d[!duplicated(d$id),c("id","sex","age_group"),drop=FALSE]
  f <- integer(nrow(u))
  for (ii in split(seq_len(nrow(u)), interaction(u$sex,u$age_group,drop=TRUE)))
    f[ii] <- sample(rep(seq_len(k),length.out=length(ii)))
  f[match(d$id,u$id)]
}
cluster_resample <- function(d, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- d[!duplicated(d$id),c("id","sex","age_group"),drop=FALSE]
  rows <- split(seq_len(nrow(d)),as.character(d$id))
  drawn <- unlist(lapply(split(as.character(u$id),interaction(u$sex,u$age_group,drop=TRUE)),
    function(v) v[sample.int(length(v),length(v),replace=TRUE)]),use.names=FALSE)
  rr <- rows[drawn]; z <- d[unlist(rr,use.names=FALSE),,drop=FALSE]
  z$id <- rep(sprintf("B%07d",seq_along(rr)),lengths(rr)); rownames(z) <- NULL
  z
}
normal_draws <- function(beta,V,n_sim,seed) {
  V <- (V+t(V))/2; e <- eigen(V,symmetric=TRUE)
  if(min(e$values) < -1e-7*max(abs(e$values),1e-12)) stop("Covariance matrix is not positive semidefinite")
  set.seed(seed)
  sweep(e$vectors %*% (sqrt(pmax(e$values,0))*matrix(rnorm(length(beta)*n_sim),nrow=length(beta))),1,beta,"+")
}
linear_bands <- function(beta,V,L,n_sim=5000L,seed=1L) {
  estimate <- drop(L%*%beta); se <- sqrt(pmax(rowSums((L%*%V)*L),0))
  draws <- normal_draws(rep(0,length(beta)),V,n_sim,seed)
  nz <- se > sqrt(.Machine$double.eps)*max(1,max(se))
  crit <- qnorm(.975)
  if(any(nz)) {

    maxima <- numeric(n_sim)
    for (ii in split(seq_len(n_sim), ceiling(seq_len(n_sim)/500)))
      maxima[ii] <- apply(abs((L[nz,,drop=FALSE]%*%draws[,ii,drop=FALSE])/se[nz]),2,max)
    crit <- unname(quantile(maxima,.95))
  }
  data.frame(estimate,se,lower=estimate-qnorm(.975)*se,upper=estimate+qnorm(.975)*se,
    sim_lower=estimate-crit*se,sim_upper=estimate+crit*se,sim_critical=crit,
    p=ifelse(se>0,2*pnorm(-abs(estimate/se)),NA_real_))
}
fit_gee <- function(d, formula, cfg=list(min_clusters=30L), corstr="exchangeable") {
  vars <- unique(c(all.vars(formula),"id","side"))
  z <- droplevels(d[complete.cases(d[,vars,drop=FALSE]),,drop=FALSE])
  z <- z[order(z$id,z$side),,drop=FALSE]
  if(length(unique(z$id)) < (cfg$min_clusters %||% 30L)) stop("Insufficient complete participant clusters")
  if(nlevels(z$side)==1L && "side" %in% all.vars(formula)) formula <- update(formula,.~.-side)
  X <- model.matrix(formula,z)
  if(qr(X)$rank<ncol(X)) stop("Rank deficient GEE design")
  z$cluster_index <- match(z$id,unique(z$id))
  if(max(table(z$id))==1L) corstr <- "independence"
  f <- geepack::geeglm(formula,data=z,id=cluster_index,corstr=corstr,family=gaussian(),
    std.err="san.se",control=geepack::geese.control(maxit=100))
  if(any(!is.finite(coef(f)))||any(!is.finite(f$geese$vbeta))||f$geese$error!=0) stop("GEE convergence failure")
  V <- f$geese$vbeta; dimnames(V) <- list(names(coef(f)),names(coef(f)))
  list(fit=f,beta=coef(f),V=V,data=z,n_children=length(unique(z$id)),n_feet=nrow(z),
    corstr=corstr,alpha=f$geese$alpha)
}
design <- function(obj,nd) {
  X <- model.matrix(delete.response(terms(obj$fit)),nd,contrasts.arg=obj$fit$contrasts,xlev=obj$fit$xlevels)
  X[,names(obj$beta),drop=FALSE]
}
wald_block <- function(obj, indices) {
  if(!length(indices)) return(c(statistic=NA_real_,df=0,p=NA_real_))
  b <- obj$beta[indices]; V <- obj$V[indices,indices,drop=FALSE]
  e <- eigen((V+t(V))/2,symmetric=TRUE); keep <- e$values>max(e$values)*1e-8
  if(!any(keep)) return(c(statistic=NA_real_,df=0,p=NA_real_))
  w <- sum(drop(crossprod(e$vectors[,keep,drop=FALSE],b))^2/e$values[keep])
  c(statistic=w,df=sum(keep),p=pchisq(w,sum(keep),lower.tail=FALSE))
}
required_packages <- function() c("mgcv","geepack","gamlss","gamlss.dist","sandwich")
check_packages <- function() {
  p <- required_packages(); missing <- p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Missing R packages. Run INSTALL_PACKAGES.R first: ",paste(missing,collapse=", "))
  invisible(TRUE)
}
