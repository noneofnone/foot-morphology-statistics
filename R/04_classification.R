cls_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
  if (!length(x)) return(data.frame())
  nm <- unique(unlist(lapply(x, names)))
  x <- lapply(x,function(z) {for (key in setdiff(nm,names(z))) z[[key]] <- NA; z[,nm,drop=FALSE]})
  ans <- do.call(rbind, x); rownames(ans) <- NULL; ans
}
cls_config <- function(cfg, key, default) {
  if (is.null(cfg[[key]])) default else cfg[[key]]
}
cls_groups <- function(d, by_age = TRUE) {
  ans <- list(list(data = d, scope = 'overall', sex = 'All', age_group = NA_character_))
  for (s in unique(as.character(d$sex))) {
    ds <- d[as.character(d$sex) == s, , drop = FALSE]
    ans[[length(ans) + 1L]] <- list(data = ds, scope = 'sex', sex = s, age_group = NA_character_)
    if (by_age) for (a in unique(as.character(ds$age_group))) {
      if (is.na(a)) next
      ans[[length(ans) + 1L]] <- list(data = ds[as.character(ds$age_group) == a, , drop = FALSE],
                                    scope = 'age_sex', sex = s, age_group = a)
    }
  }
  ans
}
cls_meta <- function(tab, g) {
  if (!nrow(tab)) return(tab)
  tab$scope <- g$scope; tab$sex <- g$sex; tab$age_group <- g$age_group
  tab
}
cls_contrib <- function(id, x) {
  x <- as.matrix(x); storage.mode(x) <- 'double'
  if (!length(id)) return(matrix(numeric(), 0L, ncol(x)))
  rowsum(x, group = factor(as.character(id), levels = unique(as.character(id))), reorder = FALSE)
}
cls_strata <- function(d) {
  if (anyNA(d$sex) || anyNA(d$age_group)) stop('Bootstrap stratification requires complete sex and age_group.')
  paste(as.character(d$sex), as.character(d$age_group), sep='::')
}
cls_child_strata <- function(id, strata, child_ids) {
  if (is.null(strata)) return(rep('single_stratum',length(child_ids)))
  if (length(strata) != length(id) || anyNA(strata)) stop('Invalid bootstrap strata.')
  id <- as.character(id); first <- !duplicated(id)
  child_strata <- as.character(strata[match(child_ids,id)])
  if (any(as.character(strata) != child_strata[match(id,child_ids)]))
    stop('Sex or age-group stratum varies within a child.')
  child_strata
}
cls_boot_weights <- function(child_strata,n_rep) {
  ans <- matrix(0L,length(child_strata),n_rep)
  for (ix in split(seq_along(child_strata),child_strata))
    ans[ix,] <- stats::rmultinom(n_rep,length(ix),rep.int(1/length(ix),length(ix)))
  ans
}
cls_ci <- function(x, estimate, boundary = FALSE) {
  x <- x[is.finite(x)]
  if (!is.finite(estimate)) return(c(lower = NA_real_, upper = NA_real_, n_boot_valid = length(x)))
  if (boundary || length(x) < 20L) return(c(lower = NA_real_, upper = NA_real_, n_boot_valid = length(x)))
  c(lower = unname(stats::quantile(x, .025, type = 7)),
    upper = unname(stats::quantile(x, .975, type = 7)), n_boot_valid = length(x))
}

cls_ratios <- function(id, numerator, denominator, labels, cfg, seed, strata = NULL) {
  numerator <- as.matrix(numerator); denominator <- as.matrix(denominator)
  if (!length(id)) return(data.frame())
  if (ncol(denominator) == 1L && ncol(numerator) > 1L)
    denominator <- denominator[, rep(1L, ncol(numerator)), drop = FALSE]
  stopifnot(identical(dim(numerator), dim(denominator)), ncol(numerator) == length(labels),
            !anyNA(numerator), !anyNA(denominator), all(numerator >= 0),
            all(numerator <= denominator))
  k <- ncol(numerator); cn <- cls_contrib(id, numerator); cd <- cls_contrib(id, denominator)
  nchild <- nrow(cn); child_strata <- cls_child_strata(id,strata,rownames(cn))
  n <- colSums(numerator); den <- colSums(denominator)
  est <- n / den; est[den == 0] <- NA_real_
  b <- as.integer(cls_config(cfg, 'prop_n_boot', 2000L))
  chunk <- as.integer(cls_config(cfg, 'bootstrap_chunk_size', 200L))
  boot <- matrix(NA_real_, nrow = b, ncol = k)
  set.seed(seed)
  if (nchild >= 2L && b > 0L) for (start in seq.int(1L, b, by = chunk)) {
    ix <- start:min(start + chunk - 1L, b)
    w <- cls_boot_weights(child_strata,length(ix))
    nb <- crossprod(w, cn); db <- crossprod(w, cd)
    z <- nb / db; z[db == 0] <- NA_real_; boot[ix, ] <- z
  }
  boundary <- den > 0 & (n == 0 | n == den)
  ci <- t(vapply(seq_len(k), function(j) cls_ci(boot[, j], est[j], boundary[j]), numeric(3)))
  ans <- data.frame(measure = labels, numerator_n = n, denominator_n = den,
                    proportion = est, percent = 100 * est, lower = ci[,1], upper = ci[,2],
                    lower_percent = 100 * ci[,1], upper_percent = 100 * ci[,2],
                    n_children = colSums(cd > 0), n_boot_valid = ci[,3],
                    ci_status = ifelse(den == 0, 'no_eligible_observations',
                      ifelse(nchild < 2L, 'fewer_than_two_children',
                        ifelse(boundary, 'boundary_degenerate',
                          ifelse(ci[,3] < 20, 'insufficient_valid_bootstrap', 'ok')))),
                    stringsAsFactors = FALSE)
  ans$small_stratum <- ans$n_children < cls_config(cfg, 'min_clusters', 30L)
  ans$ci_method <- 'child_cluster_percentile_bootstrap_fixed_scores'
  ans$bootstrap_strata <- if(is.null(strata)) 'single_stratum' else 'sex_x_age_group_fixed_counts'
  ans
}
cls_bins <- function(p) {

  z <- rep(NA_integer_, length(p)); good <- is.finite(p) & p >= 0 & p <= 1
  z[good] <- ifelse(p[good] < .03, 1L, ifelse(p[good] < .25, 2L,
                ifelse(p[good] <= .75, 3L, ifelse(p[good] <= .97, 4L, 5L))))
  z
}
cls_bin_labels <- function() c('<P3', 'P3-P25', 'P25-P75', 'P75-P97', '>P97')
cls_score_wide <- function(d, scores, source) {
  ans <- d
  for (m in c('AI','AS_ratio','Conditional_AV')) {
    z <- scores[scores$source == source & scores$metric == m, , drop = FALSE]
    if (anyDuplicated(z$row_id)) stop('Duplicate row_id in ', source, ' scores for ', m)
    match_rows <- match(ans$row_id,z$row_id)
    ans[[paste0('score_present_',m)]] <- !is.na(match_rows)
    p <- z$percentile[match_rows]
    support <- if ('in_support' %in% names(z)) z$in_support[match_rows] else rep(NA,nrow(ans))
    ans[[paste0('support_',m)]] <- support

    if (m == 'Conditional_AV') p[!is.na(support) & !support] <- NA_real_
    ans[[paste0('p_', m)]] <- p
  }
  ans
}
cls_spearman <- function(x,y) {
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x,y, method='spearman'))
}
cls_rank_setup <- function(x) {
  ux <- sort(unique(x)); group <- match(x, ux)
  list(group=group, n=length(ux))
}

cls_weighted_spearman <- function(xsetup, ysetup, w) {
  wx <- rowsum(w, xsetup$group, reorder=TRUE)
  wy <- rowsum(w, ysetup$group, reorder=TRUE)
  rx <- cumsum(wx) - as.numeric(wx)/2
  ry <- cumsum(wy) - as.numeric(wy)/2
  rx <- rx[xsetup$group]; ry <- ry[ysetup$group]
  total <- sum(w)
  if (total < 3) return(NA_real_)
  mx <- sum(w*rx)/total; my <- sum(w*ry)/total
  vx <- sum(w*(rx-mx)^2); vy <- sum(w*(ry-my)^2)
  if (vx <= 0 || vy <= 0) return(NA_real_)
  sum(w*(rx-mx)*(ry-my))/sqrt(vx*vy)
}
cls_agreement <- function(z) {
  z <- matrix(z, 5L, 5L); n <- sum(z)
  if (n <= 0) return(c(exact_agreement=NA_real_, within_one_bin=NA_real_, quadratic_weighted_kappa=NA_real_))
  p <- z/n; delta <- outer(1:5,1:5,'-'); weight <- delta^2/16
  observed <- sum(weight*p); expected <- sum(weight*outer(rowSums(p),colSums(p)))
  c(exact_agreement=sum(diag(p)), within_one_bin=sum(p[abs(delta)<=1]),
    quadratic_weighted_kappa=if (expected > 0) 1-observed/expected else NA_real_)
}
cls_pair <- function(d, xname, yname, source, orientation, cfg, seed) {
  xp <- d[[paste0('p_',xname)]]; yp <- d[[paste0('p_',yname)]]
  if (orientation == 'arch_direction_AI_reversed' && xname == 'AI') xp <- 1-xp
  if (orientation == 'arch_direction_AI_reversed' && yname == 'AI') yp <- 1-yp
  valid <- is.finite(xp) & xp >= 0 & xp <= 1 & is.finite(yp) & yp >= 0 & yp <= 1
  n_candidate <- nrow(d)
  valid_x <- is.finite(xp) & xp >= 0 & xp <= 1
  valid_y <- is.finite(yp) & yp >= 0 & yp <= 1
  coverage <- data.frame(x_metric=xname,y_metric=yname,source=source,
    n_candidate_feet=n_candidate,n_valid_x=sum(valid_x),n_valid_y=sum(valid_y),
    n_missing_x=sum(!valid_x),n_missing_y=sum(!valid_y),n_pairwise_complete_feet=sum(valid),
    n_pairwise_missing_feet=sum(!valid),pairwise_coverage_proportion=mean(valid),
    n_candidate_children=length(unique(d$id)),n_children_with_complete_pair=length(unique(d$id[valid])),
    n_AV_outside_support=if('Conditional_AV' %in% c(xname,yname)) sum(!d$support_Conditional_AV,na.rm=TRUE) else 0L,
    coverage_denominator='all_input_feet_in_stratum',stringsAsFactors=FALSE)
  d <- d[valid,,drop=FALSE]; xp <- xp[valid]; yp <- yp[valid]
  if (!nrow(d)) return(list(cells=data.frame(), agreement=data.frame(),coverage=coverage))
  xb <- cls_bins(xp); yb <- cls_bins(yp); cell <- xb+5L*(yb-1L)
  count <- matrix(0, nrow(d),25L); count[cbind(seq_len(nrow(d)),cell)] <- 1
  cc <- cls_contrib(d$id,count); nchild <- nrow(cc); cn <- colSums(cc)
  child_strata <- cls_child_strata(d$id,cls_strata(d),rownames(cc))
  b <- as.integer(cls_config(cfg,'prop_n_boot',2000L)); chunk <- as.integer(cls_config(cfg,'bootstrap_chunk_size',200L))
  boot_count <- matrix(NA_real_,b,25L); boot_agree <- matrix(NA_real_,b,4L)
  xsetup <- cls_rank_setup(xp); ysetup <- cls_rank_setup(yp)
  child_index <- match(as.character(d$id),rownames(cc))
  set.seed(seed)
  if (nchild >= 2L && b > 0L) for (start in seq.int(1L,b,by=chunk)) {
    ix <- start:min(start+chunk-1L,b)
    w <- cls_boot_weights(child_strata,length(ix))
    bc <- crossprod(w,cc); boot_count[ix,] <- bc
    boot_agree[ix,1:3] <- t(apply(bc,1L,cls_agreement))
    boot_agree[ix,4] <- vapply(seq_along(ix),function(j)
      cls_weighted_spearman(xsetup,ysetup,w[child_index,j]),numeric(1))
  }
  n <- sum(cn); cellmat <- matrix(cn,5L,5L)
  rowden <- rep(rowSums(cellmat),5L); colden <- rep(colSums(cellmat),each=5L)
  boot_total <- rowSums(boot_count)
  boot_row <- matrix(0,b,5L); boot_col <- matrix(0,b,5L)
  for (j in 1:5) {
    boot_row[,j] <- rowSums(boot_count[,seq.int(j,25L,5L),drop=FALSE])
    boot_col[,j] <- rowSums(boot_count[,(5L*(j-1L)+1L):(5L*j),drop=FALSE])
  }
  cells <- data.frame(x_metric=xname,y_metric=yname,source=source,orientation=orientation,
                      x_bin=rep(1:5,5L),y_bin=rep(1:5,each=5L),
                      x_label=rep(cls_bin_labels(),5L),y_label=rep(cls_bin_labels(),each=5L),
                      count=cn,overall_denominator=n,row_denominator=rowden,column_denominator=colden,
                      overall_percent=100*cn/n,row_percent=100*cn/rowden,column_percent=100*cn/colden,
                      n_children=nchild,n_candidate_feet=n_candidate,n_pairwise_complete_feet=n,
                      n_pairwise_missing_feet=n_candidate-n,pairwise_coverage_proportion=n/n_candidate,
                      bootstrap_strata='sex_x_age_group_fixed_counts',
                      small_stratum=nchild < cls_config(cfg,'min_clusters',30L),stringsAsFactors=FALSE)
  cells$row_percent[rowden==0] <- NA_real_; cells$column_percent[colden==0] <- NA_real_
  for (type in c('overall','row','column')) {
    den <- switch(type,overall=rep(n,25L),row=rowden,column=colden)
    bd <- switch(type,overall=matrix(boot_total,b,25L),row=boot_row[,rep(1:5,5L),drop=FALSE],
                 column=boot_col[,rep(1:5,each=5L),drop=FALSE])
    bp <- boot_count/bd; bp[bd==0] <- NA_real_
    ci <- t(vapply(1:25,function(j) cls_ci(bp[,j],cn[j]/den[j],
                       den[j]>0 && (cn[j]==0 || cn[j]==den[j])),numeric(3)))
    cells[[paste0(type,'_lower_percent')]] <- 100*ci[,1]
    cells[[paste0(type,'_upper_percent')]] <- 100*ci[,2]
    cells[[paste0(type,'_n_boot_valid')]] <- ci[,3]
    cells[[paste0(type,'_ci_status')]] <- ifelse(den==0,'no_eligible_observations',
      ifelse(nchild<2,'fewer_than_two_children',ifelse(cn==0|cn==den,'boundary_degenerate',
        ifelse(ci[,3]<20,'insufficient_valid_bootstrap','ok'))))
  }
  est <- c(cls_agreement(cn),spearman_percentile=cls_spearman(xp,yp))
  ci <- t(vapply(1:4,function(j) cls_ci(boot_agree[,j],est[j],
         j<=2 && is.finite(est[j]) && est[j] %in% c(0,1)),numeric(3)))
  agreement <- data.frame(x_metric=xname,y_metric=yname,source=source,orientation=orientation,
      statistic=names(est),estimate=as.numeric(est),lower=ci[,1],upper=ci[,2],
      n_boot_valid=ci[,3],n_feet=n,n_children=nchild,n_candidate_feet=n_candidate,
      n_pairwise_missing_feet=n_candidate-n,pairwise_coverage_proportion=n/n_candidate,
      bootstrap_strata='sex_x_age_group_fixed_counts',
      small_stratum=nchild<cls_config(cfg,'min_clusters',30L),
      ci_status=ifelse(!is.finite(est),'statistic_undefined',
        ifelse(nchild<2,'fewer_than_two_children',
          ifelse(seq_along(est)<=2 & est %in% c(0,1),'boundary_degenerate',
            ifelse(ci[,3]<20,'insufficient_valid_bootstrap','ok')))),
      ci_method='child_cluster_percentile_bootstrap_fixed_scores',stringsAsFactors=FALSE)
  list(cells=cells,agreement=agreement,coverage=coverage)
}

run_classification <- function(d,ref,av,cfg,out) {
  required <- c('row_id','id','sex','age','age_group','side','AI')
  if (!all(required %in% names(d))) stop('Classification data missing: ',paste(setdiff(required,names(d)),collapse=', '))
  if (anyDuplicated(d$row_id)) stop('Classification needs unique row_id.')
  if (anyNA(d$id) || anyNA(d$sex)) stop('Classification requires nonmissing child ID and sex.')
  if (anyDuplicated(paste(d$id,d$side,sep='\r'))) stop('Classification needs at most one row per child and side.')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  scores <- cls_bind(list(ref$scores,av$scores))
  if (!all(c('row_id','metric','percentile','source') %in% names(scores))) stop('Reference score schema is incomplete.')
  invalid_p <- is.finite(scores$percentile) & (scores$percentile<0 | scores$percentile>1)
  if (any(invalid_p)) stop('Reference percentiles must be probabilities between 0 and 1.')
  seed <- as.integer(cls_config(cfg,'seed',20260928L)); tick <- 0L
  next_seed <- function() {tick <<- tick+1L; seed+tick}
  traditional <- foot_positive <- children <- overlap <- cells <- agreements <- assignments <- coverage <- pair_coverage <- list()
  for (g in cls_groups(d,TRUE)) {
    z <- g$data; z <- z[is.finite(z$AI),,drop=FALSE]
    if (nrow(z)) {
      category <- cbind(z$AI<.21,z$AI>=.21 & z$AI<=.26,z$AI>.26)
      tab <- cls_ratios(z$id,category,matrix(1,nrow(z),1),
        c('AI_lt_0.21','AI_0.21_to_0.26_inclusive','AI_gt_0.26'),cfg,next_seed(),cls_strata(z))
      tab$denominator_definition <- 'feet_with_valid_AI'
      traditional[[length(traditional)+1L]] <- cls_meta(tab,g)
      foot_positive[[length(foot_positive)+1L]] <- cls_meta(tab[tab$measure=='AI_gt_0.26',,drop=FALSE],g)
    }
  }
  child_split <- split(seq_len(nrow(d)),as.character(d$id))
  child_data <- cls_bind(lapply(child_split,function(ix) {
    z <- d[ix,,drop=FALSE]
    if ('side_code' %in% names(z)) {
      side_values <- as.character(z$side_code); expected_sides <- c('1','2')
    } else {
      expected_sides <- as.character(cls_config(cfg,'side_labels',c('1'='Side1','2'='Side2')))
      if (length(expected_sides)!=2L || anyNA(expected_sides) || anyDuplicated(expected_sides))
        stop('cfg$side_labels must contain two distinct labels.')
      side_values <- as.character(z$side)
    }
    valid <- all(is.finite(z$AI)) && nrow(z)==2L && setequal(side_values,expected_sides)
    if (length(unique(z$sex)) != 1L || length(unique(z$age_group)) != 1L) stop('Inconsistent sex/age_group within child ',z$id[1])
    data.frame(id=as.character(z$id[1]),sex=as.character(z$sex[1]),age_group=as.character(z$age_group[1]),
               complete_two_AI=valid,n_positive=if(valid) sum(z$AI>.26) else NA_integer_,stringsAsFactors=FALSE)
  }))
  for (g in cls_groups(child_data,TRUE)) {
    z <- g$data[g$data$complete_two_AI,,drop=FALSE]
    if (!nrow(z)) next
    tab <- cls_ratios(z$id,cbind(z$n_positive==0,z$n_positive==1,z$n_positive==2,z$n_positive>=1),
      matrix(1,nrow(z),1),c('neither_foot_AI_gt_0.26','unilateral_AI_gt_0.26',
                          'bilateral_AI_gt_0.26','at_least_one_AI_gt_0.26'),cfg,next_seed(),cls_strata(z))
    tab$denominator_definition <- 'children_with_two_distinct_sides_and_two_valid_AI'
    tab$n_children_total <- nrow(g$data); tab$n_children_ineligible <- nrow(g$data)-nrow(z)
    children[[length(children)+1L]] <- cls_meta(tab,g)
  }
  available_sources <- c('OOF','apparent')
  for (src in available_sources) {
    wide <- cls_score_wide(d,scores,src)
    for (g in cls_groups(wide,TRUE)) for (m in c('AI','AS_ratio','Conditional_AV')) {
      z <- g$data; p <- z[[paste0('p_',m)]]; valid <- is.finite(p) & p>=0 & p<=1
      coverage[[length(coverage)+1L]] <- cls_meta(data.frame(metric=m,source=src,
        n_input_feet=nrow(z),n_valid_percentile=sum(valid),n_missing_percentile=sum(!valid),
        n_absent_score_rows=sum(!z[[paste0('score_present_',m)]]),
        n_outside_support=sum(!z[[paste0('support_',m)]],na.rm=TRUE),
        coverage_proportion=mean(valid),n_input_children=length(unique(z$id)),
        n_children_with_valid_score=length(unique(z$id[valid])),
        coverage_denominator='all_input_feet_in_stratum',stringsAsFactors=FALSE),g)
    }
    a <- wide[,c('row_id','id','sex','age','age_group','side','p_AI','p_AS_ratio','p_Conditional_AV')]
    a$source <- src
    for (m in c('AI','AS_ratio','Conditional_AV')) a[[paste0('bin_',m)]] <- cls_bins(a[[paste0('p_',m)]])
    a$p_AI_arch_direction <- 1-a$p_AI; a$bin_AI_arch_direction <- cls_bins(a$p_AI_arch_direction)
    assignments[[length(assignments)+1L]] <- a
    for (g in cls_groups(wide,TRUE)) {
      z <- g$data; z <- z[is.finite(z$AI)&is.finite(z$p_AI),,drop=FALSE]
      if (!nrow(z)) next
      for (tail in c('lower_AI_tail','upper_AI_tail')) {
        f <- if(tail=='lower_AI_tail') z$AI<.21 else z$AI>.26
        r <- if(tail=='lower_AI_tail') z$p_AI<.25 else z$p_AI>.75
        num <- cbind(f&r,f&!r,!f&r,!f&!r,f,r,f&r,f&r,f&!r,!f&r,xor(f,r))
        den <- cbind(1,1,1,1,1,1,f,r,f,r,1)+matrix(0,nrow(z),11L)
        labels <- c('both_of_all','fixed_only_of_all','reference_only_of_all','neither_of_all',
                    'fixed_positive_of_all','reference_positive_of_all','both_of_fixed_positive',
                    'both_of_reference_positive','fixed_only_of_fixed_positive',
                    'reference_only_of_reference_positive','discordant_of_all')
        tab <- cls_ratios(z$id,num,den,labels,cfg,next_seed(),cls_strata(z))
        tab$source <- src; tab$tail <- tail
        tab$fixed_rule <- if(tail=='lower_AI_tail') 'AI < 0.21' else 'AI > 0.26'
        tab$reference_rule <- if(tail=='lower_AI_tail') 'AI percentile < 0.25' else 'AI percentile > 0.75'
        tab$denominator_definition <- c(rep('feet_with_valid_AI_and_AI_percentile',6),
          'fixed_positive_feet','reference_positive_feet','fixed_positive_feet','reference_positive_feet',
          'feet_with_valid_AI_and_AI_percentile')
        overlap[[length(overlap)+1L]] <- cls_meta(tab,g)
      }
    }
    for(g in cls_groups(wide,FALSE)) for(pair in list(c('AI','AS_ratio'),c('AI','Conditional_AV'),c('AS_ratio','Conditional_AV'))) {
      orientations <- if('AI' %in% pair) c('original_metric_direction','arch_direction_AI_reversed') else 'original_metric_direction'
      pair_seed <- next_seed()
      for(orientation in orientations) {
        z <- cls_pair(g$data,pair[1],pair[2],src,orientation,cfg,pair_seed)
        cells[[length(cells)+1L]] <- cls_meta(z$cells,g)
        agreements[[length(agreements)+1L]] <- cls_meta(z$agreement,g)
        if(orientation=='original_metric_direction') pair_coverage[[length(pair_coverage)+1L]] <- cls_meta(z$coverage,g)
      }
    }
  }
  result <- list(traditional=cls_bind(traditional),foot_AI_gt026=cls_bind(foot_positive),
     children_AI_gt026=cls_bind(children),child_eligibility=child_data,
     fixed_reference_overlap=cls_bind(overlap),percentile_cells=cls_bind(cells),
     agreement=cls_bind(agreements),percentile_assignments=cls_bind(assignments),
     score_coverage=cls_bind(coverage),pair_coverage=cls_bind(pair_coverage),
     fixed_threshold_reference_percentiles=if(is.data.frame(ref$threshold_map)) ref$threshold_map else data.frame())
  filenames <- c(traditional='figure3_AI_traditional_proportions.csv',
    foot_AI_gt026='figure3_AI_gt026_foot_proportions.csv',children_AI_gt026='figure3_AI_gt026_child_unilateral_bilateral.csv',
    child_eligibility='figure3_child_pair_eligibility.csv',fixed_reference_overlap='figure3_fixed_reference_overlap.csv',
    percentile_cells='figure4_percentile_5x5_cells.csv',agreement='figure4_percentile_agreement.csv',
    percentile_assignments='figure4_individual_percentile_assignments.csv',
    score_coverage='classification_score_coverage.csv',pair_coverage='figure4_pairwise_score_coverage.csv',
    fixed_threshold_reference_percentiles='figure3_fixed_threshold_reference_percentiles.csv')
  for(nm in names(filenames)) {
    if(exists('write_table',mode='function')) write_table(result[[nm]],file.path(out,filenames[[nm]]))
    else utils::write.csv(result[[nm]],file.path(out,filenames[[nm]]),row.names=FALSE,na='')
  }
  result$methods <- c(
    'Traditional AI: <0.21, [0.21,0.26], >0.26. Threshold equality belongs to central category.',
    'Fixed/reference overlap uses paired valid AI and AI percentiles; lower <P25; upper >P75.',
    'Child unilateral/bilateral classification requires two distinct sides and both AI values valid.',
    'Five bands: <P3; [P3,P25); [P25,P75]; (P75,P97]; >P97.',
    'AI arch-direction score = 1 - AI percentile; AS_ratio and Conditional_AV retain their percentile direction.',
    'OOF scores are the main internal cross-fitting analysis. Apparent scores are sensitivity analysis.',
    'Bootstrap sampling unit is the child, stratified by sex x age_group with the observed child count fixed in every stratum; all feet and ties are retained. Spearman ranks use the same draw and are recomputed with bootstrap weights.',
    'Out-of-support Conditional_AV scores are excluded. Coverage tables retain all input feet as their denominator, report absent/missing scores and outside-support exclusions. Missing OOF scores are never filled from apparent scores.',
    'Bootstrap intervals condition on precomputed reference/AV scores and omit model-estimation uncertainty.',
    'Boundary-degenerate percentile intervals are reported as NA; counts and denominators remain available.',
    'Small strata are flagged, not discarded. Percentile agreement describes morphology, with no diagnostic interpretation.')
  writeLines(result$methods,file.path(out,'classification_methods.txt'),useBytes=TRUE)
  result
}
