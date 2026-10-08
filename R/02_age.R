age_basis <- function(d,cfg) {
  b <- splines::ns(d$age,knots=cfg$age_knots,Boundary.knots=cfg$age_boundary)
  for(j in seq_len(ncol(b))) d[[paste0("a",j)]] <- b[,j]
  d
}
age_mean_design <- function(obj,age,sex,cfg) {
  X <- lapply(levels(obj$data$side),function(side) {
    z <- data.frame(age=age,sex=factor(rep(sex,length(age)),levels=c("Boys","Girls")),
      side=factor(rep(side,length(age)),levels=levels(obj$data$side)))
    design(obj,age_basis(z,cfg))
  })
  Reduce(`+`,X)/length(X)
}
progress_summary <- function(age,curve,draws=NULL) {
  span <- tail(curve,1)-curve[1]
  tol <- 1e-6*max(abs(curve),1e-6)
  valid <- is.finite(span)&&abs(span)>tol
  den <- if(!is.null(draws)) draws[nrow(draws),]-draws[1,] else NULL
  span_lo <- span_hi <- NA_real_
  if(length(den)) {
    ci <- quantile(den,c(.025,.975),na.rm=TRUE); span_lo<-unname(ci[1]);span_hi<-unname(ci[2])
    valid <- valid && !(span_lo<=0&&span_hi>=0)
  }
  est <- if(valid) 100*(curve-curve[1])/span else rep(NA_real_,length(age))
  lo<-hi<-rep(NA_real_,length(age))
  if(valid&&!is.null(draws)) {
    ok <- is.finite(den)&abs(den)>tol&sign(den)==sign(span)
    if(mean(ok)>=.95) {
      pp <- 100*sweep(sweep(draws[,ok,drop=FALSE],2,draws[1,ok],"-"),2,den[ok],"/")
      lo <- apply(pp,1,quantile,.025,na.rm=TRUE); hi<-apply(pp,1,quantile,.975,na.rm=TRUE)
    }
  }
  nonmonotone <- if(valid) any(diff(est)< -1e-5) else NA
  data.frame(age,estimate=est,lower=lo,upper=hi,endpoint_span=span,span_lower=span_lo,span_upper=span_hi,
    progress_defined=valid,nonmonotone=nonmonotone,
    interpretation="Percent of fitted change across observed reference age span; endpoint 100 is imposed by normalization")
}
run_age <- function(d,cfg,out) {
  z<-age_basis(d,cfg); rhs<-paste0("a",seq_len(length(cfg$age_knots)+1L),collapse="+")
  curves<-gradients<-progress<-heatmap<-models<-ledger<-list()
  for(metric in c("FL","AI","AS_ratio")) {
    f<-as.formula(paste0(metric,"~(",rhs,")*sex+side"))
    start<-proc.time()[3]; obj<-fit_gee(z,f,cfg); models[[metric]]<-obj
    allX<-allD<-meta<-list()
    for(sx in c("Boys","Girls")) {
      range_age<-range(obj$data$age[obj$data$sex==sx]); grid<-cfg$age_grid[cfg$age_grid>=range_age[1]&cfg$age_grid<=range_age[2]]
      if(length(grid)<3L) stop("Insufficient age support for ",metric," ",sx)
      lo<-pmax(grid-.005,range_age[1]);hi<-pmin(grid+.005,range_age[2])
      allX[[sx]]<-age_mean_design(obj,grid,sx,cfg)
      allD[[sx]]<-(age_mean_design(obj,hi,sx,cfg)-age_mean_design(obj,lo,sx,cfg))/(hi-lo)
      meta[[sx]]<-data.frame(metric,sex=sx,age=grid)
    }
    X<-do.call(rbind,allX); D<-do.call(rbind,allD); m<-bind_rows_base(meta)
    cb<-linear_bands(obj$beta,obj$V,X,cfg$n_sim,cfg$seed+match(metric,c("FL","AI","AS_ratio")))
    gb<-linear_bands(obj$beta,obj$V,D,cfg$n_sim,cfg$seed+100+match(metric,c("FL","AI","AS_ratio")))
    curves[[metric]]<-cbind(m,cb,n_children=obj$n_children,n_feet=obj$n_feet)
    gradients[[metric]]<-cbind(m,gb,units=switch(metric,FL="mm/year",AI="AI/year",AS_ratio="AS/S/year"),
      band_family="all supported ages and both sexes within metric",estimator="GEE population-mean cross-sectional age gradient")
    draws<-normal_draws(obj$beta,obj$V,cfg$n_sim,cfg$seed+200)
    for(sx in names(allX)) {
      ii<-which(m$sex==sx); yy<-drop(allX[[sx]]%*%obj$beta)
      pp<-progress_summary(m$age[ii],yy,allX[[sx]]%*%draws)
      progress[[paste(metric,sx)]]<-cbind(metric,sex=sx,pp)
      heatmap[[paste(metric,sx)]]<-data.frame(metric,sex=sx,age=m$age[ii],
        normalized_gradient=if(pp$progress_defined[1])gb$estimate[ii]/pp$endpoint_span[1]*100 else NA_real_,
        definition="percentage of observed-span fitted change per year",raw_gradient=gb$estimate[ii],
        raw_sim_lower=gb$sim_lower[ii],raw_sim_upper=gb$sim_upper[ii])
    }
    ledger[[metric]]<-data.frame(metric,formula=paste(deparse(f),collapse=""),n_children=obj$n_children,
      n_feet=obj$n_feet,elapsed_seconds=proc.time()[3]-start,working_correlation=obj$corstr,
      alpha=as.numeric(obj$alpha)[1],status="fitted")
  }
  ans<-list(curves=bind_rows_base(curves),gradients=bind_rows_base(gradients),progress=bind_rows_base(progress),
    heatmap=bind_rows_base(heatmap),models=models,ledger=bind_rows_base(ledger))
  for(n in c("curves","gradients","progress","heatmap","ledger")) write_table(ans[[n]],file.path(out,"results/age",paste0(n,".csv")))
  atomic_rds(models,file.path(out,"models/age_models.rds"));ans
}
sustained_onset <- function(g,min_tail=1,margin=NA_real_) {
  g<-g[order(g$age),,drop=FALSE]
  if(!nrow(g)) return(data.frame(onset=NA_real_,status="no_supported_ages",tail_years=NA_real_))
  usable<-is.finite(g$sim_lower)&is.finite(g$sim_upper)
  ok<-usable&g$sim_lower<=0&g$sim_upper>=0
  if(is.finite(margin)) ok<-usable&g$sim_lower>= -margin&g$sim_upper<=margin
  sustained<-rev(cumprod(rev(as.integer(ok))))==1
  ii<-which(sustained&g$age<=max(g$age)-min_tail)
  data.frame(onset=if(length(ii))g$age[ii[1]]else NA_real_,
    status=if(length(ii))if(is.finite(margin))"within_prespecified_equivalence_margin"else"sustained_nondetectable_change" else "not_identified_in_observed_range",
    tail_years=if(length(ii))max(g$age)-g$age[ii[1]]else NA_real_)
}
run_stability <- function(d,age,av,cfg,out) {
  av_g<-av$gradients
  if(!is.null(av_g)&&nrow(av_g)) {

    if(!"metric"%in%names(av_g)) av_g$metric<-"AV_sizeadjusted"
    if("unit"%in%names(av_g))av_g$units<-av_g$unit
  }
  g<-bind_rows_base(list(age$gradients,av_g)); ans<-hist<-list()
  for(key in unique(paste(g$metric,g$sex,sep="::"))) {
    ii<-paste(g$metric,g$sex,sep="::")==key; zz<-g[ii,,drop=FALSE]
    metric<-zz$metric[1];sx<-zz$sex[1]
    r<-sustained_onset(zz,cfg$min_tail_years)
    margin<-unname(cfg$equivalence_margin[metric]);if(!length(margin))margin<-NA_real_
    eq<-sustained_onset(zz,cfg$min_tail_years,margin)
    a<-data.frame(metric,sex=sx,r,equivalence_margin=margin,
      equivalence_onset=if(is.finite(margin))eq$onset else NA_real_,
      equivalence_status=if(is.finite(margin))eq$status else "margin_not_prespecified")
    ans[[key]]<-a
    if(metric=="AV_sizeadjusted"&&!is.null(av$gamma)) {
      hh<-av$gamma$histogram_data
      if(is.data.frame(hh)&&nrow(hh))hist[[key]]<-hh[hh$sex==sx,,drop=FALSE]
      next
    }
    if(is.finite(r$onset)) {

      dd<-d[d$sex==sx&d$age>=r$onset&d$age<=max(zz$age),,drop=FALSE]
      variables<-if(metric=="AV_sizeadjusted")c("AV","AV_FL3")else metric
      for(v in variables) if(v%in%names(dd)) hist[[paste(key,v)]]<-data.frame(
        row_id=dd$row_id,id=dd$id,sex=as.character(dd$sex),age=dd$age,selection_metric=metric,
        distribution_variable=v,value=dd[[v]],onset=r$onset,
        interpretation="Foot-level distribution in a data-selected non-detectable-change interval")
      if(metric=="AV_sizeadjusted"&&!is.null(av$scores)&&nrow(av$scores)) {
        ss<-av$scores[av$scores$source=="OOF"&av$scores$row_id%in%dd$row_id,,drop=FALSE]
        if(nrow(ss))hist[[paste(key,"percentile")]]<-data.frame(row_id=ss$row_id,id=ss$id,sex=ss$sex,age=ss$age,
          selection_metric=metric,distribution_variable="Conditional_AV_percentile_OOF",value=ss$percentile,
          onset=r$onset,interpretation="Conditional morphology percentile in a data-selected interval")
      }
      if(metric=="AV_sizeadjusted"&&nrow(dd)&&!is.null(av$models)) {
        refFL<-median(dd$FL,na.rm=TRUE)
        adjusted<-av_adjust_to_FL(dd,av$models[sx],setNames(refFL,sx))
        adjusted<-adjusted[adjusted$reference_FL_in_support,,drop=FALSE]
        if(nrow(adjusted))hist[[paste(key,"size_adjusted_volume")]]<-data.frame(
          row_id=adjusted$row_id,id=adjusted$id,sex=adjusted$sex,age=adjusted$age,
          selection_metric=metric,distribution_variable="AV_at_stable_window_reference_FL",
          value=adjusted$AV_at_reference_FL,onset=r$onset,reference_FL=adjusted$reference_FL,
          interpretation="mm^3; each observed AV adjusted to median FL in the selected stable-age window; reference FL within local support")
      }
    }
  }
  result<-list(stability=bind_rows_base(ans),histogram_data=bind_rows_base(hist),gradients=g)
  result$progress<-bind_rows_base(list(age$progress,av$progress))
  result$progress$display_eligible<-ifelse(result$progress$metric=="AV_sizeadjusted",
    result$progress$progress_reliable,result$progress$progress_defined)
  av_heat<-list()
  if(!is.null(av$gamma)) {
    av_heat<-list(av$gamma$heatmap)
  } else if(!is.null(av$progress)&&nrow(av$progress)&&!is.null(av_g)&&nrow(av_g))for(sx in unique(av_g$sex)) {
    gg<-av_g[av_g$sex==sx,,drop=FALSE]; pp<-av$progress[av$progress$sex==sx,,drop=FALSE]
    if(nrow(pp)&&"cumulative_log_change"%in%names(pp)) {
      denominator<-tail(pp$cumulative_log_change,1)
      reliable<-if("progress_reliable"%in%names(pp))all(pp$progress_reliable)else is.finite(denominator)&&abs(denominator)>1e-8
      av_heat[[sx]]<-data.frame(metric="AV_sizeadjusted",sex=sx,age=gg$age,
        normalized_gradient=if(reliable)gg$estimate/denominator*100 else NA_real_,
        definition="percentage of integrated partial-logAV age change per year",raw_gradient=gg$estimate,
        raw_sim_lower=gg$sim_lower,raw_sim_upper=gg$sim_upper)
    }
  }
  result$heatmap<-bind_rows_base(c(list(age$heatmap),av_heat))
  write_table(result$stability,file.path(out,"results/age/stability_onsets.csv"))
  write_table(result$gradients,file.path(out,"results/age/all_four_gradients.csv"))
  write_table(result$progress,file.path(out,"results/age/all_four_progress.csv"))
  write_table(result$heatmap,file.path(out,"results/age/all_four_heatmap.csv"))
  atomic_rds(result$histogram_data,file.path(out,"models/stable_interval_histogram_data.rds"))
  if(isTRUE(cfg$export_individual))write_table(result$histogram_data,file.path(out,"results/age/stable_interval_histogram_data.csv"))
  result
}
