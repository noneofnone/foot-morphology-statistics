pipeline_run <- function(cfg,code_dir) {
  check_packages()
  if(!cfg$stage%in%c("check","quick","full")) stop("stage must be check, quick or full")
  if(cfg$stage=="quick") {
    cfg$n_boot<-min(cfg$n_boot,20L);cfg$av_n_boot<-min(cfg$av_n_boot,20L)
    cfg$prop_n_boot<-min(cfg$prop_n_boot,200L);cfg$n_sim<-min(cfg$n_sim,1000L)
    cfg$av_age_n_sim<-min(cfg$av_age_n_sim %||% 5000L,1000L)
  }
  cfg$international_dir<-cfg$international_dir %||% file.path(code_dir,"inputs")
  if(!file.exists(cfg$data_file))stop("Raw CSV not found: ",cfg$data_file,"\nEdit 00_config.R.")
  files<-c(file.path(code_dir,"00_config.R"),list.files(file.path(code_dir,"R"),pattern="\\.R$",full.names=TRUE))
  package_versions<-vapply(required_packages(),function(p)as.character(packageVersion(p)),character(1))
  international_cfg<-cfg[["international",exact=TRUE]]
  resolved_international_dir<-if(!is.null(international_cfg$input_dir))international_cfg$input_dir else cfg$international_dir
  international_files<-file.path(resolved_international_dir,c("international_series_REVIEW.csv","international_values_extracted.csv"))
  international_md5<-if(isTRUE(cfg$run_international))unname(tools::md5sum(international_files))else NULL
  hash_config<-cfg;hash_config$resume<-NULL
  key<-object_hash(list(data=unname(tools::md5sum(cfg$data_file)),code=unname(tools::md5sum(files)),
    config=hash_config,international=international_md5,packages=package_versions,R=as.character(getRversion())))
  out<-file.path(cfg$output_dir,paste0("run_",cfg$stage,"_",substr(key,1,12)))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  writeLines(normalizePath(out,winslash="/"),file.path(cfg$output_dir,"LATEST_RUN.txt"))
  for(x in c("qc","models","results","checkpoints"))dir.create(file.path(out,x),showWarnings=FALSE)
  dput(cfg,file.path(out,"configuration_used.R"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  stage_log<-list()
  run_stage<-function(name,fun,key) {
    file<-file.path(out,"checkpoints",paste0(name,".rds"))
    t0<-proc.time()[3]; cached<-FALSE
    if(isTRUE(cfg$resume)&&file.exists(file)) {
      old<-tryCatch(readRDS(file),error=function(e)NULL)
      if(!is.null(old)&&identical(old$key,key)) { value<-old$value;cached<-TRUE }
    }
    if(!cached) {
      if (isTRUE(cfg$verbose)) message("Running ",name," ...")
      value<-fun();atomic_rds(list(key=key,value=value),file)
    } else if (isTRUE(cfg$verbose)) message("Resuming cached ",name)
    stage_log[[name]]<<-data.frame(stage=name,elapsed_seconds=proc.time()[3]-t0,cached=cached)
    write_table(bind_rows_base(stage_log),file.path(out,"stage_timings.csv"));value
  }
  writeLines("RUNNING",file.path(out,"RUN_STATUS.txt"))
  failed<-TRUE
  on.exit(if(failed)writeLines("FAILED: inspect the R error and saved stage diagnostics. Earlier completed stages are retained.",file.path(out,"RUN_STATUS.txt")),add=TRUE)
  d<-read_input(cfg,out)
  write_table(data.frame(file=c(basename(cfg$data_file),basename(files)),md5=unname(tools::md5sum(c(cfg$data_file,files)))),file.path(out,"INPUT_CODE_MANIFEST.csv"))
  if(cfg$stage=="check") {failed<-FALSE;writeLines("INPUT_CHECK_COMPLETED",file.path(out,"RUN_STATUS.txt"));return(invisible(list(data=d)))}
  writeLines(key,file.path(out,"checkpoints","run_fingerprint.txt"))
  age<-run_stage("age",function()run_age(d,cfg,out),key)
  ref<-run_stage("reference",function()run_reference(d,cfg,out),key)
  av_reference<-run_stage("av_reference",function()run_av(d,cfg,out),key)
  av_gamma<-run_stage("av_gamma",function()run_av_gamma(d,av_reference,cfg,out),key)
  av<-combine_av_models(av_reference,av_gamma,cfg,out)
  cl<-run_stage("classification",function()run_classification(d,ref,av,cfg,file.path(out,"results/classification")),key)
  ass<-if(isTRUE(cfg$run_associations))run_stage("associations",function()run_associations(d,cfg,file.path(out,"results/associations")),key)else NULL
  st<-run_stage("stability",function()run_stability(d,age,av,cfg,out),key)
  intl<-if(isTRUE(cfg$run_international))run_stage("international",function()run_international(d,ref,age,cfg,out),key)else NULL
  qc<-list()
  if(!is.null(ref$status)&&nrow(ref$status))for(j in seq_len(nrow(ref$status))) {
    rr<-ref$status[j,];qc[[length(qc)+1L]]<-data.frame(component=paste("reference",rr$metric,rr$sex,sep="_"),
      check="fit_OOF_and_bootstrap",status=if(isTRUE(rr$model_success)&&isTRUE(rr$OOF_complete)&&isTRUE(rr$CI_available))"pass"else"review",
      detail=as.character(rr$error %||% ""))
  }
  if(!is.null(ass))qc[[length(qc)+1L]]<-data.frame(component="associations",check="Wald_models",
    status=if(all(ass$tests$status=="ok"))"pass"else"review",detail=paste("failed",sum(ass$tests$status!="ok")))
  if(!is.null(av$progress_qc))qc[[length(qc)+1L]]<-data.frame(component="AV_progress",check="endpoint_and_monotonicity",
    status="review",detail="Inspect av_progress_qc.csv; unreliable normalized progress is suppressed")
  if(!is.null(intl))qc[[length(qc)+1L]]<-data.frame(component="international",check="external_comparability",
    status=as.character(intl$status),detail="See international/00_status.csv and source-review metadata")
  write_table(bind_rows_base(qc),file.path(out,"QC_SUMMARY.csv"))
  bundle<-list(version=cfg$version,stage=cfg$stage,data=d,age=age,reference=ref,av=av,classification=cl,
    associations=ass,stability=st,international=intl,config=cfg,fingerprint=key,output_dir=out)
  atomic_rds(bundle,file.path(out,"STATISTICAL_RESULTS.rds"))
  plot_manifest<-data.frame(
    figure=c("1","2","2","2","2","3","4","4","5"),
    content=c("FL AI AAR sex-specific reference centiles","FL AI AAR age gradients","AV partial age gradient and normalized integral",
      "Observed-span progress and gradient heatmap","Stable interval histograms","Traditional AI and age-sex references",
      "OOF 5-by-5 percentile cross-classification","BMI and hallux adjusted morphology curves","External comparison"),
    location=c("results/reference/reference_centiles.csv","results/age/gradients.csv","results/av/av_gradients.csv; results/av/av_progress.csv",
      "results/age/all_four_progress.csv; results/age/all_four_heatmap.csv","models/stable_interval_histogram_data.rds",
      "results/classification/figure3_*.csv","results/classification/figure4_percentile_*.csv",
      "results/associations/figure4_primary_six_curves.csv","international/05_figure5_external_comparison_interface.csv"))
  write_table(plot_manifest,file.path(out,"PLOTTING_INPUT_MANIFEST.csv"))
  paths<-list.files(out,recursive=TRUE,full.names=TRUE);paths<-paths[file.info(paths)$isdir%in%FALSE]
  write_table(data.frame(file=substring(paths,nchar(out)+2),bytes=file.info(paths)$size),file.path(out,"OUTPUT_MANIFEST.csv"))

  notes<-c(if(cfg$stage=="full")"FULL_EXECUTION_COMPLETED_REVIEW_MODEL_QC"else"QUICK_EXECUTION_COMPLETED_TEST_OUTPUT_ONLY")
  writeLines(notes,file.path(out,"RUN_STATUS.txt"));failed<-FALSE
  if (isTRUE(cfg$verbose)) message("Statistical outputs: ",normalizePath(out,winslash="/"))
  invisible(bundle)
}
