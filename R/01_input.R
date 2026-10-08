read_input <- function(cfg,out) {
  if(!file.exists(cfg$data_file)) stop("Raw CSV not found: ",cfg$data_file,"\nEdit data_file in 00_config.R.")
  raw <- read.csv(cfg$data_file,check.names=FALSE,colClasses="character",fileEncoding=cfg$encoding,
    stringsAsFactors=FALSE,na.strings=c("","NA","N/A","NULL","NaN"))
  names(raw) <- trimws(sub("^\ufeff","",names(raw)))
  if(anyDuplicated(names(raw))) stop("Duplicated column names")
  need <- c("participant_id","gender","age","height","weight","BMI","side","FL","AI","hallux","AS","AV","S","AS/S")
  if(length(setdiff(need,names(raw)))) stop("Missing columns: ",paste(setdiff(need,names(raw)),collapse=", "))

  raw <- raw[,need,drop=FALSE]; raw$participant_id <- trimws(raw$participant_id)
  if(anyNA(raw$participant_id)||any(!nzchar(raw$participant_id))) stop("Missing participant_id")
  d <- data.frame(row_id=seq_len(nrow(raw)),id=sprintf("A%06d",match(raw$participant_id,sort(unique(raw$participant_id)))))
  bad <- list()
  for (n in setdiff(need,"participant_id")) {
    s <- trimws(raw[[n]]); x <- suppressWarnings(as.numeric(s)); ii <- which(!is.na(s)&(!is.finite(x)))
    if(length(ii)) bad[[n]] <- data.frame(row_id=ii,variable=n,raw_value=s[ii])
    d[[n]] <- x
  }
  write_table(bind_rows_base(bad),file.path(out,"qc/numeric_parse_failures.csv"))
  if(length(bad)) stop("Non-numeric input: inspect qc/numeric_parse_failures.csv")
  if(anyNA(d[,c("gender","age","side")])) stop("Missing gender, age or side")
  if(any(!d$gender%in%c(1,2))||any(!d$side%in%c(1,2))) stop("Expected gender 1=boys,2=girls and recorded side 1/2")
  if(any(d$age<cfg$age_boundary[1]|d$age>=cfg$age_boundary[2])) stop("Age outside configured interval")
  groups <- split(seq_len(nrow(d)),d$id)
  pair <- lapply(groups,function(ii) data.frame(id=d$id[ii[1]],n=length(ii),side1=sum(d$side[ii]==1),side2=sum(d$side[ii]==2)))
  pair <- bind_rows_base(pair); write_table(pair,file.path(out,"qc/participant_pairs.csv"))
  if(any(pair$n!=2|pair$side1!=1|pair$side2!=1)) stop("Require exactly one record per recorded side per participant")
  conflicts <- list()
  for(n in c("gender","age","height","weight","BMI")) {
    bad_id <- names(Filter(function(ii) {
      x <- d[[n]][ii]; (anyNA(x)&&!all(is.na(x))) || (sum(is.finite(x))>1&&diff(range(x,na.rm=TRUE))>1e-6)
    },groups))
    if(length(bad_id)) conflicts[[n]] <- data.frame(id=bad_id,variable=n)
  }
  write_table(bind_rows_base(conflicts),file.path(out,"qc/demographic_conflicts.csv"))
  if(length(conflicts)) stop("Within-participant demographic inconsistencies")
  nchild <- length(groups)
  write_table(data.frame(n_children=nchild,n_feet=nrow(d)),file.path(out,"qc/sample_size.csv"))
  if(isTRUE(cfg$require_expected_n)&&(nchild!=cfg$expected_children||nrow(d)!=cfg$expected_feet)) stop("Audited sample-size mismatch")
  invalid <- list()
  for(n in c("height","weight","FL","S","AS","AV","AI")) {
    x <- d[[n]]; ii <- which(is.finite(x)&if(n%in%c("AS","AV","AI")) x<0 else x<=0)
    if(n=="AI") ii <- union(ii,which(x>1))
    if(length(ii)) invalid[[n]] <- data.frame(row_id=ii,id=d$id[ii],variable=n,value=x[ii])
  }
  write_table(bind_rows_base(invalid),file.path(out,"qc/invalid_values.csv"))
  if(length(invalid)) stop("Mathematically invalid measurements; inspect qc/invalid_values.csv")
  d$BMI_file <- d$BMI; d$AS_ratio_file <- d[["AS/S"]]; d[["AS/S"]] <- NULL
  d$BMI <- d$weight/(d$height/100)^2
  d$AS_ratio <- d$AS/d$S; d$AV_FL3 <- d$AV/d$FL^3; d$AS_FL2 <- d$AS/d$FL^2
  check <- data.frame(row_id=d$row_id,BMI_file=d$BMI_file,BMI_recomputed=d$BMI,
    AS_ratio_file=d$AS_ratio_file,AS_ratio_recomputed=d$AS_ratio,
    BMI_mismatch=abs(d$BMI-d$BMI_file)>cfg$bmi_tolerance,
    AS_ratio_mismatch=abs(d$AS_ratio-d$AS_ratio_file)>cfg$ratio_tolerance)
  write_table(check,file.path(out,"qc/derived_variable_agreement.csv"))
  if(!isTRUE(cfg$allow_derived_discrepancies)&&any(check$BMI_mismatch|check$AS_ratio_mismatch,na.rm=TRUE))
    stop("Recomputed BMI or AS/S differs from file; inspect derived_variable_agreement.csv")
  d$sex <- factor(d$gender,levels=c(1,2),labels=c("Boys","Girls"))
  d$side_code <- d$side; d$side <- factor(as.character(d$side),levels=c("1","2"),labels=cfg$side_labels[c("1","2")])
  d$age_group <- floor(d$age+1e-8); d$gender <- NULL
  d <- d[order(d$id,d$side),,drop=FALSE]; rownames(d) <- NULL
  count <- aggregate(row_id~sex+age_group,d,length); names(count)[3]<-"n_feet"; count$n_children<-count$n_feet/2
  write_table(count,file.path(out,"qc/sample_by_age_sex.csv"))
  miss <- lapply(c("FL","AI","AS","AV","S","AS_ratio","hallux","BMI"),function(n) data.frame(
    variable=n,n_missing=sum(!is.finite(d[[n]])),n_zero=sum(d[[n]]==0,na.rm=TRUE)))
  write_table(bind_rows_base(miss),file.path(out,"qc/missingness_and_zeros.csv"))

  limits <- list(height=c(60,230),weight=c(5,250),BMI=c(7,70),FL=c(80,400),hallux=c(-40,80),AS_ratio=c(0,1))
  flags <- lapply(names(limits),function(n){ii<-which(d[[n]]<limits[[n]][1]|d[[n]]>limits[[n]][2]);
    data.frame(row_id=d$row_id[ii],variable=rep(n,length(ii)),value=d[[n]][ii])})
  write_table(bind_rows_base(flags),file.path(out,"qc/plausibility_review.csv"))
  atomic_rds(d,file.path(out,"models/analysis_data.rds"))
  d
}
