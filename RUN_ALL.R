local({
  frames<-sys.frames();files<-Filter(function(x)!is.null(x)&&length(x)==1L,lapply(frames,function(x)x$ofile))
  cli<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  self<-if(length(files))tail(files,1)[[1]]else if(length(cli))sub("^--file=","",cli[1])else "RUN_ALL.R"
  root<-dirname(normalizePath(self,mustWork=TRUE))
  e<-new.env(parent=globalenv())
  sys.source(file.path(root,"00_config.R"),envir=e)
  for(p in sort(list.files(file.path(root,"R"),pattern="\\.R$",full.names=TRUE)))sys.source(p,envir=e)
  e$pipeline_run(e$foot_config(root),root)
})
