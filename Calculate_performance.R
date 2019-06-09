CalcPerf <- function(ref, type, excl_ct, info, excl = NULL) {
  if(!is.null(excl)) ref <- ref[!ref %in% excl]; type <- type[names(ref)]
  # % of no.ref unclassified
  norefcells <- type[ref %in% excl_ct]
  noref <- sum(grepl('Node|Unassigned', norefcells))
  noref_cl <- sum(!grepl('Node|Unassigned', norefcells))
  #% of yes.ref correct
  yesref <- type[!(ref %in% excl_ct)]
  yrref <- ref[names(yesref)]
  corr <- sum(yesref == yrref)
  #% of yes.ref incorrect{
  incorr <- sum(yesref != yrref & !grepl("Node|Unassigned", yesref))
  ## % of yes.ref in Node1 (new Node0 == Unassigned)
  yesrefi <- yesref[grepl('Node|Unassigned', yesref)]
  yrref <- ref[names(yesrefi)]
  incunass <- sum(grepl('Unassigned', yesrefi))
  ## % of yes.ref in correct branch
  perc_gb <- 0 ## percentage in good branch,
  for (i in unique(yrref)) {
    Nodes <- paste0('Node', which(unlist(lapply(info$nodetypes, function (x) i %in% names(x)))) - 1)
    Nodes <- Nodes[!grepl("Node0", Nodes)]
    perc_gb <- (perc_gb + sum(yesrefi[yrref == i] %in% Nodes))
  }
  perc_gb <- perc_gb
  ## % of yes.ref in incorrect branch
  perc_wb <- length(yesref) - incunass - perc_gb - corr - incorr
  out <- c(perc_wb + perc_gb, incunass, corr, incorr, noref, noref_cl)/length(ref)
  names(out) <- c('intermediate', 'unassigned', 'correct', 
                  'incorrect', 'correctly unclassified', 'erroneously classified')
  return(out)
}