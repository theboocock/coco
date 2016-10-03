#' Prepare dataset coco.
#'
#' \code{prep_dataset_coco} Common processing function for checking and parsing a dataset.
#' @title estimate_vary
#' @author James Boocock
#' @date 2 Aug 2016
#' @title stepwise_conditional
#' @param data_set data.frame. Coco formatted input data.
#' @param ld_matrix matrix. LD matrix ordered by postiion.
#' @param var_y float. Phenotypic variance for Y, if not specified it will be estimated using \code{\link{estimate_vary}}.
#' @param hwe_variance boolean. Force the use of Hardy-weinberg equillibrium when calculating the genotypic variance.
#' @param exact boolean. Perform an exact analysis, extra data checks used in \code{\link{prep_data_set_common}}
#' 
#' @examples 
#' data(coco_dataset)
#' data(coco_ld_matrix)
#' coco_data = prep_datset_coco(coco_dataset,coco_ld_matrix)
#' coco_data$var_y
#' 
#' @export
prep_dataset_coco = function(data_set,ld_matrix,var_y,hwe_variance=F,exact=F,use_info=T){
  if(sum(class(data_set) %in% "data.table") > 0){
      data_set = setDF(data_set)
  }
  if(missing(data_set)){
    stop("Must specify a data_set argument")
  }
  if(missing(ld_matrix)){
    stop("Must specify an ld_matrix argument on the command-line")
  }
  message("Preparing the dataset for coco")
  id = which(grepl("^SNP$|^rsid$", names(data_set), ignore.case = T))
  if(length(id) == 0){
    stop("RSID column not found")
  }
  colnames(data_set)[id] = "rsid"
  pos = which(grepl("^pos$|^bp|^bp_hg19$|^position$|^chrom_start$", names(data_set),  ignore.case = T))
  if(length(pos) == 0){
    stop("Position column not found")
  }
  colnames(data_set)[pos] = "pos"
  beta = which(grepl("^beta$|^effect$|^logor$|^b$", names(data_set),  ignore.case = T))[1]
  if(is.na(beta)){
    beta =  which(grepl("^or$", names(data_set),  ignore.case = T))[1]
    if(is.na(beta)){
      stop("BETA/OR column not found")
    }
    data_set[,beta] = log(data_set[,beta])
  }
  colnames(data_set)[beta] = "b"
  af = which(grepl("^FREQ1$|^F$|^freq$|^FRQ$|^MAF$|^af$", names(data_set),  ignore.case = T))[1]
  # sometimes have 2 (one for each allele), doesn't matter whcih you take for our applications (fGWAS and coloc)
  if(is.na(af)){
    stop("AF column not found")
  }
  colnames(data_set)[af] = "af"
  se = which(grepl("^se$|^StdErr$|BMIadjMainSE|^SE$", names(data_set),  ignore.case = F))
  if(length(se) == 0){
    stop("SE column not found")
  }
  #print(colnames(data_set[se]))
  colnames(data_set)[se] = "se"
  n = which(grepl("^N$|^n$|^sum_n$|^TotalSampleSize$|^SAMPLE_SIZE$", names(data_set),  ignore.case = T))
  if(length(n) == 0){
    stop("n column not found")
  }
  colnames(data_set)[n] = "n"
  z = which(grepl("^Z$|zscore", names(data_set),  ignore.case = T))
  if(length(z) ==0){
    data_set$z = data_set$b/data_set$se
  } else {
    colnames(data_set)[z] = "z"
  }
  
  p = which(grepl("^P$|^Pval", names(data_set),  ignore.case = T))
  if(length(p) ==0){
    data_set$p = 2*pnorm(abs(data_set$z),lower.tail=F)
  } else {
    colnames(data_set)[p] = "p"
  }
  var = which(grepl("^GVAR$", names(data_set), ignore.case= T))[1]
  if(is.na((var))|| hwe_variance){
    var_flag_hwe = T
    message("GVAR column not found or hwe_variance flag set. Assuming HWE genotypic variance")
    data_set$var = 2*data_set$af*(1-data_set$af)
  }else{
    var_flag_hwe = F
    data_set$var = data_set[,var]
  }
  if(missing(var_y)){
    message("Var_y missing, will estimate using cojo method")
    var_y = estimate_vary(data_set)
  }
  if(is.na(var) & use_info){
    info = which(grepl("INFO|RSQ", names(data_set),  ignore.case = T))
    if(is.na(info)){
      if(var_flag_hwe){
        message("Info column not found. Assuming HWE genotypic variance")
      }
    }else{
      if(var_flag_hwe){
        colnames(data_set)[info] = "info"
#        print(data_set$var)
        message("Info column found. Adjusting HWE var to estimate genotypic variance")
        data_set$var = data_set$var * data_set$info
      }
    }
  }

  if(exact){
    # Check for the columns needed by excat. 
    gvar_present = !is.na(var)
    if(!gvar_present){
      stop("Can not perform an exact analysis if GVAR is missing.")
    }
    message("Performing an exact analysis, assuming the genotypic variance has been calculated
            from the hard-called or dossage-based LD matrix that was used to generate the association statistics. \n
            GVAR column is assumed to have been calculated directly from the data.")
  }
  #var_y = 4.52020627198678
  #order the dataset on the position column
  data_set = data_set[order(data_set$pos),]
  # data_set$b = data_set$b * sqrt(data_set$var)
  #  data_set$se = data_set$se * sqrt(data_set$var)
  #hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$info * data_set$n)
  
 # (_jma_Vp - _MSX[i] * _beta[i] * _beta[i]) / (_MSX[i] * _beta_se[i] * _beta_se[i]
  data_set$neff = (var_y - data_set$var *data_set$b^2)/(data_set$var *data_set$se^2) + 1
  #hwe_diag =  ( data_set$n - 1) * data_set$var
  # hwe_diag =  ( data_set$n -1)
  #data_set$neff = (var_y) / (hwe_diag *data_set$se^2)  - (data_set$b^2) / (data_set$se^2) +1
  #data_set$neff = (var_y* (data_set$n-1)) / (hwe_diag *data_set$se^2)  - (data_set$b^2) / (data_set$se^2) +1
  #data_set$neff = (var_y * (data_set$n-1)- hwe_diag*data_set$b ^2)/ (hwe_diag*data_set$se^2) + 1
  
  #data_set$neff = (median(var_y*hwe_diag)) / (hwe_diag *data_set$se^2)  - (data_set$b^2) / (data_set$se^2) +1
  #hwe_diag =  data_set$n
  if(exact){
    #data_set$n = data_set$n -1
    hwe_diag =  (data_set$neff -1 ) * data_set$var
    # idxs = which(!(data_set$neff > (mean(data_set$neff) + 6 * sd(data_set$neff)) | data_set$neff < (mean(data_set$neff) - 6 * sd(data_set$neff))))
    #  ld_matrix = ld_matrix[idxs,idxs]
    #  data_set = data_set[idxs,]
    #hwe_diag_outside =  (2*data_set$af * ( 1- data_set$af) * data_set$info)
    hwe_diag_outside = data_set$var * (data_set$neff -1 )
    data_set$neff = data_set$neff -1
  }else{
    hwe_diag =  (data_set$neff) * data_set$var
    # idxs = which(!(data_set$neff > (mean(data_set$neff) + 6 * sd(data_set$neff)) | data_set$neff < (mean(data_set$neff) - 6 * sd(data_set$neff))))
    #  ld_matrix = ld_matrix[idxs,idxs]
    #  data_set = data_set[idxs,]
    #hwe_diag_outside =  (2*data_set$af * ( 1- data_set$af) * data_set$info)
    hwe_diag_outside = data_set$var 
  }
  #print(data_set$neff)
  return(structure(list(data_set=data_set,ld_matrix=ld_matrix,var_y=var_y,hwe_diag=hwe_diag,hwe_diag_outside=hwe_diag_outside,exact=exact),class="coco_data"))
}

##' Stepwise estimation.
##'
##' Generate conditional signals from a conditional dataset.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title stepwise_conditional
##' @param data_set
##' @param p_value_threshold.
##'
##'


