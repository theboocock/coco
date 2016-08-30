##'
##' Estimate the variance in Y using the cojo method.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title estimate_vary
##' @param data_set data.frame. Contains the relevant
##' @param data_set
##'
##' Note: Below equation 8 in the cojo paper.
##'
##' Estimate trait variance.
##'

estimate_vary = function(data_set){
  vars=data_set$var *(data_set$n) * data_set$se ^2 * (data_set$n -1)  +  data_set$var * (data_set$n) * data_set$b^2
  print(paste("Var y = ",median(vars/(data_set$n -1))))
  # return(median(vars/(data_set$n-1)))
  #print(vars)
  return(median(vars/(data_set$n-1)))
}
##' Prepare dataset common
##'
##' Common utilities for checking and parsing a dataset
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title stepwise_conditional
##' @param data_set
##'
prep_dataset_common = function(data_set,ld_matrix,var_y, ld_noise=0,hwe_variance=F,exact=F){
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
  af = which(grepl("^F$|freq|FRQ|MAF", names(data_set),  ignore.case = T))[1]
  if(length(af) == 0){
    stop("AF column not found")
  }# sometimes have 2 (one for each allele), doesn't matter whcih you take for our applications (fGWAS and coloc)
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
    data_set$var = 2*data_set$af*(1-data_set$af)
  }else{
    var_flag_hwe = F
    data_set$var = data_set[,var]
  }
  if(missing(var_y)){
    message("Var_y missing, will estimate using cojo method")
    var_y = estimate_vary(data_set)
  }
  info = which(grepl("INFO|RSQ", names(data_set),  ignore.case = T))
  if(length(info) == 0){
    if(var_flag_hwe){
      message("Info column not found assuming HWE genotypic variance")
      data_set$var = data_set$var 
    }
  }else{
    if(var_flag_hwe){
      colnames(data_set)[info] = "info"
      print(data_set$var)
      message("Info column not found using combination with HWE to estimate genotypic variance")
      data_set$var = data_set$var * data_set$info
    }
  }
  if(exact){
    # Check for the columns needed by excat. 
    gvar_present = !is.na(var)
    if(!gvar_present){
      stop("Can not perform an exact analysis if GVAR is missing.")
    }
    message("   Performing an exact analysis, assuming the genotypic variance has been calculated
            from the hard-called or dossage-based LD matrix that was used to generate the association statistics. \n
            GVAR column is assumed to have been calculated directly from the data.")
  }
  #var_y = 4.52020627198678
  #order the dataset on the position column
  data_set = data_set[order(data_set$pos),]
  # data_set$b = data_set$b * sqrt(data_set$var)
  #  data_set$se = data_set$se * sqrt(data_set$var)
  #hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$info * data_set$n)
  hwe_diag =  ( data_set$n - 1) * data_set$var
  # hwe_diag =  ( data_set$n -1)
  #data_set$neff = (var_y) / (hwe_diag *data_set$se^2)  - (data_set$b^2) / (data_set$se^2) +1
  data_set$neff = (var_y* (data_set$n-1)) / (hwe_diag *data_set$se^2)  - (data_set$b^2) / (data_set$se^2) +1
  #data_set$neff = (median(var_y*hwe_diag)) / (hwe_diag *data_set$se^2)  - (data_set$b^2) / (data_set$se^2) +1
  #hwe_diag =  data_set$n
  if(exact){
    #data_set$n = data_set$n -1
    hwe_diag =  (data_set$n -1 ) * data_set$var
    # idxs = which(!(data_set$neff > (mean(data_set$neff) + 6 * sd(data_set$neff)) | data_set$neff < (mean(data_set$neff) - 6 * sd(data_set$neff))))
    #  ld_matrix = ld_matrix[idxs,idxs]
    #  data_set = data_set[idxs,]
    #hwe_diag_outside =  (2*data_set$af * ( 1- data_set$af) * data_set$info)
    hwe_diag_outside = data_set$var * (data_set$n -1 )
    data_set$neff = data_set$n -1
  }else{
    hwe_diag =  (data_set$neff) * data_set$var
    # idxs = which(!(data_set$neff > (mean(data_set$neff) + 6 * sd(data_set$neff)) | data_set$neff < (mean(data_set$neff) - 6 * sd(data_set$neff))))
    #  ld_matrix = ld_matrix[idxs,idxs]
    #  data_set = data_set[idxs,]
    #hwe_diag_outside =  (2*data_set$af * ( 1- data_set$af) * data_set$info)
    hwe_diag_outside = data_set$var * (data_set$neff)
  }
  #print(data_set$neff)
  if(ld_noise != 0){
    ld_matrix[cbind(1:nrow(ld_matrix),1:nrow(ld_matrix))]  = ld_matrix[cbind(1:nrow(ld_matrix),1:nrow(ld_matrix))] + ld_noise
  }
  return(list(data_set=data_set,ld_matrix=ld_matrix,var_y=var_y,hwe_diag=hwe_diag,hwe_diag_outside=hwe_diag_outside))
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
stepwise_run = function(res_preparation,p_value_threshold=0.0001,colinear_threshold=0.9,joint=T){
  #extract variance from the dataset
  #hwe_diag = (2*freq_af$af * ( 1- data_set$af) * data_set$n)
  # Extract the effective sample size for a SNP.
  hwe_diag = res_preparation$hwe_diag
  data_set = res_preparation$data_set
  var_y = res_preparation$var_y
  hwe_diag_outside = res_preparation$hwe_diag_outside
  ld_matrix = res_preparation$ld_matrix
  var_y = res_preparation$var_y
  #print(nrow(ld_matrix))
  idx_top_tmp = c(which(max(abs(data_set$z)) == abs(data_set$z)))
  conditional_df = data_set[idx_top_tmp,]
  #  print(idx_top_tmp)
  message(paste("Top SNP for region = ", conditional_df$rsid, ", with Z-score ", conditional_df$z), " index ",idx_top_tmp)
  current_best_p = 2 * pnorm(abs(conditional_df$z), lower.tail = F)
  idx_cond = c()
  message(paste("Conditioning ..."))
  out_all_buts = data.frame(rsid=conditional_df$rsid, pos=conditional_df$pos, beta_old=conditional_df$b,se_old=conditional_df$se,
                            z_old=conditional_df$z, p_old=current_best_p, beta_new=conditional_df$b,se_new=conditional_df$se,z_new=conditional_df$z,p_new=current_best_p)
  rownames(out_all_buts) = c(idx_top_tmp)
  while(current_best_p < p_value_threshold){
    idx_cond = c(idx_cond, idx_top_tmp)
    
    if(length(idx_cond) == nrow(ld_matrix)){
      message("All SNPs in joint model")
      break
    }
    res_step= get_joint_betas(idx_joint = idx_cond, ld_matrix = ld_matrix,
                              hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside,data_set = data_set,
                              colinear_threshold = colinear_threshold, var_y=var_y)
    res_step$z_new = res_step$beta_new/res_step$se_new
    idx_top_tmp = which(max(abs(res_step$z_new),na.rm=T) == abs(res_step$z_new))
    if(length(idx_top_tmp) > 1){
      idx_top_tmp = idx_top_tmp[1]
    }
    best_cond_row = res_step[idx_top_tmp,]
    current_best_p = best_cond_row$p_new
    
    message(paste("Best SNP is = ", best_cond_row$rsid, "with P-value ", current_best_p, " Index = ", idx_top_tmp))
    message(paste("Beta = ", best_cond_row$beta_new, " Z = ",best_cond_row$z_new))
    if(current_best_p < p_value_threshold){
      # best_cond_row$p = current_best_p
      out_all_buts = rbind(out_all_buts,best_cond_row)
    }
  }
  return(out_all_buts)
}


