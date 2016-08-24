##' corrocond - cojo module of coco.
##'
##' Coco's functions for implementing the GCTA-cojo in R.
##' 
##' Currently only includes the joint modelling approach. 
##' 
##' TODO: implement the conditional selection approach in the GCTA paper.
##' 
##' 
##' 
##' @author James Boocock
##' @date Aug 19th 2017
##' @citation Yang, Jian, et al. "Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants 
##'           influencing complex traits." Nature genetics (2012)



##' Stepwise conditional wrapper
##'
##' Generate conditional signals from a specific set of SNPs
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title conditional_from_ids
##' @param data_set data.frame. Contains the relevant
##' @param ld_matrix

conditional_from_ids = function(rsids, res_preparation, idxs,return_only_these=F){
  hwe_diag = res_preparation$hwe_diag
  data_set = res_preparation$data_set
  hwe_diag_outside = res_preparation$hwe_diag_outside
  ld_matrix = res_preparation$ld_matrix
  var_y = res_preparation$var_y
  if(!is.null(rsids)){
    rsid_idx = which(data_set$rsid %in% rsids)
  } else if(!is.null(idxs)){
    rsid_idx = idxs
  }else{
    stop("Neither indexes of rsids specified")
  }
  if(length(rsid_idx) == 0){
    stop("No rsids or indexes to conditional on were found")
  }
  if(return_only_these){
      out_tmp = step_conditional(betas= data_set$b[rsid_idx],ld_matrix = ld_matrix[rsid_idx,rsid_idx],
                        neffs=data_set$neff[rsid_idx],hwe_diag_outside = hwe_diag_outside[rsid_idx],
                        var_y = var_y,return_entire_beta_set = T,hwe_diag = hwe_diag[rsid_idx])
      to_return = cbind(data_set[rsid_idx,], out_tmp/sqrt(data_set$var[rsid_idx]),2*pnorm(abs(out_tmp[,1]/out_tmp[,2]),lower.tail = F))
      colnamse(to_return) = c(colnames(data_set), "beta_new","se_new","p_new")
      return(to_return)
  }else{
    res_step= get_joint_betas(idx_joint = rsid_idx, ld_matrix = ld_matrix,
                  hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside, data_set = data_set,
                  colinear_threshold = 0.9,var_y=var_y)
      res_step = list(res_step=res_step, main_hit=NA,
       conditional_on=c(as.character(data_set$rsid[rsid_idx])))
    return(res_step)
  }
}
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
  return(median(vars)/(data_set$n-1))
}


##'
##' Generate conditional signals from a conditional dataset.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title stepwise_conditional
##' @param data_set data.frame. Contains the relevant
##' @param ld_matrix
##' @param p_value_threshold.
##' @param ld_noise
##'
##'
stepwise_conditional_wrapper = function(data_set ,ld_matrix,p_value_threshold = 1e-3,var_y = 1.6421, ld_noise=0){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y)
  stepwise_results = stepwise_conditional_run(res_preparation = res_preparation,p_value_threshold = p_value_threshold)
  all_but_one_df = all_but_one(res_preparation=res_preparation,stepwise_results = stepwise_results)
  return (all_but_one_df)
}

##'
##' Generate conditional signals from a conditional dataset.
##'
##'
##' @author Eli A. Stahl
##'  @data 2 Aug 2016
##' @title stepwise_conditional
##' @param data_set data.frame. Contains the relevant
##' @param ld_matrix
##' @param p_value_threshold.
##' @param ld_noise
##'
##'
conditional_from_ids_wrapper = function(data_set ,ld_matrix, rsids=NULL, idxs=NULL,var_y = 1, ld_noise=0){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y)
  stepwise_results = conditional_from_ids(rsids=rsids, idxs=idxs, res_preparation = res_preparation)
  return (stepwise_results)
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
prep_dataset_common = function(data_set,ld_matrix,var_y, ld_noise=0,hwe_variance=F){
    if(missing(data_set)){
      stop("Must specify a data_set argument")
    }
    if(missing(ld_matrix)){
      stop("Must specify an ld_matrix argument on the command-line")
    }

    message("Preparing the dataset for a  conditional analysis")
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

    info = which(grepl("INFO|RSQ", names(data_set),  ignore.case = T))
    if(length(info) == 0){
      message("Info column not found assuming HWE genotypic variance")
      data_set$info = 1
    }
    colnames(data_set)[info] = "info"

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
      data_set$var = 2*data_set$af*(1-data_set$af)*data_set$info
    }else{
      print(var)
      data_set$var = data_set[,var]
    }
    if(missing(var_y)){
      message("Var_y missing, will estimate using cojo method")
      var_y = estimate_vary(data_set)
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
    hwe_diag =  (data_set$neff - 1) * data_set$var
   # idxs = which(!(data_set$neff > (mean(data_set$neff) + 6 * sd(data_set$neff)) | data_set$neff < (mean(data_set$neff) - 6 * sd(data_set$neff))))
  #  ld_matrix = ld_matrix[idxs,idxs]
  #  data_set = data_set[idxs,]
  #hwe_diag_outside =  (2*data_set$af * ( 1- data_set$af) * data_set$info)
    hwe_diag_outside = data_set$var * (data_set$neff -1)
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
stepwise_conditional_run = function(res_preparation,p_value_threshold=0.0001,colinear_threshold=0.9){
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

##' Step conditional
##'
##' Run a joint step of the analysis.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title step_conditional
##' @param betas vector of betas that will be adjusted in the matrix operations.
##' @param ld_matrix LD matrix for the subset of the region.
##' @param neff effective sample size for each SNP in the dataset
##' @param var_y variance of of the phenotype
##' @param hwe_diag_outside Diagonal matrix containing the genotypic variance
##' @param hwe_diag Diagonal matrix containing the genotypic variance *n
##' @return Joint betas and standard errors.

step_conditional = function(betas, ld_matrix,neffs,var_y, hwe_diag_outside,hwe_diag,return_entire_beta_set=F){
  inside =  ld_matrix
  n_betas = length(betas)
  outside = sqrt(diag(hwe_diag_outside)) %*% inside %*% sqrt(diag(hwe_diag_outside))
  #print(outside)
  for(j in 1:ncol(outside)){
    for(k in (j):ncol(outside)){
      if(k > ncol(outside)){
        break
      }
      outside[j,k] = min(neffs[c(j,k)]) *  outside[j,k] 
      if(k!=j){
        outside[k,j] =  min(neffs[c(j,k)]) * outside[k,j] 
      }
    }
  }
  #print( sqrt(diag(hwe_diag)) %*% inside %*% sqrt(diag(hwe_diag)))
#  beta_inv = chol2inv(chol(outside))
  outside = sqrt(diag(hwe_diag_outside)) %*% inside %*% sqrt(diag(hwe_diag_outside))
  #print(outside)
  
  beta_inv = chol2inv(chol(outside))
  beta_inv= (solve(outside))
  #print(betas)
  new_betas = beta_inv %*% diag(hwe_diag) %*% betas
  #print(beta_inv)
  # We only really care about the results from the first SNP
  #because that's our SNP we are adding to the model
  vars = (var_y * (median(neffs) -1)- t(new_betas) %*%  diag(hwe_diag) %*% ((betas))) / (median(neffs) - n_betas - 1)
  #print(vars)
  ses = sqrt(diag(vars[1] * beta_inv))
  #print(ses)
  if(return_entire_beta_set){
    return(cbind(new_betas[,1],ses))  
  }else{
  return(c(new_betas[1,1],ses[1]))
  }
}


get_conditional_betas = function(idx_joint,data_set,ld_matrix, hwe_diag,hwe_diag_outside,var_y,colinear_threshold=.9){

}


get_joint_betas = function(idx_joint,data_set,ld_matrix, hwe_diag,hwe_diag_outside,var_y,colinear_threshold=.9){
  res_step = data.frame(rsid=data_set$rsid,pos=data_set$pos,beta_old=data_set$b,se_old=data_set$se, z_old=data_set$z, p_old=data_set$p,
                        beta_new=rep(NA,nrow(ld_matrix)),
                        se_new=rep(NA,nrow(ld_matrix)), z_new=rep(NA,nrow(ld_matrix)),p_new=rep(NA,nrow(ld_matrix)))
  n_snps_skipped = 0
  for(j in 1:nrow(data_set)){
    if(j %in% idx_joint){
#        message("Skipping SNP as already being conditioned on")
        next
    }

    of_interest = c(j, idx_joint)
    betas = data_set$b[of_interest]
    tmp_ld_matrix = ld_matrix[of_interest,of_interest]
    neffs = data_set$neff[of_interest]
    if(any(abs(tmp_ld_matrix[1,2:ncol(tmp_ld_matrix)]) > colinear_threshold)){
      # message("Skipping SNP as co-linear with top_snps.")
      n_snps_skipped = n_snps_skipped + 1
      next
    }
    
    hwe_diag_outside_tmp = hwe_diag_outside[of_interest]
    hwe_diag_tmp = hwe_diag[of_interest]
    cond_results = step_conditional(betas, tmp_ld_matrix, neffs,var_y, hwe_diag_outside_tmp,hwe_diag_tmp)

    res_step$beta_new[j]= cond_results[1] 
    res_step$se_new[j]= cond_results[2] 
    res_step$z_new[j]= res_step$beta_new[j]/res_step$se_new[j]
    res_step$p_new[j] = 2 * pnorm(abs(res_step$beta_new[j]/res_step$se_new[j]), lower.tail = F)
  }
  message(paste("Skipped",n_snps_skipped,"SNPs as colinear"))
  return(res_step)
}

##' All but one analysis
##'
##' Run the all but one analysis of the dataset.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title all_but_one
##' @param data_set input file format
##' @param ld_matrix LD matrix for the subset of the region.
##' @param stepwise_results - results from a stepwise analysis.
##' @return all_but_ones  betas and standard error for the entire region.
##'
all_but_one = function(res_preparation, stepwise_results,p_value_threshold=1e-6,colinear_threshold=0.9){
  #extract variance from the dataset
  hwe_diag = res_preparation$hwe_diag
  data_set = res_preparation$data_set
  var_y = res_preparation$var_y
  hwe_diag_outside = res_preparation$hwe_diag_outside
  ld_matrix = res_preparation$ld_matrix
  stepwise_results = stepwise_results[stepwise_results$p_new < p_value_threshold,]
  if(nrow(stepwise_results) == 0){
    message("No conditional SNPs remaining after apply p value threshold.")
  }
  idx_joint = match(stepwise_results$rsid,data_set$rsid)
  idx_joint = rev(idx_joint[!is.na(idx_joint)])

  if(length(idx_joint) == 0){
    stop("No SNPs to perform an all but one conditional analysis on")
  }
  idxs = 1:length(idx_joint)
  if(length(idx_joint) == 1){
   combinations=matrix(c(1))
  }else{
    combinations = combn(length(idx_joint),length(idx_joint) - 1)
    print(combinations)
  }
  message(paste("Runnnig the all but one conditional analysis for ",length(idx_joint)," snps"))
  #message(paste("SNPs in model ... "))
  #message(as.character(data_set$rsid[idx_joint]), sep="\n")
  all_but_one_res = list()
  for(i in 1:ncol(combinations)){
  #  idx_conditional = c(tmp_model_idx)
    tmp_model_idx = idx_joint[combinations[,i]]
    res_step= get_joint_betas(idx_joint = tmp_model_idx, ld_matrix = ld_matrix,
                              hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside,data_set = data_set,
                              colinear_threshold = 0.9,var_y = var_y)
    hit = which(!(1:length(idx_joint) %in% combinations[,i]))
    all_but_one_res[[i]] = list(res_step=res_step, main_hit=as.character(data_set$rsid[idx_joint[hit]]),
                                             conditional_on=c(as.character(data_set$rsid[idx_joint[combinations[,i]]])))
  }
  return(all_but_one_res)
}
