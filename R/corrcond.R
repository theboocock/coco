##' Stepwise conditional wrapper
##'
##' Generate joint signals from a specific set of SNPs
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title conditional_from_ids
##' @param data_set data.frame. Contains the relevant
##' @param ld_matrix

joint_from_ids = function(rsids, res_preparation, idxs,return_only_these=F, joint=T){
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
    if(joint){
      out_tmp = step_joint(betas= data_set$b[rsid_idx],ld_matrix = ld_matrix[rsid_idx,rsid_idx],
                        neffs=data_set$neff[rsid_idx],hwe_diag_outside = hwe_diag_outside[rsid_idx],
                        var_y = var_y,return_entire_beta_set = T,hwe_diag = hwe_diag[rsid_idx])
      to_return = cbind(data_set[rsid_idx,], out_tmp/sqrt(data_set$var[rsid_idx]),2*pnorm(abs(out_tmp[,1]/out_tmp[,2]),lower.tail = F))
      colnamse(to_return) = c(colnames(data_set), "beta_new","se_new","p_new")
      return(to_return)
    }else{
      out_tmp = step_conditional(betas= data_set$b[rsid_idx],ld_matrix = ld_matrix[rsid_idx,rsid_idx],
                           neffs=data_set$neff[rsid_idx],hwe_diag_outside = hwe_diag_outside[rsid_idx],
                           var_y = var_y,return_entire_beta_set = T,hwe_diag = hwe_diag[rsid_idx])
      to_return = cbind(data_set[rsid_idx,], out_tmp/sqrt(data_set$var[rsid_idx]),2*pnorm(abs(out_tmp[,1]/out_tmp[,2]),lower.tail = F))
      colnamse(to_return) = c(colnames(data_set), "beta_new","se_new","p_new")
      return(to_return)
    }
  }else{
    res_step= get_betas(idx_joint = rsid_idx, ld_matrix = ld_matrix,
                  hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside, data_set = data_set,
                  colinear_threshold = 0.9,var_y=var_y,joint=joint)
      res_step = list(res_step=res_step, main_hit=NA,
       conditional_on=c(as.character(data_set$rsid[rsid_idx])))
    return(res_step)
  }
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
stepwise_wrapper = function(data_set ,ld_matrix,p_value_threshold = 1e-3,var_y = 1.6421, ld_noise=0, joint=T){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y)
  stepwise_results = stepwise_conditional_run(res_preparation = res_preparation,p_value_threshold = p_value_threshold, joint=joint)
  all_but_one_df = all_but_one(res_preparation=res_preparation,stepwise_results = stepwise_results, joint=joint)
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
conditional_from_ids_wrapper = function(data_set ,ld_matrix, rsids=NULL, idxs=NULL,var_y = NULL, ld_noise=0, exact=F,hwe_variance=T,joint=T){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y,hwe_variance = hwe_variance,exact = exact)
  stepwise_results = joint_from_ids(rsids=rsids, idxs=idxs, res_preparation = res_preparation, joint=joint)
  return (stepwise_results)
}




##' Step joint
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

step_joint = function(betas, ld_matrix,neffs,var_y, hwe_diag_outside,hwe_diag,return_entire_beta_set=F){
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
  #print(beta_two_inv)
  #print(betas)
  new_betas = beta_inv %*% diag(hwe_diag) %*% betas
  #print(beta_inv)
  # We only really care about the results from the first SNP
  #because that's our SNP we are adding to the model
  vars = (var_y * (median(neffs))- t(new_betas) %*%  diag(hwe_diag) %*% ((betas))) / (median(neffs) - n_betas )
  ses = sqrt(diag(vars[1] * beta_inv))
  if(return_entire_beta_set){
    return(cbind(new_betas[,1],ses))  
  }else{
  return(c(new_betas[1,1],ses[1]))
  }
}



##' Step joint
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
  # todo add the new math.
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
  beta_two = outside[1,1]
  beta_one = outside[2:ncol(outside),2:ncol(outside)]
  beta_one_inv = chol2inv(chol(beta_one))
  beta_two_inv = chol2inv(chol(beta_two))
  beta_two_cond_lhs = betas[1]
  x2_x1 = outside[1,2:ncol(outside)]
  if((length(betas) - 1) == 1 ){
    diag_one = hwe_diag[2:length(hwe_diag)]
  }else{
    diag_one = diag(hwe_diag[2:length(hwe_diag)])
  }
  beta_two_cond_rhs = beta_two_inv %*% x2_x1 %*% beta_one_inv %*% diag_one %*% betas[2:length(betas)]
  new_beta = beta_two_cond_lhs - beta_two_cond_rhs
  #print(new_beta)
  #print(diag_one)
  #print(beta_one_inv)
  #print(beta_two_inv)
  new_beta_joint = beta_one_inv %*% diag_one %*% betas[2:length(betas)]
  vars = (var_y * (median(neffs))- t(new_beta_joint) %*%  diag_one %*% ((betas[2:length(betas)])) - t(new_beta) %*% beta_two %*% new_beta) / (median(neffs) - n_betas )
  ses = sqrt(diag(vars[1] *beta_two_inv))
  if(return_entire_beta_set){
    return(cbind(new_beta[,1],ses))  
  }else{
    return(c(new_beta[1,1],ses[1]))
  }
}


get_betas = function(idx_joint,data_set,ld_matrix, hwe_diag,hwe_diag_outside,var_y,colinear_threshold=.9, joint=T){
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
    if(joint){
      cond_results = step_joint(betas, tmp_ld_matrix, neffs,var_y, hwe_diag_outside_tmp,hwe_diag_tmp)
    }else{
      cond_results = step_conditional(betas, tmp_ld_matrix, neffs,var_y, hwe_diag_outside_tmp,hwe_diag_tmp)
    }
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

all_but_one = function(res_preparation, stepwise_results,p_value_threshold=1e-6,colinear_threshold=0.9,joint=T){
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
    res_step= get_betas(idx_joint = tmp_model_idx, ld_matrix = ld_matrix,
                        hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside,data_set = data_set,
                        colinear_threshold = 0.9,var_y = var_y, joint=joint)
    hit = which(!(1:length(idx_joint) %in% combinations[,i]))
    all_but_one_res[[i]] = list(res_step=res_step, main_hit=as.character(data_set$rsid[idx_joint[hit]]),
                                conditional_on=c(as.character(data_set$rsid[idx_joint[combinations[,i]]])))
  }
  return(all_but_one_res)
}

