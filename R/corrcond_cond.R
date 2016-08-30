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
    res_step= get_betas(idx_joint = rsid_idx, ld_matrix = ld_matrix,
                              hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside, data_set = data_set,
                              colinear_threshold = 0.9,var_y=var_y, joint=F)
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
stepwise_conditional_wrapper = function(data_set ,ld_matrix,p_value_threshold = 1e-3,var_y = 1.6421, ld_noise=0, exact=F){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y, exact=exact)
  stepwise_results = stepwise_run(res_preparation = res_preparation,p_value_threshold = p_value_threshold, joint=F)
  all_but_one_df = conditional_all_but_one(res_preparation=res_preparation,stepwise_results = stepwise_results)
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
conditional_from_ids_wrapper = function(data_set ,ld_matrix, rsids=NULL, idxs=NULL,var_y = NULL, ld_noise=0, exact=F,hwe_variance=T){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y,hwe_variance = hwe_variance,exact = exact)
  stepwise_results = conditional_from_ids(rsids=rsids, idxs=idxs, res_preparation = res_preparation, joint=F)
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
  
  beta_inv = chol2inv(chol(outside))
  beta_inv= (solve(outside))
  #print(betas)
  new_betas = beta_inv %*% diag(hwe_diag) %*% betas
  #print(beta_inv)
  # We only really care about the results from the first SNP
  #because that's our SNP we are adding to the model
  
  vars = (var_y * (median(neffs) -1)- t(new_betas) %*%  diag(hwe_diag) %*% ((betas))) / (median(neffs) - n_betas -1)
  ses = sqrt(diag(vars[1] * beta_inv))
  
  if(return_entire_beta_set){
    return(cbind(new_betas[,1],ses))  
  }else{
    return(c(new_betas[1,1],ses[1]))
  }
}




