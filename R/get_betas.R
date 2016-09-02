##' Internal function get betas.
##'
##' Run the all but one analysis of the dataset.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title get_betas.
##' @param data_set input file format
##' @param ld_matrix LD matrix for the subset of the region.
##' @param stepwise_results - results from a stepwise analysis.
##' @return all_but_ones  betas and standard error for the entire region.
##'
##' @export
get_betas = function(idx_joint,data_set,ld_matrix, hwe_diag,hwe_diag_outside,var_y,colinear_threshold=.9, joint=T, exact=F){
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
      cond_results = step_joint(betas, tmp_ld_matrix, neffs,var_y, hwe_diag_outside_tmp,hwe_diag_tmp, exact=exact)
    }else{
      cond_results = step_conditional(betas, tmp_ld_matrix, neffs,var_y, hwe_diag_outside_tmp,hwe_diag_tmp, exact=exact)
    }
    res_step$beta_new[j]= cond_results[1] 
    res_step$se_new[j]= cond_results[2] 
    res_step$z_new[j]= res_step$beta_new[j]/res_step$se_new[j]
    res_step$p_new[j] = 2 * pnorm(abs(res_step$z_new[j]), lower.tail = F)
  }
  message(paste("Skipped",n_snps_skipped,"SNPs as colinear"))
  return(res_step)
}
