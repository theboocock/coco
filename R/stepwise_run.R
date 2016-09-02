#' Forward model selection in coco.
#'
#' \code{stepwise_coco}. Forward model selection in coco.
#'
#' @author James Boocock
#' @date 2 Aug 2016
#' @title get_betas.
#' @param data_set input file format
#' @param ld_matrix LD matrix for the subset of the region.
#' @param stepwise_results - results from a stepwise analysis.
#' @return all_but_ones  betas and standard error for the entire region.
#'
#' @export
stepwise_coco = function(res_preparation,p_value_threshold=0.0001,colinear_threshold=0.9,joint=T,max_iter=100, return_all_betas=F){
  #extract variance from the dataset
  #hwe_diag = (2*freq_af$af * ( 1- data_set$af) * data_set$n)
  # Extract the effective sample size for a SNP.
  
  if(class(res_preparation) != "coco_data"){
    stop("res_preparation object must be of type coco_data")
  }
  hwe_diag = res_preparation$hwe_diag
  data_set = res_preparation$data_set
  var_y = res_preparation$var_y
  hwe_diag_outside = res_preparation$hwe_diag_outside
  ld_matrix = res_preparation$ld_matrix
  var_y = res_preparation$var_y
  exact = res_preparation$exact
  print(exact)
  #print(nrow(ld_matrix))
  idx_top_tmp = c(which(max(abs(data_set$z)) == abs(data_set$z)))
  conditional_df = data_set[idx_top_tmp,]
  #  print(idx_top_tmp)
  res_betas = list()
  message(paste("Top SNP for region = ", conditional_df$rsid, ", with Z-score ", conditional_df$z), " index ",idx_top_tmp)
  current_best_p = 2 * pnorm(abs(conditional_df$z), lower.tail = F)
  idx_cond = c()
  message(paste("Conditioning ..."))
  out_all_buts = data.frame(rsid=conditional_df$rsid, pos=conditional_df$pos, beta_old=conditional_df$b,se_old=conditional_df$se,
                            z_old=conditional_df$z, p_old=current_best_p, beta_new=conditional_df$b,se_new=conditional_df$se,z_new=conditional_df$z,p_new=current_best_p)
  rownames(out_all_buts) = c(idx_top_tmp)
  iters = 1
  last_p = 0
  while(current_best_p < p_value_threshold && iters < max_iter ){
    idx_cond = c(idx_cond, idx_top_tmp)
    
    if(length(idx_cond) == nrow(ld_matrix)){
      message("All SNPs in joint model")
      break
    }
    res_step= get_betas(idx_joint = idx_cond, ld_matrix = ld_matrix,
                              hwe_diag = hwe_diag, hwe_diag_outside = hwe_diag_outside,data_set = data_set,
                              colinear_threshold = colinear_threshold, var_y=var_y,exact=exact,joint=joint)
    res_step$z_new = res_step$beta_new/res_step$se_new
    idx_top_tmp = which(max(abs(res_step$z_new),na.rm=T) == abs(res_step$z_new))
    if(length(idx_top_tmp) > 1){
      idx_top_tmp = idx_top_tmp[1]
    }
    best_cond_row = res_step[idx_top_tmp,]
    current_best_p = best_cond_row$p_new
    
    message(paste("Best SNP is = ", best_cond_row$rsid, "with P-value ", current_best_p, " Index = ", idx_top_tmp))
    message(paste("Beta = ", best_cond_row$beta_new, " Z = ",best_cond_row$z_new))
    if(return_all_betas){
      res_betas[[iters]] = list(res_step=res_step, main_hit=NA,conditional_on=as.character(data_set$rsid)[idx_cond])
    }
    if(current_best_p < p_value_threshold){
      # best_cond_row$p = current_best_p
      out_all_buts = rbind(out_all_buts,best_cond_row)
    }
    iters = iters + 1
  }
  if(return_all_betas){
    return(structure(list(stepwise_summary=out_all_buts, step_betas=res_betas),class="coco_stepwise"))
  }
  return(structure(list(stepwise_summary=out_all_buts), class="coco_stepwise"))
}


