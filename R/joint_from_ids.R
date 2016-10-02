#' Joint from ids conditional wrapper
#'
#' \code{joint_from_idxs} Performs a joint or conditional signals from a specific set of SNPs, or indexes.
#'
#' @author James Boocock and Eli Stahl
#' @date 2 Aug 2016
#' @title conditional_from_ids
#' @param res_prepration coco_data - the result of running the coco prep_dataset_common command.
#' @param rsids character vector - SNPIDs in res preparation to perform CoCo analysis on.
#' @param idxs integer vector - postiion in input data to condition on.
#' @param return_only_these boolean
#' \itemize{
#' \item (T) = return the joint or conditional statistics for only SNPs in idxs or rsids vectors.
#' \item (F) = return the joint or conditional statistics for SNPs not in the idxs or rsids vectors.
#' }
#'                                  
#' @param joint boolean
#' \itemize{
#'  \item(T) = perform joint estimation.
#'  \item(F) = perform conditional estimation
#'  }      
#' 
#' @examples 
#' 
#' data(coco_nek1)
#' res = joint_from_ids(rsids="rs10520157",coco_nek1,exact=T)
#' # Second hit rs4417927
#' res$res_step[res$res_step$rsid == "rs4417927",]
#' 
#' @export
joint_from_ids = function(res_preparation, rsids=NULL, idxs=NULL,return_only_these=F, joint=T,exact=F){
  if(class(res_preparation) != "coco_data"){
    stop("res_preparation object must be of type coco_data")
  }
  
  hwe_diag = res_preparation$hwe_diag
  data_set = res_preparation$data_set
  hwe_diag_outside = res_preparation$hwe_diag_outside
  ld_matrix = res_preparation$ld_matrix
  var_y = res_preparation$var_y
  exact = res_preparation$exact
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
                           var_y = var_y,return_entire_beta_set = T,hwe_diag = hwe_diag[rsid_idx], exact=exact)
      to_return = cbind(data_set[rsid_idx,], out_tmp/sqrt(data_set$var[rsid_idx]),2*pnorm(abs(out_tmp[,1]/out_tmp[,2]),lower.tail = F))
      colnames(to_return) = c(colnames(data_set), "beta_new","se_new","p_new")
      return(to_return)
    }else{
      out_tmp = step_conditional(betas= data_set$b[rsid_idx],ld_matrix = ld_matrix[rsid_idx,rsid_idx],
                                 neffs=data_set$neff[rsid_idx],hwe_diag_outside = hwe_diag_outside[rsid_idx],
                                 var_y = var_y,return_entire_beta_set = T,hwe_diag = hwe_diag[rsid_idx], exact=exact)
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