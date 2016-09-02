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

all_but_one = function(res_preparation, stepwise_results,p_value_threshold=1e-6,colinear_threshold=0.9,joint=T, exact=F){
  #extract variance from the dataset
  if(class(res_preparation) != "coco_data"){
    stop("res_preparation object must be of type coco_data")
  }
  if(class(stepwise_results) != "coco_stepwise"){
    stop("stepwise results object must be of type coco_stepwise")
  }
  hwe_diag = res_preparation$hwe_diag
  data_set = res_preparation$data_set
  var_y = res_preparation$var_y
  hwe_diag_outside = res_preparation$hwe_diag_outside
  ld_matrix = res_preparation$ld_matrix
  stepwise_results = stepwise_results$stepwise_summary[stepwise_results$stepwise_summary$p_new < p_value_threshold,]
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
                        colinear_threshold = colinear_threshold,var_y = var_y, joint=joint, exact=exact)
    hit = which(!(1:length(idx_joint) %in% combinations[,i]))
    all_but_one_res[[i]] = list(res_step=res_step, main_hit=as.character(data_set$rsid[idx_joint[hit]]),
                                conditional_on=c(as.character(data_set$rsid[idx_joint[combinations[,i]]])))
  }
  return(structure(all_but_one_res, class="coco_all_but_one"))
}