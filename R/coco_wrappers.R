#' All but one coco conditional analysis wrappper.
#' 
#' \code{all_but_one_coco} Performs a full all_but_one conditional analysis from the raw data.
#' @title all_but_one_coco
#' @date 30 Aug 2016
#' 
#' @author James Boocock and Eli Stahl
#' @date 2 Aug 2016
#' @title all_but_one_coco
#' @param data_set data.frame. Coco formatted input data.
#' @param ld_matrix matrix. LD matrix ordered by postiion
#' @param stepwise_p float. P-value threshold for forward-stepwise step.
#' @param var_y float. Phenotypic variance for Y, if not specified it will be estimated using \code{\link{estimate_vary}}.
#' @param hwe_variance boolean. Force the use of Hardy-weinberg equillibrium when calculating the genotypic variance.
#' @param exact boolean. Perform an exact analysis, extra data checks used in \code{\link{prep_data_set_common}}
#' @param joint_step boolean
#' \itemize{
#'  \item(T) = perform joint estimation when doing stepwise selection.
#'  \item(F) = perform conditional estimation when doing stepwise selection.
#'  }  
#' @param joint_all boolean
#' \itemize{
#'  \item(T) = perform joint estimation when performing all but one analysis
#'  \item(F) = perform conditional estimation when performing all but one analysis.
#'  }  
#' @export
all_but_one_coco = function(data_set ,ld_matrix,stepwise_p = 1e-3,var_y = 1.6421, joint_step=T,joint_all=T, reported_p = 1e-3){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix,ld_noise=ld_noise, var_y = var_y)
  stepwise_results = stepwise_run(res_preparation = res_preparation,p_value_threshold = stepwise_p, joint=joint_step)
  all_but_one_df = all_but_one(res_preparation=res_preparation,stepwise_results = stepwise_results, joint=joint_all ,p_value_threshold = reported_p)
  return (all_but_one_df)
}

#' Stepwise coco parameter estimation wrapper.
#'
#' \code{stepwise_coco} wraps \code{\link{joint_from_ids}} for convinience.
#' 
#' @author James Boocock and Eli Stahl
#' @date 2 Aug 2016
#' @title stepwise_conditional
#' @param data_set data.frame. Coco formatted input data.
#' @param ld_matrix matrix. LD matrix ordered by postiion
#' @param p_value_threshold float P-value threshold for forward-stepwise selection.
#' @param rsids character vector. SNPIDs in res preparation to perform CoCo analysis on.
#' @param idxs integer vector - postiion in input data to condition on.
#' @param var_y float. Phenotypic variance for Y, if not specified it will be estimated using \code{\link{estimate_vary}}.
#' @param hwe_variance boolean. Force the use of Hardy-weinberg equillibrium when calculating the genotypic variance.
#' @param joint boolean
#' \itemize{
#'  \item(T) = perform joint estimation when doing stepwise selection.
#'  \item(F) = perform conditional estimation when doing stepwise selection.
#'  } 
#'  @export

stepwise_coco = function(data_set ,ld_matrix, rsids=NULL, idxs=NULL,var_y = NULL, exact=F,hwe_variance=T,joint=T){
  res_preparation = prep_dataset_common(data_set = data_set,ld_matrix= ld_matrix, var_y = var_y,hwe_variance = hwe_variance,exact = exact)
  stepwise_results = joint_from_ids(rsids=rsids, idxs=idxs, res_preparation = res_preparation, joint=joint, exact=excat)
  return (stepwise_results)
}