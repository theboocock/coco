#" Internal function to estimate the variance in Y using only the summary statistics.
#'
#' \code{estimate_vary} Estimate trait variance from summary statistics.
#' @title estimate_vary
#' @author James Boocock and Eli Stahl
#' @date 30 Aug 2016
#' @title estimate_vary
#' @param data_set data.frame. The summary statistics data.frame.
#'
#' Note: Below equation 8 in the cojo paper, is how we estimate the trait variance. \cr
#'   y'y = D_j S_j^2 (n-1) + D_j B_j^2 \cr
#'   Trait variance taken as the median of the equation above.
#'
#' @export
estimate_vary = function(data_set){
  vars=data_set$var *(data_set$n) * data_set$se ^2 * (data_set$n -1)  +  data_set$var * (data_set$n) * data_set$b^2
  return(median(vars/(data_set$n-1)))
}