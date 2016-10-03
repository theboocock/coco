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

step_conditional = function(betas, ld_matrix,neffs,var_y, hwe_diag_outside,hwe_diag,return_entire_beta_set=F,exact=F,ses=NULL){
  # todo add the new math.
  inside =  ld_matrix
  n_betas = length(betas)
  if(n_betas == 1){
    return(cbind(betas,ses))
  }
  outside = sqrt(diag(hwe_diag_outside)) %*% inside %*% sqrt(diag(hwe_diag_outside))
  #print(outside)

  if(!exact){
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
  }
  #print( sqrt(diag(hwe_diag)) %*% inside %*% sqrt(diag(hwe_diag)))
  #  beta_inv = chol2inv(chol(outside))
 # outside = sqrt(diag(hwe_diag_outside)) %*% inside %*% sqrt(diag(hwe_diag_outside))
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
  if(exact){
    vars = (var_y * (median(neffs))- t(new_beta_joint) %*%  diag_one %*% ((betas[2:length(betas)]))- new_beta^2 * hwe_diag[1]) / (median(neffs) -2 )
  }else{
    vars = c(var_y)
  }
  
  
  ses = sqrt(diag(vars[1] *beta_two_inv))
  if(return_entire_beta_set){
    return(cbind(new_beta[,1],ses))  
  }else{
    return(c(new_beta[1,1],ses[1]))
  }
}