##' Internal function joint_step.
##'
##' \code{step_joint}. Performs a joint modelling step to estimate new beta values.
##'
##'
##' @author James Boocock
##' @title step_conditional
##' @param betas vector of betas that will be adjusted in the matrix operations.
##' @param ld_matrix LD matrix for the subset of the region.
##' @param neff effective sample size for each SNP in the dataset
##' @param var_y variance of of the phenotype
##' @param hwe_diag_outside Diagonal matrix containing the genotypic variance
##' @param hwe_diag Diagonal matrix containing the genotypic variance *n
##' @return Joint betas and standard errors.
##' 
##' @export
step_joint = function(betas, ld_matrix,neffs,var_y, hwe_diag_outside,hwe_diag,return_entire_beta_set=F,exact=F){
  inside = ld_matrix
  n_betas = length(betas)
  outside = sqrt(diag(hwe_diag_outside)) %*% inside %*% sqrt(diag(hwe_diag_outside))
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
  #print(outside)

  beta_inv = chol2inv(chol(outside))
  #print(beta_two_inv)
  #print(betas)
  new_betas = beta_inv %*% diag(hwe_diag) %*% betas
  #print(beta_inv)
  # We only really care about the results from the first SNP
  #because that's our SNP we are adding to the model
  if(exact){
    vars = (var_y * (median(neffs))- t(new_betas) %*%  diag(hwe_diag) %*% ((betas))) / (median(neffs) - n_betas )
  }else{
    vars = c(var_y)
  }
    ses = sqrt(diag(vars[1] * beta_inv))
  
  if(return_entire_beta_set){
    return(cbind(new_betas[,1],ses))  
  }else{
    return(c(new_betas[1,1],ses[1]))
  }
}
