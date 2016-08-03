##' Stepwise estimation.
##'
##' Generate conditional signals from a conditional dataset.
##'
##'
##' @author James Boocock
##' @date 2 Aug 2016
##' @title stepwise_conditional
##' @param data_set
##' @param p_value_threshold.
##'
##'

stepwise_conditional_run = function(data_set,ld_matrix ,p_value_threshold=0.00001,colinear_threshold=0.9,var_y = 1.6421){
  #extract variance from the dataset
  #hwe_diag = (2*freq_af$af * ( 1- data_set$af) * data_set$n)
  # Extract the effective sample size for a SNP.
  # Generate hwe_diagonal
  hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$n )
#  hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$n * data_set$info)
  data_set$neff = (var_y * data_set$n) / (hwe_diag *data_set$se^2  - data_set$b) / (data_set$se^2 +1)
  # Remove HWE diagonal
  hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$neff)
#  hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$neff * data_set$info)  ## ??
  # Get hwe D matrix

  # Get hwe D matri x without sample size, needed to generate B matrix.
  hwe_diag_outside=  (2*data_set$af * ( 1- data_set$af))
  message("Running a stepwise conditional analysis")
  if (!("Z" %in% names(data_set))){
    data_set$Z = data_set$b / data_set$se
  }
#  print(names(data_set))
  # ids of the conditional hits.
  # First find the maximum SNP.
#  print(max(abs(data_set$Z)))
  idx_top_tmp = c(which(max(abs(data_set$Z)) == abs(data_set$Z)))
  conditional_df = data_set[idx_top_tmp,]
#  print(idx_top_tmp)
  message(paste("Top SNP for region = ", conditional_df$rsid, ", with Z-score ", conditional_df$Z))
  current_best_p = 2 * pnorm(abs(conditional_df$Z), lower.tail = F)
  idx_cond = c()
  message(paste("Conditioning ..."))
  out_all_buts = data.frame(rsid=conditional_df$rsid, beta_old=conditional_df$b, beta_new=NA
                            ,se_old=conditional_df$se, se_new =NA,Znew=NA,p=current_best_p)
  while(current_best_p < p_value_threshold){
    idx_cond = c(idx_cond, idx_top_tmp)
    res_step = data.frame(rsid=data_set$rsid,beta_old=data_set$b, beta_new=rep(NA,nrow(ld_matrix)), se_old=data_set$se, se_new=rep(NA,nrow(ld_matrix)))
    for(i in 1:nrow(ld_matrix)){
        if(i %in% idx_cond){
          message("Skipping SNP as already being conditioned on")
          next
        }
        of_interest = c(i, idx_cond)
        betas = data_set$b[of_interest]
        tmp_ld_matrix = ld_matrix[of_interest,of_interest]
        neffs = data_set$neff[of_interest]
        if(any(abs(tmp_ld_matrix[1,2:ncol(tmp_ld_matrix)]) > colinear_threshold)){
          message("Skipping SNP as co-linear with top_snps.")
          next
        }

        hwe_diag_outside_tmp = hwe_diag_outside[of_interest]
        hwe_diag_tmp = hwe_diag[of_interest]
        cond_results = step_conditional(betas, tmp_ld_matrix, neffs,var_y, hwe_diag_outside_tmp,hwe_diag_tmp)
        res_step$beta_new[i]= cond_results[1]
        res_step$se_new[i]= cond_results[2]
    }
    res_step$Znew = res_step$beta_new/res_step$se_new
    idx_top_tmp = which(max(abs(res_step$Znew),na.rm=T) == abs(res_step$Znew))
    best_cond_row = res_step[idx_top_tmp,]
    current_best_p = 2 * pnorm(abs(best_cond_row$Znew), lower.tail = F)
    message(paste("Best SNP is = ", best_cond_row$rsid, "with P-value ", current_best_p, " Index = ", idx_top_tmp))
    message(paste("Beta = ", best_cond_row$beta_new, " Z = ",best_cond_row$Znew))
    if(current_best_p < p_value_threshold){
      best_cond_row$p = current_best_p
      out_all_buts = rbind(out_all_buts,best_cond_row)
    }

  }
  return(out_all_buts)
}

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

step_conditional = function(betas, ld_matrix,neffs,var_y, hwe_diag_outside,hwe_diag){
  inside =  ld_matrix
  n_betas = length(betas)
  outside = sqrt(diag(hwe_diag_outside)) %*% inside %*% sqrt(diag(hwe_diag_outside))

  # outside = sqrt(diag(hwe_diag[idx])) %*% inside %*% sqrt(diag(hwe_diag[idx]))
  #print(outside)
  for(j in 1:ncol(outside)){
    for(k in (j):ncol(outside)){
      if(k > ncol(outside)){
        break
      }
      outside[j,k] = min(neffs[c(j,k)]) *  outside[j,k]
      if(k!=j){
        outside[k,j] = min(neffs[c(j,k)]) * outside[k,j]
      }
    }
  }
  beta_inv = chol2inv(chol(outside))
  new_betas = beta_inv %*% diag(hwe_diag) %*% betas
  # We only really care about the results from the first SNP
  #because that's our SNP we are adding to the model
  neff_var_y = neffs[1] * var_y
  vars = (neff_var_y - t(new_betas) %*%  diag(hwe_diag) %*% ((betas))) / (neffs[1] - length(n_betas))
  ses = sqrt(diag(vars[1] * beta_inv))
  return(c(new_betas[1,1],ses[1]))
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
all_but_one = function(data_set,ld_matrix, stepwise_results,p_value_threshold=0.00001,colinear_threshold=0.9,var_y = 1.6421){
  #extract variance from the dataset
  #hwe_diag = (2*freq_af$af * ( 1- data_set$af) * data_set$n)
  # Extract the effective sample size for a SNP.
  # Generate hwe_diagonal
  hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$n)
  # print(hwe_diag)
  data_set$neff = (var_y * data_set$n) / (hwe_diag *data_set$se^2  - data_set$b) / (data_set$se^2 +1)
  # Remove HWE diagonal
  hwe_diag =  (2*data_set$af * ( 1- data_set$af) * data_set$neff)
  # Get hwe D matrix
  if(missing(stepwise_results)){
    stop("Require the results of a stepwise analysis to perform an all but one beta generation")
  }
  # Get hwe D matri x without sample size, needed to generate B matrix.
  hwe_diag_outside=  (2*data_set$af * ( 1- data_set$af))
  message("Running a stepwise conditional analysis")
  if (!("Z" %in% names(data_set))){
    data_set$Z = data_set$b / data_set$se
  }
  # List containing the results for all the all_but_one analyses.
  all_but_one_res = list()

  ##TODO make the all but one analysis happen. Basically it is just the joint model leaving out all the independent effects.

}
