##' Stepwise conditional wrapper
##'
##' Generate conditional signals from a conditional dataset.
##'
##'
##' @author James Boocock and Eli Stahl
##' @date 2 Aug 2016
##' @title estimate_betas
##' @param data_set data.frame. Contains the relevant
##' @param ld_matrix
##' @param p_value_threshold.
##' @param ld_noise
##'
##' citation. Pirinen, Matti, Peter Donnelly, and Chris CA Spencer. "Efficient computation with a linear mixed model on large-scale data sets
##' with applications to genetic studies." The Annals of Applied Statistics 7.1.

logistic_to_linear = function(b, se, phi, ref_af, mu){
  mu=phi # ??
  B =  (0.5*(1-2*phi) * (1 - 2*ref_af) * b -1)
  A = - (0.084 + 0.9 * phi * (1-2 * phi) * ref_af * (1-ref_af)) * b
  C = b
  x1 = (-B+ sqrt( B^2 - 4 * A * C  ))/(2*A)
  x2 = (-B - sqrt( B^2 - 4 * A * C  ))/(2*A)
  ## Always x2 ??
  #beta = x2 # inverse of GWAS_approximation
  alpha = log(mu/(1-mu))
  beta =  x2 * exp(alpha)/((1+exp(alpha))^2) # FOA
  ## Always x2 ??
  ## Always x2 ??
  alpha = log(mu/(1-mu))
  se = se^2 * ((0.5*(1-2*phi) * (1 - 2*ref_af) * x2) -((0.084 + 0.9 * phi * (1-2 * phi) * ref_af * (1-ref_af)) * x2^2) +  1) ^2 
  se =  sqrt(se* (exp(alpha)/((1+exp(alpha))^2))^2) # FOA
  return(cbind(beta,se))
}



linear_to_logistic =  function(b, se, gamma,phi,theta){
  # mu == phi ??
  mu = phi
  gamma_flat = b/(phi*(1-phi))
  
  gamma = b/(((phi*(1-phi)) + 0.5*(1-2*phi) * (1 - 2*theta) * b) - (0.084 + 0.9 * phi * (1-2 * phi) * theta * (1-theta)) /(phi*(1-phi))* b^2)
  se_fist = se^2* (1/(mu*(1-mu)))^2
  se_second = sqrt(se_fist * ((0.5*(1-2*phi) * (1 - 2*theta) * gamma_flat) -((0.084 + 0.9 * phi * (1-2 * phi) * theta * (1-theta)) *gamma_flat * 2 * gamma)  + 1 ) ^2 )
  se_fist = sqrt(se_fist)
  return(cbind(gamma,se_fist,se_second))
}


