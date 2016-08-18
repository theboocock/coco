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
##'

logistic_to_linear = function(b, z, phi, ref_af, mu){
  mu=phi # ??
  B =  (0.5*(1-2*phi) * (1 - 2*ref_af) * b -1)
  A = - (0.084 + 0.9 * phi * (1-2 * phi) * ref_af * (1-ref_af)) * b
  C = b
  x1 = (-B+ sqrt( B^2 - 4 * A * C  ))/(2*A)
  x2 = (-B - sqrt( B^2 - 4 * A * C  ))/(2*A)
  ## Always x2 ??
  beta = x2 # inverse of GWAS_approximation
  print(x2[1])
  alpha = log(mu/(1-mu))
  beta =  beta * exp(alpha)/((1+exp(alpha))^2) # FOA
  ## Always x2 ??
  se = beta / z
  return(cbind(beta,se))
}

linear_to_logistic =  function(b, se, gamma,phi,theta){
  # mu == phi ??
  gamma = b/(((phi*(1-phi)) + 0.5*(1-2*phi) * (1 - 2*theta) * b) - (0.084 + 0.9 * phi * (1-2 * phi) * theta * (1-theta)) /(phi*(1-phi))* b^2)
  se = se/(((phi*(1-phi)) + 0.5*(1-2*phi) * (1 - 2*theta) * se) - (0.084 + 0.9 * phi * (1-2 * phi) * theta * (1-theta)) /(phi*(1-phi))* se^2)
  return(c(gamma,se))
}


