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

r_inverse= function(gamma,phi,ref_af){
  B =  (0.5*(1-2*phi) * (1 - 2*ref_af) * gamma -1)
  A = - (0.084 + 0.9 * phi * (1-2 * phi) * ref_af * (1-ref_af)) * gamma
  C = gamma
  x1 = (-B+ sqrt( B^2 - 4 * A * C  ))/(2*A)
  x2 = (-B - sqrt( B^2 - 4 * A * C  ))/(2*A)
  return(cbind(x1,x2))
}
r =  function(gamma,phi,theta){
  r = gamma/(((phi*(1-phi)) + 0.5*(1-2*phi) * (1 - 2*theta) * gamma) - (0.084 + 0.9 * phi * (1-2 * phi) * theta * (1-theta)) /(phi*(1-phi))* gamma^2)
  r
}
