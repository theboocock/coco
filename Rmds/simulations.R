require(bindata)

simulate_phenotypes = function(x){

}

##' Simulate genotypes using a random LD matrix
##'
##'
##'
##' @author James Boocock
##' @date 4 Aug 2016
##' @title stepwise_conditional
##' @param data_set
##'

x= function (commonprob)
{
  retval <- TRUE
  message <- character(0)
  nm <- 0
  if ((any(commonprob < 0)) || (any(commonprob > 1))) {
    retval <- FALSE
    message[nm <- nm + 1] <- "Not all probabilities are between 0 and 1."
  }
  n <- dim(commonprob)[1]
  if (n != dim(commonprob)[2]) {
    retval <- FALSE
    message[nm <- nm + 1] <- "Matrix of common probabilities is not quadratic."
  }
  for (i in 1:(n - 1)) {
    for (j in 2:n) {
      ul <- min(commonprob[i, i], commonprob[j, j])
      ll <- max(commonprob[i, i] + commonprob[j, j] - 1,
                0)
      print(ul)
      print(ll)
      print(commonprob)
      if ((commonprob[i, j] > ul) || (commonprob[i, j] <
                                      ll)) {
        retval <- FALSE
        message((commonprob[i, j] > ul))
        message(((commonprob[i, j] < ll)))
        message[nm <- nm + 1] <- paste("Error in Element (",
                                       i, ",", j, "): Admissible values are in [",
                                       ll, ",", ul, "].")
      }
    }
  }
  if (n > 2)
    for (i in 1:(n - 2)) for (j in (i + 1):(n - 1)) for (k in (j +
                                                               1):n) {
      l <- commonprob[i, i] + commonprob[j, j] + commonprob[k,
                                                            k] - 1
      if (commonprob[i, j] + commonprob[i, k] + commonprob[j,
                                                           k] < l) {
        retval <- FALSE
        message[nm <- nm + 1] <- paste("The sum of the common probabilities of",
                                       i, ",", j, ",", k, "must be at least", l, ".")
      }
    }
  attr(retval, "message") <- message
  retval
}



##' Get allowable correlations
##'
##'
##'
##' @author James Boocock
##' @date 4 Aug 2016
##' @title get_allowable_correlations
##' @param probs - marginal probabilities.
##'

get_allowable_correlations = function(probs,corr){
  # For all pairs of SNPs.
  # Convert probs to af.
  # This should also make the matrix symmetric.
  afs = ifelse(probs > 0.5, 1-probs,afs)
  for(i in 1:(length(probs)-1)){
    for(j in (i+1):length(probs)){
        min_af = min(afs[c(i,j)])
        maximum_corr = (min_af - afs[i] * afs[j])/ sqrt(afs[i] * (1-afs[i])* afs[i] * (1-afs[j]))
        corr[j,i] = maximum_corr * corr[i, j]
        corr[i,j] = maximum_corr * corr[i, j]
    }
  }
  return(corr)
}

rgeno = function(n,probs,corr,min_maf=0.05,size=10000){
  if(missing(probs)){
    probs = runif(n,0.05,0.95)
  }
  if(missing(corr)){
    corr = matrix(runif(n*n,0,1), nrow=n,ncol=n)
    corr[c(1:ncol(corr)),c(1:nrow(corr))] = 1
    corr= get_allowable_correlations(probs,corr)
  }
  h1 = rmvbin(n = size,margprob = probs,bincorr = corr)
  h2 = rmvbin(n = size,margprob = probs,bincorr = corr)
  cbind(h1 + h2)
}
