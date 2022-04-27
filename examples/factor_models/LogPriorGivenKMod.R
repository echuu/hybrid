



# function calculating the log prior of theta(beta, Sigma )given k
####################################################################



ldthetaMod = function(beta, Sigma , C0, nu, s){
  # beta is a m by k, lower block triangular matrix
  # Sigma is a m vector
  
  
  if(is.na(beta)[1]){
    
    loglike = 0 
    m = length(Sigma)  
    for(j in 1:m) loglike = loglike + log(densigamma(Sigma[i], alpha = nu/2, beta = nu*s^2/2))
    return(loglike)
  }
  
  
  m = dim(beta)[1]
  k = dim(beta)[2]
  loglike = 0 
  
  # deal with beta first
  for (i in 1:m){
    
    for(j in 1: min(k, i)){  # j is column index
      if(i == j){
        loglike = loglike +  y3(1/C0*beta[i, j]^2, df = k - i + 1, log = TRUE) + 
          log(2) + log(beta[i, j])- log(C0)
        # loglikelihood contributed by diagonal entries
      }else{
        
        # loglikelihood contributed by other entries
        loglike = loglike+ dnorm(beta[i, j], mean = 0, 
                                 sd = sqrt(C0), log = TRUE)
        
      }
    }
  }
  
  
  #deal with Sigma next
  for (i in 1:m){
    loglike = loglike + log(densigamma(Sigma[i], alpha = nu/2, beta = nu*s^2/2))
  }
  return(loglike)
  
}
##################################################
