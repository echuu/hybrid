# function calculating log density of data given k and theta (beta + Sigma)
###################################################################



logQsubK = function(beta, Sigma, a, b, v, Bsubk, bsubk, y,  nu , s){
  # beta is a m by k matrix
  #Sigma is a m vector
  # bsubk is a vector of dimension same as the number of !=0 entries in beta 
  # Bsubk is a matrix of length(bsubk) by length(bsubk)
  
  T= dim(y)[1]
  if(is.na(beta)[1]) {
    loglike = 0 
    m = length(Sigma)
    for(i in 1:m) loglike = loglike + log(densigamma(Sigma[i], alpha = (nu + T)/2,
                                                     beta =  (sum(y[,i]^2)+ nu*s^2)/2))
    return(loglike)
  }
  
  
  m= dim(beta)[1]
  k = dim(beta)[2]
  
  
  ###indicator for making beta into a vector
  mat.ind = matrix(1, nr = m, nc = k)
  if(k != 1) for(row in 1:(k-1)) {mat.ind[row, (row+1):k ] =  rep(0, k-row)}
  ind= which(c(mat.ind)==1)
  ############################
  
  
  ###indicator for which entries in the vector form of beta is positive
  lower.vec = rep(-Inf, length(ind))
  
  ind.vec = numeric()
  for(t in 1:k){
    ind.vec = c(ind.vec , c(1, rep(0, m-t)))
  }
  lower.vec[which(ind.vec==1)] = 0
  #####################################
  
  loglike = 0
  for(i in 1:m){ loglike = loglike + log(densigamma(Sigma[i], alpha = a, beta = a*v[i]))}
  loglike = loglike + dtmvnorm(x = c(beta)[ind], mean = bsubk, sigma = b*Bsubk, 
                               lower = lower.vec , log = TRUE)
  return(loglike)
}
#############################################################
