# a function to generate truncated multivariate normal (negative truncations)


# input: mean, covariance, and position vector

# indn vector is an indicator vector pointing to where we want to truncate in our random vector

# use library msm, use the rtnorm function there

my.rtmvnorm= function(mean, Sigma, ind, m, k){
  if(m <= k) {
    print("something goes wrong")
  return()}
  
  
  d = length(mean)
  
  no.trun = sum(ind)  # non of truncated coordinates
  
  out.vec = numeric(d)
  

  pos = which(ind ==1)
  
  for(i in 1: no.trun){
    
    
    if(i ==1){
      out.vec[1] = rtnorm (n =1,mean = mean[1], sd = sqrt(Sigma[1, 1]) ,lower = 0) 
      if(no.trun != 1){
      out.vec[2:(pos[2]-1)] = rmvnorm(n =1, 
                                  mean = mean[2:(pos[2]-1)] + (1/Sigma[1, 1])*(out.vec[1]-mean[1])*Sigma[2:(pos[2]-1), 1], 
                                  sigma = Sigma[2:(pos[2]-1), 2:(pos[2]-1)] - (1/Sigma[1, 1])*tcrossprod(Sigma[2:(pos[2]-1), 1]))
      }else{
      out.vec[2:d] = rmvnorm(n=1, 
                             mean = mean[2:d] + (1/Sigma[1, 1])*(out.vec[1]-mean[1])*Sigma[2:d, 1], 
                             sigma = Sigma[2:d, 2:d] - (1/Sigma[1, 1])*tcrossprod(Sigma[2:d, 1]))  
      }
      
    }else if (i>1 && i < no.trun ){
      
      
  req.mean = mean[pos[i]:(pos[i+1]-1)]+
            Sigma[ pos[i]:(pos[i+1]-1), 1:(pos[i]-1)]%*%
            solve(Sigma[1:(pos[i]-1), 1:(pos[i]-1)])%*%(out.vec[1:(pos[i]-1)]- mean[1:(pos[i]-1)])
    
  req.cov = Sigma[ pos[i]:(pos[i+1]-1), pos[i]:(pos[i+1]-1)] -
            Sigma[ pos[i]:(pos[i+1]-1) ,1:(pos[i]-1)]%*%
            solve(Sigma[1:(pos[i]-1),  1:(pos[i]-1)])%*%
            Sigma[1:(pos[i]-1),  pos[i]:(pos[i+1]-1)]
    
  sub.d = length(req.mean)  
    
  out.vec[pos[i]] = rtnorm (n =1,mean =req.mean[1] , sd = sqrt(req.cov[1, 1]) ,lower = 0) # univariate truncated normal
 
  out.vec[(pos[i]+1):(pos[i+1]-1)]  = 
    rmvnorm(n =1, 
            mean = req.mean[2:sub.d] + req.cov[2:sub.d, 1]*(1/req.cov[1, 1])*(out.vec[pos[i]] - req.mean[1]), 
            sigma= req.cov[2:sub.d, 2:sub.d]- (1/req.cov[1, 1])*tcrossprod(req.cov[2:sub.d, 1], req.cov[1, 2:sub.d])
            ) # mult.var conditonal normal
    
    
  }else{  # last column
    
    
    req.mean = mean[pos[i]:d]+
      
      Sigma[ pos[i]:d ,1:(pos[i]-1)]%*%
      solve(Sigma[1:(pos[i]-1), 1:(pos[i]-1)]) %*%
      (out.vec[1:(pos[i]-1)]- mean[1:(pos[i]-1)])
    
    req.cov = Sigma[ pos[i]:d, pos[i]:d] -
      Sigma[ pos[i]:d ,1:(pos[i]-1)]%*%
      solve(Sigma[1:(pos[i]-1),  1:(pos[i]-1)])%*%
      Sigma[1:(pos[i]-1),  pos[i]:d]
    
    sub.d = length(req.mean)
    
    
    out.vec[pos[i]] = rtnorm (n =1,mean =req.mean[1] , sd = sqrt(req.cov[1, 1]) ,lower = 0) # univariate truncated normal
    
    out.vec[(pos[i]+1):d]  = 
      rmvnorm(n =1, 
              mean = req.mean[2:sub.d] + req.cov[2:sub.d, 1]*(1/req.cov[1, 1])*(out.vec[pos[i]] - req.mean[1]), 
              sigma= req.cov[2:sub.d, 2:sub.d]- (1/req.cov[1, 1])*tcrossprod(req.cov[2:sub.d, 1], req.cov[1, 2:sub.d])
            ) # mult.var conditonal normal
    
    
  }  
    
    
  }
  
  
  return(out.vec)
}

