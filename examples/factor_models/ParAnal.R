#### Function performing the parallell MCMC analysis

#input: K: vector of  model dimensions
#       y: the data matrix
#       a, b: positive scale parameters to be specified
#       And all other arguments needed for factGibbs function
#       beta0.list is a list of init beta0 for model of different dimension
#


#note that all the entries in the vector K must be at least 1


ParAnal.MCMC<- function(K, y,  beta0.list, 
                        Sigma0.list, nu, mu0,
                        No.Iter, C0, s, thin, burnin, seed ){
  
  library(tmvtnorm)
  library(mvtnorm)
  library(modeest)
  library(msm)
  library(parallel)
  library(pscl)
  library(MCMCpack)
  source("GibbsFact.R")
  source("BetaMatrixBack.R")
  source("LogPriorGivenK.R")
  source("LogQsubK.R")
  source("RJMCMC.R")
  source("MyTrunNorm.R")
  
  
  m= dim(y)[2]
  
  #### all the list objects use to store parallel anal info
  Gibbs.Store = vector(mode = "list", length = length(K))
  b.list= vector(mode = "list", length = length(K))
  B.list= vector(mode = "list", length = length(K))
  v.list= vector(mode = "list", length = length(K))
  output = list()
  ##############################3
  
  extract = function(L){return(L$M)} #L is object of type mlv
  
#  for(k in K){
#     
#     
#   Gibbs.Store.Temp= factGibbs(y, k, beta0.list[[k]], Sigma0, nu, mu0, C0, s, No.Iter)
#   Gibbs.Store[[k]] = Gibbs.Store.Temp  # store the posterior results for the k factor model
#   
#   
#   ###These are the indicators for the non zero entries of beta
#   
#   mat.ind = matrix(1, nr = m, nc = k)  # a m by k matrix with all 1
#   
#   if(k != 1)
#   for(row in 1:(k-1)) {mat.ind[row, (row+1):k ] =  rep(0, k-row)}
#   
#   
#   ind= which(c(mat.ind)==1)
#   
#   
#   
#   ######
#   
#   ## note we only return truncated entires here in the following, as a vector
#   b.list[[k]] = rowMeans(do.call(cbind, lapply(Gibbs.Store.Temp$beta, c)))[ind]
#   B.list[[k]] = var (t( do.call(cbind, lapply(Gibbs.Store.Temp$beta, c))[ind, ] ) )
#   ######
#   # grab the posterior MODE!
#   mode.list= apply(do.call(cbind, Gibbs.Store.Temp$Sigma), MARGIN = 1, FUN = mlv, 
#              method = "naive")
#   
#   v.list[[k]] = sapply(mode.list, FUN = extract)
#   }
#   
#    
#   
#   
#   output$b.info = b.list
#   output$B.info = B.list
#   output$v.info = v.list
#   output$sim.info = Gibbs.Store
#   
#   return(output)
#   
  
  
  #############################################################
  
  # The following part of the code  does the same thing as before, but uses the MCMCpack
  
  
  arrange.beta <- function(vec, k, m){
    
    
    carrier = rep(0, k*m)
carrier[ (k*k + 1):(k*m)] = vec[- (1:    (k*(k+1)/2))]
    
    ind = 0 
    for(i in 1 :k){
      carrier[(k*(i- 1) + 1):(k*(i- 1) + i)] =  vec[(ind + 1): (ind+ i)] 
      ind = ind + i      
    }  
    
    beta.wanted = c(t(matrix(carrier , nr = k)))
    
    
    return(beta.wanted)
  }
  
  
  
  MCMCobj = vector(mode = "list", length = length(K))
  
  for(k in K){
    
    # a list of constraints
    constr = list()
    nconstr = list()
    
    
    for(i in 1:k){
      
      for(j in i:k){
        
        nconstr  = c(nconstr, paste("V", i, sep = ""))
   
        if(i == j)
           constr = c(constr , list(list(i, "+")))
        else
           constr = c(constr, list(list(j, 0)))
       }
    }
 
    names(constr) = nconstr
  
  MCMCobj[[k]] = MCMCfactanal(x = y,
                 factors = k, 
                 lambda.constraints = constr,
                 lambda.start = beta0.list[[k]],
                 psi.start =  Sigma0.list[[k]] ,            
                 std.var = FALSE, 
                 store.scores=FALSE,
                 l0 = mu0,
                 L0 = C0, 
                 a0 = nu, 
                 b0 = nu*s^2,
                 mcmc = No.Iter,
                 thin = thin,
                 burnin = burnin,
                 seed = seed
                  )
  
    tot.no.var = dim(MCMCobj[[k]])[[2]]
    
    
    #################### all these stuffs take out the mode of the sigma's 
    mat.of.sigma = t(MCMCobj[[k]])[(tot.no.var -(m-1)) : tot.no.var, ]
    
    mode.list= apply(  mat.of.sigma, MARGIN = 1, FUN = mlv, 
                     method = "naive")
    
    v.list[[k]] = sapply(mode.list, FUN = extract)
    
    Gibbs.Store[[k]]= list()
    Gibbs.Store[[k]]$Sigma = as.list(data.frame(mat.of.sigma))
    #####################
    
    ###############all these stuffs take out the mean and cov of beta###
    
    mat.of.beta = t(MCMCobj[[k]])[ -((tot.no.var -(m-1)) : tot.no.var), ]
        
    mat.of.beta2 = apply(mat.of.beta, MARGIN = 2, FUN = arrange.beta, k = k , m= m  )
    
    
    mat.ind = matrix(1, nr = m, nc = k)  # a m by k matrix with all 1
    if(k != 1)
      for(row in 1:(k-1)) {mat.ind[row, (row+1):k ] =  rep(0, k-row)}
    ind= which(c(mat.ind)==1)
    
    b.list[[k]] = rowMeans(mat.of.beta2)[ind]
    B.list[[k]] = var (t( mat.of.beta2[ind, ] ) )


    Gibbs.Store[[k]]$beta = list()
    Gibbs.Store[[k]]$beta = 
      lapply(as.list(data.frame(mat.of.beta2)), FUN = matrix , nr = m)
    #######################
    
  
  
  }
  
  

  output$b.info = b.list
  output$B.info = B.list
  output$v.info = v.list
  output$sim.info = Gibbs.Store
  
  
  
  
  return (output)
  
  
    
}
