#### Function performing the parallell MCMC analysis

#input: K: vector of  model dimensions
#       y: the data matrix
#       a, b: positive scale parameters to be specified
#       And all other arguments needed for factGibbs function
#       beta0.list is a list of init beta0 for model of different dimension
#


#note that all the entries in the vector K must be at least 1


ParAnal.MCMC.Mod<- function(K, y,  beta0.list, Sigma0.list,
                             nu, No.Iter, C0, 
                            s, thin, burnin, seed ){
  
  
  library(mvtnorm)
  library(modeest) # for estimaing mode
  library(parallel)
  library(msm)  # for univariate truncated normal
  library(pscl)  # for inverse gamma rv (used in factGibbsMod)
  source("GibbsFactMod.R")

  
  
  m= dim(y)[2]
  
  #### all the list objects use to store parallel anal info
  Gibbs.Store = vector(mode = "list", length = length(K)) # a list storing length(K) Gibbs sampler objects
  b.list= vector(mode = "list", length = length(K)) # a list of length(K), each component of the list being a posterior mean VECTOR of non zero entries of beta
  B.list= vector(mode = "list", length = length(K)) # a list of length(K), each component of the list being a posterior cov matrix of non zero entries of beta
  v.list= vector(mode = "list", length = length(K))  # a list of length(K), each component being a vector of length m holding the mode of the posterior sigma's
  
  # the list used to store all the other list
  output = list()
  ##############################3
  
  extract = function(L){return(L$M)} #L is object of type mlv
  
 for(k in K){
    
    
  Gibbs.Store.Temp= factGibbsMod(y, k, beta0.list[[k]], Sigma0.list[[k]], nu, C0, 
                                 s, No.Iter, burnin, thin, seed  = seed)
  Gibbs.Store[[k]] = Gibbs.Store.Temp  # store the posterior results for the k factor model
  
  
  ###These are the indicators for the non zero entries of beta
  
  mat.ind = matrix(1, nr = m, nc = k)  # a m by k matrix with all 1
  
  if(k != 1)
  for(row in 1:(k-1)) {mat.ind[row, (row+1):k ] =  rep(0, k-row)}
  
  # dismantle mat.ind and bind all the columns to form a m-by-k vector, 
  # and find the index for entries= 1
  
  ind= which(c(mat.ind)==1) 
  
  
  
  ######
  
  ## note we only return truncated entires here in the following, as a vector
  b.list[[k]] = rowMeans(  do.call(cbind, lapply(Gibbs.Store.Temp$beta, c))  )[ind]
  B.list[[k]] = var (t(    do.call(cbind, lapply(Gibbs.Store.Temp$beta, c))[ind, ]     ) )
  ######
  # grab the posterior MODE!
  mode.list= apply(  do.call(cbind, Gibbs.Store.Temp$Sigma), MARGIN = 1, FUN = mlv, 
             method = "naive"   )
  
  v.list[[k]] = sapply(mode.list, FUN = extract)
  }
  
   
  
  
  output$b.info = b.list
  output$B.info = B.list
  output$v.info = v.list
  output$sim.info = Gibbs.Store
  
 
  
  
  
  return(output)
  
  
  
  #############################################################
  
  # This part of the code  does the same thing as before, but uses the MCMCpack
  
#   
#   arrange.beta <- function(vec, k, m){
#     
#     
#     carrier = rep(0, k*m)
# carrier[ (k*k + 1):(k*m)] = vec[- (1:    (k*(k+1)/2))]
#     
#     ind = 0 
#     for(i in 1 :k){
#       carrier[(k*(i- 1) + 1):(k*(i- 1) + i)] =  vec[(ind + 1): (ind+ i)] 
#       ind = ind + i      
#     }  
#     
#     beta.wanted = c(t(matrix(carrier , nr = k)))
#     
#     
#     return(beta.wanted)
#   }
#   
#   
#   
#   MCMCobj = list()
#   
#   for(k in K){
#     
#     # a list of constraints
#     constr = list()
#     nconstr = list()
#     
#     
#     for(i in 1:k){
#       
#       for(j in i:k){
#         
#         nconstr  = c(nconstr, paste("V", i, sep = ""))
#    
#         if(i == j)
#            constr = c(constr , list(list(i, "+")))
#         else
#            constr = c(constr, list(list(j, 0)))
#        }
#     }
#  
#     names(constr) = nconstr
#   
#   MCMCobj[[k]] = MCMCfactanal(x = y,
#                  factors = k, 
#                  lambda.constraints = constr,
#                  std.var = FALSE, 
#                  store.scores=FALSE,
#                  l0 = mu0,
#                  L0 = C0, 
#                  a0 = nu, 
#                  b0 = nu*s^2,
#                  mcmc = No.Iter,
#                  thin = thin,
#                  burnin = burnin
#                   )
#   
#     tot.no.var = dim(MCMCobj[[k]])[[2]]
#     
#     
#     #################### all these stuffs take out the mode of the sigma's 
#     mat.of.sigma = t(MCMCobj[[k]])[(tot.no.var -(m-1)) : tot.no.var, ]
#     
#     mode.list= apply(  mat.of.sigma, MARGIN = 1, FUN = mlv, 
#                      method = "naive")
#     
#     v.list[[k]] = sapply(mode.list, FUN = extract)
#     
#     Gibbs.Store[[k]]= list()
#     Gibbs.Store[[k]]$Sigma = as.list(data.frame(mat.of.sigma))
#     #####################
#     
#     ###############all these stuffs take out the mean and cov of beta###
#     
#     mat.of.beta = t(MCMCobj[[k]])[ -((tot.no.var -(m-1)) : tot.no.var), ]
#         
#     mat.of.beta2 = apply(mat.of.beta, MARGIN = 2, FUN = arrange.beta, k = k , m= m  )
#     
#     
#     mat.ind = matrix(1, nr = m, nc = k)  # a m by k matrix with all 1
#     if(k != 1)
#       for(row in 1:(k-1)) {mat.ind[row, (row+1):k ] =  rep(0, k-row)}
#     ind= which(c(mat.ind)==1)
#     
#     b.list[[k]] = rowMeans(mat.of.beta2)[ind]
#     B.list[[k]] = var (t( mat.of.beta2[ind, ] ) )
# 
# 
#     Gibbs.Store[[k]]$beta = list()
#     Gibbs.Store[[k]]$beta = 
#       lapply(as.list(data.frame(mat.of.beta2)), FUN = matrix , nr = m)
#     #######################
#     
#   
#   
#   }
#   
#   
# 
#   output$b.info = b.list
#   output$B.info = B.list
#   output$v.info = v.list
#   output$sim.info = Gibbs.Store
#   
#   
#   
#   
#   return (output)
  
  
    
}
