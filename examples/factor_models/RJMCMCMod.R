### This is a RJMCMC function



#input:




RJMCMC.Mod = function(y, P.anal, K, a, b, k.start,No.Iter.Anal, RJ.iter, K.prob, tran.prob,
                  nu, s,  C0){
  
  
  m = dim(y)[2]
  T = dim(y)[1]
  
  # to store the RJMCMC parameter values
  theta.beta = list()
  theta.Sigma = list()
  dim.rec = numeric()
  #################
  

  k.cur = k.start  # set the current value for k 
  
  dim.rec[1] = k.cur
  


  # make a draw with k.start
  Gibbs.Out.start = factGibbsMod(y= y, k= k.cur, 
                              beta0 = P.anal$sim.info[[k.cur]]$beta[[No.Iter.Anal]], 
                              Sigma0 = P.anal$sim.info[[k.cur]]$Sigma[[No.Iter.Anal]], 
                              nu = nu, 
                              C0 = C0, 
                              s = s, 
                              No.Iter = 1,
                              burnin = 0, 
                              thin = 1   )
  
  # initialize the parameter store list
  
  theta.beta[[1]]=Gibbs.Out.start$beta[[1]] 
  beta.cur = Gibbs.Out.start$beta[[1]] #set a current value for beta
  
  
  theta.Sigma[[1]] = Gibbs.Out.start$Sigma[[1]]
  Sigma.cur = Gibbs.Out.start$Sigma[[1]] # set current value for Sigma
  
  
  
  
  
  
  ####################################
  
  count = 1
  
  while(count < RJ.iter){  
    
    draw.prob = tran.prob[which(rownames(tran.prob) == k.cur), ]
    
    k.cand = K[sample(length(K), size =1 , prob = draw.prob)] 
    # candidate for the next k
    
    # candidate for new beta
    
    if(k.cand == 0){
      
    beta.cand = NA
      
    }else{
    
    ##### set lower bound for truncated mult. normal when drawing beta candidate
    lower.vec  = rep(- Inf, length(P.anal$b.info[[k.cand]])) 
    
    ind.vec = numeric()
    for(t in 1:k.cand){
      ind.vec = c(ind.vec , c(1, rep(0, m-t)))
    }
    lower.vec[which(ind.vec==1)] = 0
    ########
    beta.cand = beta.back(m = m, 
                          k = k.cand, 
                          my.rtmvnorm(mean= P.anal$b.info[[k.cand]], 
                                      Sigma= b* P.anal$B.info[[k.cand]],
                                      ind = ind.vec,
                                      m = m,
                                      k= k.cand)
    )
    
    
    
    #beta.cand = beta.back(m = m, 
    #                      k = k.cand, 
    #                      rtmvnorm(n = 1, 
    #                               mean =P.anal$b.info[[k.cand]], 
    #                               sigma = b* P.anal$B.info[[k.cand]] , 
    #                               lower =lower.vec,
    #                               algorithm = "gibbs", 
    #                               burn.in.samples = 10000)
    #)
    }
    
    # candidate for new Sigma
    
    Sigma.cand = numeric()
    
    if(k.cand == 0){
      
    for(j in 1:m)  Sigma.cand[j] = rigamma(n =1, alpha = (nu + T)/2, 
                                             beta = (sum(y[,j]^2)+ nu*s^2)/2)
      
    }else{
    for(j in 1:m) Sigma.cand[j] = rigamma(n = 1, alpha = a, 
                                             beta = a*P.anal$v.info[[k.cand]][j])
    }
    # calculate log accept ratio
    
    
    
    #log of numerator
    if(k.cand ==0){ 
      log.num = sum(dmvnorm(x = y, sigma =diag(Sigma.cand), log = TRUE))
      }else{
      log.num = sum(dmvnorm(x = y, sigma = beta.cand%*%t(beta.cand) + diag(Sigma.cand), log = TRUE))
      }
    
    
    log.num= 
      log.num + 
      ldthetaMod(beta.cand, Sigma.cand ,  C0, nu, s) + 
      
      log(K.prob[which(K == k.cand)])+
      
      logQsubK(beta = beta.cur, Sigma = Sigma.cur, a, b, Bsubk =  P.anal$B.info[[k.cur]], 
              bsubk = P.anal$b.info[[k.cur]], v = P.anal$v.info[[k.cur]] , y,  nu , s) + 
      
      log(tran.prob[as.character(k.cand),as.character(k.cur) ]) 
    
    #log of denominator
    
    if(k.cur == 0) {log.den = sum(dmvnorm(x = y, sigma =  diag(Sigma.cur), log = TRUE) )
    }else {
      log.den = sum(dmvnorm(x = y, sigma = beta.cur%*%t(beta.cur) + diag(Sigma.cur), log = TRUE) )
    }
    
    log.den = 
      log.den+
      ldthetaMod (beta.cur, Sigma.cur,  C0, nu, s)+
      log(K.prob[which(K == k.cur)]) +
      logQsubK(beta = beta.cand, Sigma = Sigma.cand, a, b, Bsubk = P.anal$B.info[[k.cand]], 
              bsubk = P.anal$b.info[[k.cand]], v =  P.anal$v.info[[k.cand]], y,  nu , s) +   
      log( tran.prob[as.character(k.cur), as.character(k.cand)] )
    
    # calculate the log acceptance ratio
    
    log.ratio = log.num - log.den
    
    
    if (log(runif(1)) <= log.ratio) { #we accept
      
      k.cur= k.cand
      
      
      count = count +1

      dim.rec[[count]]=k.cur
      ######## Run an additional within model step
      
      Gibb.obj = factGibbsMod(y, k.cur, beta0=beta.cand, 
                           Sigma0 = Sigma.cand, nu,  C0, s, No.Iter= 1,
                              burnin = 0, 
                              thin = 1  )  
      
      
      
      beta.cur = Gibb.obj$beta[[1]]
      Sigma.cur =  Gibb.obj$Sigma[[1]]
      
      
      theta.beta[[count]] = beta.cur
      theta.Sigma[[count] ]= Sigma.cur
      
      
    }else{
      
      count= count +1
      dim.rec[[count]]=k.cur
      
      
      Gibb.obj = factGibbsMod(y, k.cur, beta0  =beta.cur, 
                           Sigma0 = Sigma.cur, nu, C0, s, No.Iter= 1,
                              burnin = 0, 
                              thin = 1  )  
      
      
      beta.cur = Gibb.obj$beta[[1]]
      Sigma.cur =  Gibb.obj$Sigma[[1]]
      theta.beta[[count]] = beta.cur
      theta.Sigma[[count]]= Sigma.cur
      
      
    }
    
    
    print( dim.rec[[count]])
    
  } # end of while loop 
  
  
  output = list()
  output$beta = theta.beta
  output$Sigma = theta.Sigma
  output$dim.rec = dim.rec
  return(output)
  
  
}