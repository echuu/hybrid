# function for rejection sampling

# input: k, i , m, c

# basically, we will sampling from the following prob kernel for the variable l
# 1_{l>0} l^{k - i} exp(-0.5 * (l - m)^2/c)



library(ars)


# density function for the kernel we need to sample from (mod a mult. const)
dkernel  = function(x, k, i, mu, c){
  
  return( x^{k - i} *exp(-0.5 * (x - mu)^2/c))
  
}



# log kernel density function for the kernel (mod additive const.)
ldkernel = function(x, k, i, mu, c){
  return((k - i)*log(x) - 0.5* (x - mu)^2/c)
  
}


# derivative of the log kernel density fn.
ldkernelprima = function(x, k, i, mu, c){
    return((k - i)/x - (x- mu)/c )
  
}


#We will pick the proposal density g to be gamma rv with shape = k - i + 1 and scale = 2*c

# Here we create the function  kernel/ g
ratiofn = function(x, k, i, mu, c){
  
  const = gamma(k - i +1)*(2*c)^(k - i +1)
  
  return(  const * exp(0.5*(x/c - (x - mu)^2/c)))
  
}

# 
# k = 7
# c= 3
# mu  = 1
# i = 3
# # 
# # 
# 
# # a bound for the ratio kernel / g
# M  = -optim(par =1 , fn = ratiofn, k  = k, i = i, mu= mu, c=c)$value + 100
# 
# 
# # accept prob = integral of kernel / M
# acceptprob= integrate (dkernel, lower = 0, upper = Inf, k= k, i =i, c=c, mu =mu)$value/M
# accept = 0
# sample  = numeric()
# count = 0
# 
# sampsize = 1000
# 
# while (accept < sampsize){
#   
#   count = count + 1
#   
#   u = runif(1)
#   y = rgamma(1, shape = k - i + 1, scale = c)
#   
#   if(u <= dkernel(y, k, i, mu, c)/(M*dgamma(y, shape = k - i + 1, scale = c))){
#     
#     sample = c(sample, y)
#     
#     accept = accept + 1
#     
#   }
#   
# }
# par(mfrow  = c(1, 2))
# 
# plot(seq(0, 10, by = .01), dkernel(seq(0, 10, by = .01), k, i, mu, c), 
#      type = "l",
#      main = "kernel",
#      ylab  = "y")
# 
# hist(sample, breaks = 100, freq = F, xlim = c(0,  10))
# 
# sampsize / count
# acceptprob

###################################################### 
#################ignore these for now###############


# # attempted to use the ars function for adaptive rejection samplingm but failed
#  arsSamp <- ars(20000, ldkernel, 
#                 ldkernelprima,x=c(0.1, 1000), m = length(c(0.1, 1000)), lb=TRUE,ub = FALSE,
#                 xlb= 0, k = k, i= i, mu= mu, c = c)
# 

##########################################################
######################################################

#the rejection sampling function

rkernel  = function(n, k, i, mu, c){
  
  ratiofn = function(x, k, i, mu, c){
    
    const = gamma(k - i +1)*(2*c)^(k - i +1)
    
    return(  const * exp(0.5*(x/c - (x - mu)^2/c)))
    
  }
  
  
  M  = ratiofn(x  = (1+ 2* mu)/2, k, i, mu, c) # analytic solution for the max of ratiofn
  acceptprob= integrate (dkernel, lower = 0, upper = Inf, k= k, i =i, c=c, mu =mu)$value/M
  
  cat("The calculated accept rate is: ", acceptprob, "\n")
  accept = 0
  sample  = numeric()
  count = 0
  
  
  
  while (accept < n){
    
    count = count + 1
    
    u = runif(1)
    y = rgamma(1, shape = k - i + 1, scale = c)
    
    if(u <= dkernel(y, k, i, mu, c)/(M*dgamma(y, shape = k - i + 1, scale = c))){
      
      sample = c(sample, y)
      
      accept = accept + 1
      
    }
    
  }

  return(sample)
}
# 
# data = rkernel (10000, k, i, mu, c)
# 
# 
# 
# par(mfrow = c(1, 2))
# hist(arsSamp, freq = FALSE, breaks = 100)
# 
# plot(seq(0.01, 11, 0.01), 
#      dkernel(seq(0.01, 11,  0.01), 
#              k, i, mu, c), type =  "l")
# 
