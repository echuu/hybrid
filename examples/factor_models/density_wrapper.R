




formatBeta = function(x) {
  x[lower.tri(x, diag = TRUE)]
}
psi = function(x, params) {
  -logpost(x, params)
}
grad = function(x, params) {
  -grad_logpost(x, params)
}
hess = function(x, params) {
  -hess_logpost(x, params)
}

grad = function(x, params) {
  pracma::grad(psi, x0 = x, params = params)
}
hess = function(x, params) {
  pracma::hessian(psi, x0 = x, params = params)
}


preprocess = function(samps, d, params) {
  psi_u = apply(samps, 1, psi, params = params)
  # (1.2) name columns so that values can be extracted by partition.R
  u_df_names = character(d + 1)
  for (i in 1:d) {
    u_df_names[i] = paste("u", i, sep = '')
  }
  u_df_names[d + 1] = "psi_u"

  # populate u_df
  u_df = cbind(samps, psi_u) # J x (D + 1)
  names(u_df) = u_df_names
  return(u_df)
}




############## code below is for testing purposes -- leave commented out #######

# log_density = function(x, params) {
#   -psi(x, params)
# }



# ldthetaMod = function(beta, Sigma , C0, nu, s){
#   # beta is a m by k, lower block triangular matrix
#   # Sigma is a m vector
#
#
#   if(is.na(beta)[1]){
#
#     loglike = 0
#     m = length(Sigma)
#     for(j in 1:m) loglike = loglike + log(densigamma(Sigma[i], alpha = nu/2, beta = nu*s^2/2))
#     return(loglike)
#   }
#
#   m = dim(beta)[1]
#   k = dim(beta)[2]
#   loglike = 0
#
#   # deal with beta first
#   for (i in 1:m){
#
#     for(j in 1: min(k, i)){  # j is column index
#       if(i == j){
#         loglike = loglike +  dchisq(1/C0*beta[i, j]^2, df = k - i + 1, log = TRUE) +
#           log(2) + log(beta[i, j])- log(C0)
#         # loglikelihood contributed by diagonal entries
#       }else{
#
#         # loglikelihood contributed by other entries
#         loglike = loglike+ dnorm(beta[i, j], mean = 0,
#                                  sd = sqrt(C0), log = TRUE)
#
#       }
#     }
#   }
#
#
#   #deal with Sigma next
#   for (i in 1:m){
#     loglike = loglike + log(densigamma(Sigma[i], alpha = nu/2, beta = nu*s^2/2))
#   }
#   return(loglike)
#
# }
