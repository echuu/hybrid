


## log posterior of x = (beta, Omega)
logpost = function(x, z) { # x = (lower.tri(beta), diag(Omega))
  # x is (m+d)-dimensional
  # beta elements are d-dimensional, d = (m-k) * k + k * (k+1) / 2
  # omega elements are m-dimensional
  d = z$d
  m = z$m
  k = z$k
  n = z$n
  y = z$y

  xbeta = x[1:d]
  xomega = x[(d+1):(m+d)]

  ## recreate beta matrix
  b = matrix(0, m, k)
  b[lower.tri(b, diag = TRUE)] = xbeta

  ## create Omega matrix
  Omega = diag(xomega)
  Sigma = Omega + b %*% t(b)
  Sigma_inv = solve(Sigma)

  q_term = 0
  for (i in 1:n) {
    q_term = q_term + c(t(y[i,]) %*% Sigma_inv %*% y[i,])
  }
  ## (1) compute loglikelihood -- this is the same value as ll_omega, ll_beta
  lp1 = -n*m/2 * log(2*pi) - n/2 * log(det(Sigma)) - 1/2 * q_term

  ## (2) compute log prior of beta
  lp2 = lp_beta(xbeta, z)

  ## (3) compute log prior of omega
  lp3 = lp_omega(xomega, z)

  return(lp1 + lp2 + lp3)
} # end logpost() function


## gradient of log posterior of x = (beta, Omega)
grad_logpost = function(x, z) {
  d = z$d
  m = z$m
  k = z$k
  xbeta = x[1:d]
  xomega = x[(d+1):(m+d)]

  ## recreate beta matrix
  beta = matrix(0, m, k)
  beta[lower.tri(beta, diag = TRUE)] = xbeta

  ## create Omega matrix
  Omega = diag(xomega)

  return(c(grad_logpost_beta(beta, Omega, z),   ## grad of log-post wrt beta
           grad_logpost_omega(beta, Omega, z))) ## grad of log-post wrt Omega
} # end grad_logpost() function


## 2nd order derivative for the cross terms of log posterior of
## x = (beta, Omega); this computes the upper right block matrix in the hessian
## matrix. THe output is a (d x m) matrix
dbeta_domega = function(beta, Omega, z) {
  d = z$d
  m = z$m
  k = z$k
  y = z$y
  n = z$n
  # xbeta = x[1:d]
  # xomega = x[(d+1):(m+d)]
  ## recreate beta matrix
  # beta = matrix(0, m, k)
  # beta[lower.tri(beta, diag = TRUE)] = xbeta
  ## create Omega matrix
  # Omega = diag(xomega)

  # Compute 2nd order mixed partial: (dBeta dOmega) --> (d x m) matrix returned
  # repeated calculations:
  Sigma = Omega + beta %*% t(beta)
  Sigma_inv = solve(Sigma)
  tmp = matrix(0, m*k, m*m)
  # SyyS = Sigma_inv %*% y %*% t(y) %*% Sigma_inv
  for (i in 1:n) {
    SyyS = Sigma_inv %*% y[i,] %*% t(y[i,]) %*% Sigma_inv
    tmp = tmp +
      ( ((t(beta) %*% Sigma_inv) %x% SyyS) +
          ((t(beta) %*% SyyS) %x% Sigma_inv) ) %*% z$K_mm
    # print(dim(tmp))
  }

  # hessian of log det term
  h_ld =  -2 * (t(beta) %*% Sigma_inv) %x% Sigma_inv
  h_q = 2 * tmp

  ## note: when we subset indices out of the 2nd order partial in matrix form
  ## we use beta_INDEX, omega_INDEX, rather than beta_ROWS, omega_ROWS
  ## b/c the mixed partial including beta and omega comes after all of the
  ## 2nd order partials wrt beta, so we started indexing from 1 to extract the
  ## relevant omega terms that correspond to the diagonal elements of omega
  # h_beta_omega = (-n/2 * h_ld - 1/2 * h_q)[beta_index, omega_index]
  # h[beta_rows, omega_rows] = h_beta_omega
  # h[omega_rows, beta_rows] = t(h_beta_omega)

  return((-n/2 * h_ld - 1/2 * h_q)[z$beta_index, z$omega_index])
} # end dbeta_domega() function


## hessian of log posterior of x = (beta, Omega) -> (m + d) x (m + d) matrix
hess_logpost = function(x, z) {

  d = z$d
  m = z$m
  k = z$k
  y = z$y
  n = z$n
  xbeta = x[1:d]
  xomega = x[(d+1):(m+d)]

  ## recreate beta matrix
  beta = matrix(0, m, k)
  beta[lower.tri(beta, diag = TRUE)] = xbeta
  ## create Omega matrix
  Omega = diag(xomega)

  # stack the hessian matrices along the diagonal to form block diagonal matrix
  h = magic::adiag(hess_logpost_beta(beta, Omega, z),
                   hess_logpost_omega(beta, Omega, z))

  # computed 2nd order mixed partials
  h[z$beta_rows, z$omega_rows] = dbeta_domega(beta, Omega, z)
  h[z$omega_rows, z$beta_rows] = t(h[z$beta_rows, z$omega_rows])

  return(h)
} # end hess_logpost() function

#### ---------------------------------------------------------------------------
#### log posterior functions in terms of the parameters: (beta, Omega) written
#### separately. These are eventually called by the highest level log posterior
#### functions

## log posterior of beta
logpost_beta = function(x, z) { ## this function just used to verify grad/hess
  lp_beta(x, z) + ll_beta(x, z)
} # end logpost_beta() function


# gradient of log posterior as function of beta
grad_logpost_beta = function(beta, Omega, z) {
  grad_logp_beta(beta, z) + grad_ll_beta(beta, Omega, z)
} # end grad_logpost_beta() function


## hessian of log posterior of beta
hess_logpost_beta = function(beta, Omega, z) {
  hess_logp_beta(beta, z) + hess_ll_beta(beta, Omega, z)
} # end hess_logpos_beta() function


## log posterior as function of omega
logpost_omega = function(x, z) {
  lp_omega(x, z) + ll_omega(x, z)
} # end of logpost_omega() function


## gradient of log posterior as function of omega
grad_logpost_omega = function(beta, Omega, z) {
  grad_logp_omega(Omega, z) + grad_ll_omega(beta, Omega, z)
} # end of grad_logpost_omega() function


## hessian of log posterior as function of omega
hess_logpost_omega = function(beta, Omega, z) {
  hess_logp_omega(Omega, z) + hess_ll_omega(beta, Omega, z)
} # end of hess_logpost_omega() function


##### log likelihood functions -------------------------------------------------

## log likelihood as a function of beta
ll_beta = function(x, z) { ## log likelihood as function of beta
  y = z$y
  n = z$n
  m = z$m
  k = z$k
  b = matrix(0, m, k)
  b[lower.tri(b, diag = TRUE)] = x
  Sigma = Omega + b %*% t(b) ## Omega is global here -- not ideal but ok for now
  Sigma_inv = solve(Sigma)
  q_term = 0
  for (i in 1:n) {
    q_term = q_term + c(t(y[i,]) %*% Sigma_inv %*% y[i,])
  }
  -n*m/2 * log(2*pi) -n/2 * log(det(Sigma)) - 1/2 * q_term
} # end ll_beta() function


## gradient of log likelihood wrt beta
grad_ll_beta = function(beta, Omega, z) {
  # b = matrix(0, z$m, z$k)
  # b[lower.tri(b, diag = TRUE)] = x
  y = z$y
  n = z$n
  m = z$m
  k = z$k
  Sigma = Omega + beta %*% t(beta)
  Sigma_inv = solve(Sigma)
  g1 =  2 * Sigma_inv %*% beta # grad for log-det term
  tmp = matrix(0, m, k)
  for (i in 1:n) {
    tmp = tmp + Sigma_inv %*% y[i,] %*% t(y[i,]) %*% Sigma_inv %*% beta
  }
  g2 = -2 * tmp # grad for quadratic term
  g = -n/2 * g1 - 1/2 * g2
  return(g[lower.tri(g, diag = TRUE)])
} # end grad_ll_beta() function


## hessian of the log likelihood wrt beta
hess_ll_beta = function(beta, Omega, z) {
  y = z$y
  n = z$n
  k = z$k
  m = z$m

  # construct the indices used to extract the free elements from beta
  # tmp_mat = matrix(1:(k*m), m, k)
  # b_ind = tmp_mat[lower.tri(tmp_mat, diag = TRUE)]
  b_ind = z$beta_index

  Sigma = Omega + beta %*% t(beta)
  Sigma_inv = solve(Sigma)
  # K_mk = commutation.matrix(m, k)
  K_mk = z$K_mk

  # hessian of log-det term
  h1 = -2 * (
    -(diag(1, z$k) %x% Sigma_inv) +
      ((t(beta) %*% Sigma_inv) %x% (Sigma_inv %*% beta)) %*% K_mk +
      (t(beta) %*% Sigma_inv %*% beta) %x% Sigma_inv
  )
  # hessian of quadratic term
  q_term = matrix(0, m*k, m*k)
  for (i in 1:n) {
    SyyS = Sigma_inv %*% y[i,] %*% t(y[i,]) %*% Sigma_inv
    q_term = q_term + t(K_mk) %*% ((Sigma_inv %*% beta) %x% (t(beta) %*% SyyS) +
                                     (SyyS %*% beta) %x% (t(beta) %*% Sigma_inv)) +
      (t(beta) %*% Sigma_inv %*% beta) %x% SyyS +
      (t(beta) %*% SyyS %*% beta) %x% Sigma_inv -
      (diag(1, k) %x% SyyS)
    # print(dim(q_term))
  }
  # h2 = 2 * (
  #   t(K_mk) %*% ((Sigma_inv %*% beta) %x% (t(beta) %*% SyyS) +
  #                  (SyyS %*% beta) %x% (t(beta) %*% Sigma_inv)) +
  #     (t(beta) %*% Sigma_inv %*% beta) %x% SyyS +
  #     (t(beta) %*% SyyS %*% beta) %x% Sigma_inv -
  #     (diag(1, z$k) %x% SyyS)
  # )
  h2 = 2 * q_term
  h = -n/2 * h1 -1/2 * h2
  return(h[b_ind, b_ind])
} # end hess_ll_beta() function


## log likelihood as a function of Omega
ll_omega = function(x, z) { ## log likelihood as function of omega (missing constants)
  y = z$y
  n = z$n
  m = z$m
  k = z$k
  Sigma = (diag(x) + beta %*% t(beta))
  Sigma_inv = solve(Sigma)
  q_term = 0
  for (i in 1:n) {
    q_term = q_term + c(t(y[i,]) %*% Sigma_inv %*% y[i,])
  }
  -n*m/2 * log(2*pi) -n/2 * log(det(Sigma)) - 1/2 * q_term
} # end ll_omega() function


## gradient of log likelihood wrt omega
grad_ll_omega = function(beta, Omega, z) {
  y = z$y
  n = z$n
  m = z$m
  k = z$k
  Sigma = Omega + beta %*% t(beta)
  Sigma_inv = solve(Sigma)
  g1 = Sigma_inv                                    # log-det term

  q_term = matrix(0, m, m)
  for (i in 1:n) {
    q_term = q_term + c(Sigma_inv %*% y[i,] %*% t(y[i,]) %*% Sigma_inv)
  }
  # g2 = -Sigma_inv %*% z$y %*% t(z$y) %*% Sigma_inv  # quadratic term
  g2 = -q_term
  return(diag(-n/2 * g1 - 1/2 * g2))
} # end grad_ll_omega() function


## hessian of log likelihood wrt omega
hess_ll_omega = function(beta, Omega, z) {
  y = z$y
  n = z$n
  m = z$m
  k = z$k
  h_index = (1:z$m - 1) * (z$m + 1) + 1 # for 0-index in C++, remove the (+1)
  Sigma = (Omega + beta %*% t(beta))
  Sigma_inv = solve(Sigma)
  h1 = - Sigma_inv %x% Sigma_inv # hess of log-det term
  # hess of quadratic term
  q_term = matrix(0, m^2, m^2)
  for (i in 1:n) {
    q_term = q_term + (Sigma_inv %*% y[i,] %*% t(y[i,]) %*% Sigma_inv) %x% Sigma_inv
  }
  h2 = 2 * q_term
  h = -n/2 * h1 -1/2 * h2
  return(h[h_index, h_index])
} # end hess_ll_omega() function

##### log prior functions ------------------------------------------------------

## log prior of beta
lp_beta_old = function(x, z) { # still takes original form of input (easier)
  m = z$m; k = z$k; C0 = z$C0;
  beta = matrix(0, m, k)
  beta[lower.tri(beta, diag = TRUE)] = x
  logprior = 0
  for (i in 1:m){ # iterate through beta matrix row-wise
    for (j in 1: min(k, i)) {  # j is column index
      if (i == j) { # loglikelihood contributed by diagonal entries
        logprior = logprior + dchisq(1/C0*beta[i,j]^2, df = k-i+1, log = TRUE) +
          log(2) + log(beta[i,j])- log(C0)
      } else { # loglikelihood contributed by other entries
        logprior = logprior + dnorm(beta[i,j], 0, sqrt(C0), log = TRUE)
      }
    }
  }
  return(logprior)
} # end lp_beta() function


lp_beta = function(x, z) {
  m = z$m; k = z$k; C0 = z$C0;
  beta = matrix(0, m, k)
  beta[lower.tri(beta, diag = TRUE)] = x
  logprior = 0
  for (i in 1:m){ # iterate through beta matrix row-wise
    for (j in 1: min(k, i)) {  # j is column index
      if (i == j) { # loglikelihood contributed by diagonal entries
        logprior = logprior - lgamma(0.5 * (k-i+1)) -
          (0.5 * (k-i+1) - 1) * log(2) - 0.5 * (k-i) * log(C0) +
          (k-i) * log(beta[i,i]) - beta[i,i]^2 / (2 * C0)
      } else { # loglikelihood contributed by other entries
        logprior = logprior - 0.5 * log(2*pi*C0) - 1 / (2 * C0) * beta[i,j]^2
      }
    }
  }
  return(logprior)
}


## gradient of log prior of beta
grad_logp_beta = function(beta, z) {
  gvec = matrix(0, z$m, z$k)
  for (i in 1:z$m){
    for(j in 1: min(z$k, i)){
      if(i == j){
        gvec[i,j] = -1/z$C0 * beta[i, j] + (z$k-i) / beta[i, j]
      } else {
        gvec[i,j] = -1/z$C0 * beta[i, j]
      }
    }
  }
  return(gvec[lower.tri(gvec, diag = TRUE)])
} # end grad_logp_beta() function


## hessian of log prior of beta
hess_logp_beta = function(beta, z) {
  hmat = matrix(0, z$d, z$d)
  for (r in 1:z$d) {
    i = z$ind_mat[r, 1]; j = z$ind_mat[r, 2];
    for (c in r:z$d) {
      kk = z$ind_mat[c, 1]; ll = z$ind_mat[c, 2];
      if ((i == j) & (kk == ll) & (i == kk)) {
        # print(paste('(', i, ', ', j, ')', ', ',
        #             '(', kk, ', ', ll, ')',sep = ''))
        hmat[r, c] = -1/z$C0  - (z$k-i) / beta[i,j]^2
        # print(-1/C_0  - (k-i) / b[i,j]^2)
      } else if ((i == kk) & (j == ll)) {
        hmat[r, c] = -1/z$C0
      }
    }
  }
  return(hmat)
} # end of hess_logp_beta() function


## log prior of omega
lp_omega = function(x, z) { # still takes original form of input (easier)
  loglike = 0
  m = z$m
  for (i in 1:m){
    loglike = loglike + log(densigamma(x[i],
                                       alpha = z$nu / 2,
                                       beta  = z$nu * z$s^2 / 2))
  }
  return(loglike)
} # end of lp_omega() function


## gradient of log prior of Omega
grad_logp_omega = function(Omega, z) {
  gvec = numeric(z$m)
  for (i in 1:z$m) {
    gvec[i] = 1/Omega[i,i] * ( 1/Omega[i,i] * z$nu*z$s^2/2 -(z$nu/2+1) )
  }
  return(gvec)
} # end of grad_logp_omega() function


## hessian of log prior of Omega
hess_logp_omega = function(Omega, z) {
  h = matrix(0, z$m, z$m)
  for (i in 1:z$m) {
    h[i,i] = 1/Omega[i,i]^2 * ( -1/Omega[i,i] * z$nu*z$s^2 + (z$nu/2+1) )
  }
  return(h)
} # end of hess_logp_omega() function
