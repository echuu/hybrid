

# load density_wrapper.R


#### true parameter settings

m  = 7        # number of columns in y, number of rows in beta
k0 = 1        # number of columns in beta (number of factors)
n  = 500      # number of rows in y (number of observations)


#### generate true value of beta
# beta0 = matrix(0, nr = m, nc = k0) # init value for beta0 used in Gibbs sampler
# for(r in 1:k0) {
#   if(r == 1){
#     beta0 [1, 1] = rtnorm(n = 1, mean = 0, sd = 1, lower = 0)
#   } else {
#     lower.vec = c(rep(-Inf, r - 1), 0)
#     beta0[r, 1:r] = rtmvnorm(n = 1, mean = rep(0, r), sigma = diag(r),
#                              lower = lower.vec)
#   }
#
# }
# for(r in (k0+1):m) {
#   beta0[r,] = rmvnorm(n = 1, mean = rep(0, k0), sigma = diag(k0))
# }
# beta0[,1] = c(0.995, 0.975, 0.949, 0.922, 0.894, 0.866, 0.837)

beta0 = c(0.995, 0.975, 0.949, 0.922, 0.894, 0.866, 0.837)
omega0 = diag(c(0.01, 0.05, 0.01, 0.05, 0.02, 0.03, 0.30))

set.seed(123)
y = rmvnorm(n = n, mean = rep(0, m) , sigma = omega0 + tcrossprod(beta0)) # n x m


# fact_fit = FactorAnalyses(numCovariates=ncol(y), maxNumFactors=3)
# sbic_results = sBIC::sBIC(scale(diff(y)), fact_fit)
# sbic_results$sBIC



set.seed(1)
n_sims = 100
k_vec = c(1:3)
pp = matrix(0, n_sims, length(k_vec))
for (i in 1:n_sims) {

  # (0) generate new parameters, data ------------------------------------------
  # true.beta3 = matrix(0, nr = m, nc = k0)
  # true.beta3[, 1] =  c(rtnorm(n = 1, lower = 0, mean = -6), rnorm(n = m -1, mean = 9))
  # true.beta3[, 2] = c(0, rtnorm(n= 1, lower = 0, mean = 5), rnorm(n = m - 2, mean = 3))
  # true.beta3[, 3][3:m] = c(rtnorm(n= 1, lower = 0, mean = 2), rnorm(n = m - 3, mean = 5, sd = 18))
  # set.seed(123)
  # y3 = rmvnorm (n = n, mean= rep(0,m) , sigma = diag(m) + tcrossprod(true.beta3))
  # y3 = scale(y3, center = TRUE,  scale = TRUE)
  # y3.swap = y3[, sample(m)]

  y = rmvnorm(n = n, mean = rep(0, m) , sigma = omega0 + tcrossprod(beta0)) # n x m

  ## fit model using a range for k = 2, 3, 4, 5
  nu  = 2.2
  C0  = 2
  s   =  sqrt(.1/2.2)
  C0  = 1; C_0 = 1;

  # (1) fit models using different values of k; compute the corresponding logML
  for (j in 1:length(k_vec)) {

    # (1a) initialize variables that will be used in hybrid
    k = k_vec[j]
    m = ncol(y)
    d = (m-k) * k + k * (k+1) / 2
    # K_mk = commutation.matrix(m, k)
    # K_mm = commutation.matrix(m, m)
    ind_mat = matrix(1:(m*k), m, k)
    lowerind = ind_mat[lower.tri(ind_mat, diag = TRUE)]
    ind_mat = cbind((lowerind-1) %% m, floor((lowerind-1) / m) ) + 1
    omega_index = (1:m - 1) * (m + 1) + 1
    beta_rows = 1:d
    omega_rows = (d+1):(d+m)
    z0 = list(y = y, n = nrow(y), nu = nu, C0 = C0, s = s, k = k, m = m, d = d,
              ind_mat = ind_mat, K_mk = K_mk, K_mm = K_mm, beta_index = lowerind,
              omega_index = omega_index, beta_rows = beta_rows,
              omega_rows = omega_rows)


    # (1b) initialize parameter values for gibbs sampler
    beta0 = matrix(0, nr = m, nc = k) # init value for beta0 used in Gibbs sampler
    for(r in 1:k) {
      if(r == 1){
        beta0 [1, 1] = rtnorm(n = 1, mean = 0, sd = 1, lower = 0)
      } else {
        lower.vec = c(rep(-Inf, r - 1), 0)
        beta0[r, 1:r] = rtmvnorm(n = 1, mean = rep(0, r), sigma = diag(r),
                                 lower = lower.vec)
      }

    }
    for(r in (k+1):m) {
      beta0[r,] = rmvnorm(n = 1, mean = rep(0, k), sigma = diag(k))
    }
    Sigma0  = rep(1, m) ## init value for Sigma0 used in Gibbs sampler

    # (1c) obtain MCMC samples of (beta, Omega) for model_k
    thin    = 1
    burnin  = 1000
    seed    = 123
    No.Iter = 5000
    ## pass in initial values: (beta0, Sigma0)
    obj = factGibbsMod(y = y, k = k, beta0, Sigma0, nu, C0, s,
                       No.Iter, burnin, thin, seed)

    D = d + m
    beta_samps = do.call(rbind, lapply(obj$beta, formatBeta))
    omega_samps = do.call(rbind, obj$Sigma)
    theta = data.frame(cbind(beta_samps, omega_samps))

    # (1d) compute approximation to log marg like for model k
    u_df = hybrid::preprocess(theta, D, z0)
    pp[i,j] = hybridml::hybml_const(u_df)$logz

  } ## end inner for() over k

  # print(paste('model prob for iter ', i, '/', n_sims, ': ', sep = ''))
  cat(paste('model prob for iter ', i, '/', n_sims, ': ', sep = ''))
  cat(find_prob(pp[i,]))
  cat('\n')

} ## end outer for() over # of replications


testp = pp[1,]
find_prob = function(x) { round(exp(x - log_sum_exp(x)), 5) }
find_outcome = function(x) { round(exp(x - log_sum_exp(x))) }
model_probs = t(apply(pp, 1, find_prob))

model_choice = t(apply(pp, 1, find_outcome))
colMeans(model_choice)

