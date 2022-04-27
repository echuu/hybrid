


m  = 7        # number of columns in y, number of rows in beta
k0 = 1        # number of columns in beta (number of factors)
n  = 100      # number of rows in y (number of observations)

beta0 = c(0.995, 0.975, 0.949, 0.922, 0.894, 0.866, 0.837)
omega0 = diag(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
y = rmvnorm(n = n, mean = rep(0, m) , sigma = omega0 + tcrossprod(beta0))

# hyperparameters
nu = 2.2
C0 = 1
s = sqrt(.1/2.2)

# MCMC settings
thin    = 1
burnin  = 1000
seed    = 123
No.Iter = 2000    # number of iterations to run Gibbs sampler
n_post = No.Iter / thin


# dimensions that don't change
p        = m*(m+1)/2         ## number of elements in lower cholesky factor
K_mm     = matrixcalc::commutation.matrix(m, m)
ind_mat  = matrix(1:m^2, m, m)
chol_ind = ind_mat[lower.tri(ind_mat, diag = TRUE)]


#### draw samples from the prior of phi to find the sample mean/covariance
J = 5000 # number of prior samples to draw

# generate priors draws of sigma
sigma_prior = sample_sigma(J, m, k0) ## this will change as we vary k (# factors)

# obtain prior draws of phi by computing the lower chol factor of sigma
phi_list = lapply(sigma_prior, getLowerChol) # phi * phi' = sigma2^{-1}

# estimate the mean, covariance of the approximating mvn distribution
prior = estimate_params(phi_list)


set.seed(1)
n_sims = 1
k_vec = c(1:3)
pp = matrix(0, n_sims, length(k_vec))
pp_ez = matrix(0, n_sims, length(k_vec))
pp_lap = matrix(0, n_sims, length(k_vec))

for (i in 1:n_sims) {

  # (0) generate new parameters, data ------------------------------------------
  tau = 0.25
  omega0 = tau * diag(1, m) # identity matrix
  beta0 = matrix(0, m, k0)
  beta0[,1] = c(rtnorm(n = 1, lower = 0, mean = -1), rnorm(n = m - 1, mean = 1))
  # beta0[,2] = c(0, rtnorm(n = 1, lower = 0, mean = 1), rnorm(n = m - 2, mean = 1))
  # beta0[,3][3:m] = c(rtnorm(n= 1, lower = 0, mean = 1), rnorm(n = m - 3, mean = 1, sd = 1))
  y = rmvnorm(n = n, mean = rep(0, m) , sigma = omega0 + tcrossprod(beta0))
  S = t(y) %*% y

  ## fit model using a range for k = 2, 3, 4, 5

  # (1) fit models using different values of k; compute the corresponding logML
  for (j in 1:length(k_vec)) {

    # (1a) initialize variables that will be used in hybrid
    k = k_vec[j]
    # (1b) initialize parameter values for gibbs sampler
    beta_init = matrix(0, nr = m, nc = k)
    for(r in 1:k) {
      if(r == 1){
        beta_init[1, 1] = rtnorm(n = 1, mean = 0, sd = 1, lower = 0)
      } else {
        lower.vec = c(rep(-Inf, r - 1), 0)
        beta_init[r, 1:r] = rtmvnorm(n = 1, mean = rep(0, r), sigma = diag(r),
                                     lower = lower.vec)
      }
    }
    for(r in (k+1):m) {
      beta_init[r,] = rmvnorm(n = 1, mean = rep(0, k), sigma = diag(k))
    }
    Sigma_init = diag(omega0) ## init value for Sigma0 used in Gibbs sampler

    ## Gibbs sampler for the regular dataset: y
    obj = factGibbsMod(y = y, k = k, beta_init, Sigma_init, nu, C0, s,
                       No.Iter, burnin, thin)

    #### START: original samples not taking choleksy factor of sigma ###########
    # m = ncol(y)
    # d = (m-k) * k + k * (k+1) / 2
    # D = d + m
    # beta_samps = do.call(rbind, lapply(obj$beta, formatBeta))
    # omega_samps = do.call(rbind, obj$Sigma)
    # theta = data.frame(cbind(beta_samps, omega_samps))
    # u_df = hybrid::preprocess(data.frame(theta), D, z0)
    #### END: original samples not taking choleksy factor of sigma #############

    #### draw samples from the prior of phi to find the sample mean/covariance
    J = 10000 # number of prior samples to draw
    # generate priors draws of sigma
    sigma_prior = sample_sigma(J, m, k) ## this will change as we vary k (# factors)
    # obtain prior draws of phi by computing the lower chol factor of sigma
    phi_list = lapply(sigma_prior, getLowerChol) # phi * phi' = sigma2^{-1}
    # estimate the mean, covariance of the approximating mvn distribution
    prior = estimate_params(phi_list)

    z0 = list(y = y, S = S,
              n = n, nu = nu, C0 = C0, s = s, k = k, m = m, p = p,
              K_mm = K_mm, chol_ind = chol_ind,
              mu_phi = prior$mu,
              sigma_phi = prior$sigma)

    ## obtain posterior samples of phi
    phi_samps = phi_post(obj, n_post, p) # (n_post x p)

    ## evaluate objective function in terms of phi
    u_df = hybrid::preprocess(data.frame(phi_samps), p, z0)

    # (1d) compute approximation to log marg like for model k
    pp[i,j] = hybrid::hybml(u_df, z0, psi = psi, grad = grad, hess = hess)$logz
    pp_ez[i,j] = hybridml::hybml_const(u_df)$logz

    MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))[1]
    ustar = u_df[MAP_LOC,1:p] %>% unname() %>% unlist()
    pp_lap[i,j] = 0.5*(p)*log(2*pi) -
      0.5*hybrid::log_det(hess(ustar, z0)) - psi(ustar, z0)

  } ## end inner for() over k

  # print(paste('model prob for iter ', i, '/', n_sims, ': ', sep = ''))
  cat(paste('model prob for iter ', i, '/', n_sims, ': ', sep = ''))
  cat(pp[i,])
  cat('\n')
  cat(find_prob(pp[i,]))
  cat('\n')
} ## end outer for() over # of replications


testp = pp[1,]
find_prob = function(x) { round(exp(x - log_sum_exp(x)), 5) }
find_outcome = function(x) { round(exp(x - log_sum_exp(x))) }
model_probs = t(apply(pp, 1, find_prob))
model_choice = t(apply(pp_ez, 1, find_outcome))

model_probs
colMeans(model_choice)




