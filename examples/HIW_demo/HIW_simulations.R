
## HIW_simulations.R
## author: Eric Chuu

## This file compares the performance of the GNORM and EPSOM-Hyb estimators.
## for high-dimensional graphs that are constructed using G_9, as described
## in Section 4.1 of the paper. Since G_9 is decomposable, the resulting
## true normalizing constant is available in closed form, so we evaluate 
## the estimates against the ground truth. 


# load dplyr for easier data manipulation
library(dplyr)

# load Atay's GNORM estimator
library(BDgraph)

# load EPSOM-Hyb estimator
library(hybrid)

# important: set working directory to source file location
setwd("C:/Users/ericc/Dropbox/jcgs_simulation_files")


# load HIW density related functions
source("HIW_helper.R")

# load faster gradient, hessian functions (C++ implementations)
Rcpp::sourceCpp("hiw.cpp")


# create G_9 as given in the paper
G_9 = matrix(  c(1,1,0,0,1,0,0,0,0,
                 1,1,1,1,1,0,0,0,0,
                 0,1,1,1,0,0,0,0,0,
                 0,1,1,1,1,1,1,0,0,
                 1,1,0,1,1,1,0,0,0,
                 0,0,0,1,1,1,1,1,1,
                 0,0,0,1,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1), 9, 9)
a = c(1, 3, 2, 5, 4, 6, 7, 8, 9)
G_9 = G_9[a, a]


################################################################################

n_G9 = 8                    # number of times to stack G_9 along block diagonal
G = diag(1, n_G9) %x% G_9   # larger graph created from G_9
p = ncol(G)                 # dimension of vertex set
V = diag(1, p)              # prior scale matrix
b = 3                       # prior degrees of freedom
D_0 = 0.5 * p * (p + 1)     # number of entries on diagonal and upper diagonal
J = 1000                    # number of MCMC samples to draw

# logical vector determining existence of edges between vertices
edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical
upperInd = G[upper.tri(G)] %>% as.logical
D_u = sum(edgeInd)

#### data generation step
set.seed(1)
true_params = HIWsim(G, b, V)
Sigma_G = true_params$Sigma
Omega_G = true_params$Omega # true precision matrix

# Generate data Y based on Sigma_G
N = 100
Y = matrix(0, N, p)
for (i in 1:N) {
  Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(p, 0, 1)) # (500 x D)
}
S = t(Y) %*% Y


# Compute HIW-related quantities
nu = rowSums(chol(Omega_G) != 0) - 1
xi = b + nu - 1
t_ind = which(chol(Omega_G) != 0, arr.ind = T)

# initialize HIW parameters
params = list(N = N, D = p, D_0 = D_0, D_u = D_u,
              testG = G, edgeInd = edgeInd,
              upperInd = upperInd, S = S, V = V, b = b,
              nu = nu, xi = xi, G = G, t_ind = t_ind)

################################################################################

#### Start EPSOM-Hyb approximation

# sample from posterior distribution
postIW = sampleHIW(J, D_u, D_0, G, b, N, V, S, edgeInd)
post_samps = postIW$post_samps                 # J x D_u

# evaluate the posterior samples using negative log posterior
u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)

# compute the global mode of the target distribution
u_star = hybrid::globalMode(u_df, params)

# compute the EPSOM-Hyb estimator
out = hybrid::hybml(u_df, params, grad = grad, hess = hess, u_0 = u_star)
out$logz

# compute the true log marginal likelihood
(LIL = logmarginal(Y, G, b, V, S))

# compute Atay's GNORM approximation to the log marginal likelihood
- 0.5 * p * N * log(2 * pi) + BDgraph::gnorm(G, b + N, V + S, J) -
  BDgraph::gnorm(G, b, V, J)



