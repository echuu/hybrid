
### updated version of the hybrid approximation as of 3/2/21

### final form of the hybrid approximation function



### the following functions need to be passed into the function
### (1) psi     :  negative log posterior
### (2) grad    :  gradient of the negative log posterior
### (3) hess    :  hessian of the negative log posterior

hybml = function(u_df, params, psi, grad, hess, u_0 = NULL, D = ncol(u_df) - 1) {

    options(scipen = 999)
    options(dplyr.summarise.inform = FALSE)

    ## fit the regression tree via rpart()
    u_rpart = rpart::rpart(psi_u ~ ., u_df)

    ## (3) process the fitted tree
    # (3.1) obtain the (data-defined) support for each of the parameters
    param_support = extractSupport(u_df, D) #

    # (3.2) obtain the partition
    u_partition = extractPartition(u_rpart, param_support)

    # return(u_partition)

    #### hybrid extension begins here ------------------------------------------

    ### (1) find global mean
    # u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean

    if (is.null(u_0)) {
        MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
        u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
        # print(u_0)
    }

    ### (2) find point in each partition closest to global mean (for now)
    # u_k for each partition
    u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

    l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
    u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

    # take min result, group_by() leaf_id
    psi_df = u_df_part %>%
        dplyr::group_by(leaf_id) %>% dplyr::filter(l1_cost == min(l1_cost)) %>%
        data.frame

    bounds = u_partition %>% dplyr::arrange(leaf_id) %>%
        dplyr::select(-c("psi_hat", "leaf_id"))
    psi_df = psi_df %>% dplyr::arrange(leaf_id)

    K = nrow(bounds)
    log_terms = numeric(K) # store terms so that we can use log-sum-exp()
    G_k = numeric(K)       # store terms coming from gaussian integral

    # lambda_k = apply(psi_df[,1:D], 1, lambda, params = params)
    # k = 1
    for (k in 1:K) {
        u_k = unname(unlist(psi_df[k,1:D]))

        # print(u_k)

        H_k = hess(u_k, params = params)
        H_k = as.matrix(Matrix::forceSymmetric(H_k))
        # H_k_inv = suppressWarnings(chol2inv(chol(H_k))) ## for gwish example
        H_k_inv = solve(H_k)
        H_k_inv = as.matrix(Matrix::forceSymmetric(H_k_inv))

        # print(paste("k = ", k, "/", K, " - Hessian matrix symmetric: ",
        #     matrixcalc::is.symmetric.matrix(H_k_inv), sep = ''))

        # lambda_k = pracma::grad(psi, u_k, params = params) # numerical
        lambda_k = grad(u_k, params)
        b_k = H_k %*% u_k - lambda_k
        m_k = H_k_inv %*% b_k

        lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
        ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

        # G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
        # print('using TN:')
        # log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])
        G_k[k] = ep(m_k, H_k_inv, lb, ub)
        # G_k[k] = epmgp_stable(m_k, H_k_inv, b_k, lb, ub)$logZ
        # print(G_k[k])
        # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

        log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
            psi_df$psi_u[k] + sum(lambda_k * u_k) -
            0.5 * t(u_k) %*% H_k %*% u_k +
            0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]
    }

    return(list(logz = log_sum_exp(log_terms),
                bounds = bounds,
                G_k = G_k))
}


hybml_map = function(u_df, u_0 = NULL, D = ncol(u_df) - 1) {
    options(scipen = 999)
    options(dplyr.summarise.inform = FALSE)

    ## fit the regression tree via rpart()
    u_rpart = rpart::rpart(psi_u ~ ., u_df)

    ## (3) process the fitted tree
    # (3.1) obtain the (data-defined) support for each of the parameters
    param_support = extractSupport(u_df, D) #

    # (3.2) obtain the partition
    u_partition = extractPartition(u_rpart, param_support)

    # return(u_partition)

    #### hybrid extension begins here ------------------------------------------

    ### (1) find global mean
    # u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean

    if (is.null(u_0)) {
        MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
        u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
        # print(u_0)
    }

    ### (2) find point in each partition closest to global mean (for now)
    # u_k for each partition
    u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

    l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
    u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

    # take min result, group_by() leaf_id
    psi_df = u_df_part %>%
        dplyr::group_by(leaf_id) %>% dplyr::filter(l1_cost == min(l1_cost)) %>%
        data.frame

    bounds = u_partition %>% dplyr::arrange(leaf_id) %>%
        dplyr::select(-c("psi_hat", "leaf_id"))
    psi_df = psi_df %>% dplyr::arrange(leaf_id)


    log_vol_vec = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>%
        log() %>% rowSums()

    logzhat = (-psi_df$psi_u + log_vol_vec) %>% log_sum_exp

    return(list(logz   = logzhat,
                bounds = bounds))
}


globalMode = function(u_df, params, D = ncol(u_df) - 1, tolerance = 0.00001,
                      maxsteps = 200, factormodel = TRUE) {

    # use the MAP as the starting point for the algorithm
    MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))[1]
    theta = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()

    numsteps = 0
    tolcriterion = 100
    step.size = 1


    while(tolcriterion > tolerance && numsteps < maxsteps){

        G = -hess(theta, params)
        invG = solve(G)
        thetaNew = theta + step.size * invG %*% grad(theta, params)

        # print(thetaNew)
        if (factormodel && any(thetaNew[params$diag_index] <= 0)) {
            ## idea is to check if diagonal entries are positive
            ## if not, then don't make the move
            ## try option of positive version of the value to see if it's a good move
            print("cannot have negative entry on diagonal of beta")
            # thetaNew[params$diag_index] = abs(thetaNew[params$diag_index])
            # print(paste("post prob of update:", -psi(thetaNew, params)))
            # print(paste("post prob of old:", -psi(theta, params)))
            thetaNew = theta # do not update theta
        }

        # if precision turns negative or if the posterior probability of
        # thetaNew becomes smaller than the posterior probability of theta
        if(-psi(thetaNew, params) < -psi(theta, params)) {
            cat('tolerance reached on log scale =', tolcriterion, '\n')
            print(paste("converged -- ", numsteps, " iters", sep = ''))
            return(theta)
        }

        tolcriterion = abs(psi(thetaNew, params)-psi(theta, params))
        theta = thetaNew
        numsteps = numsteps + 1
    }

    if(numsteps == maxsteps)
        warning('Maximum number of steps reached in Newton method.')

    print(paste("converged -- ", numsteps, " iters", sep = ''))
    return(theta)
}



ep_step = function(lb, ub, m_k, H_k_inv) {
    out <- tryCatch(
        {
            epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
        },
        error=function(cond) {
            message("integral evaluation is 0")
            message(cond)
            return(0)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            return(0)
        },
        finally={
        }
    )
    return(out)}
