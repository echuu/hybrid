
# logML_approx.R

## functions in this file
## TODO: add function documentation/description
##
##     hybrid_ml_call()
##     hybml_const()
##     bridge_approx()
##



#### hybml_const_call() --------------------------------------------------------
#### this function is called from the candidatePartition.R file to fit the
#### first stage partition
#
hybml_const_call = function(u_df) {
    options(scipen = 999)
    options(dplyr.summarise.inform = FALSE)
    # const_vec  = numeric(N_approx) # store constant approximation
    D = ncol(u_df) - 1

    ## (2) fit the regression tree via rpart()
    u_rpart = rpart::rpart(psi_u ~ ., u_df)

    ## (3) process the fitted tree
    # (3.1) obtain the (data-defined) support for each of the parameters
    param_support = extractSupport(u_df, D) #

    # (3.2) obtain the partition
    u_partition = extractPartition(u_rpart, param_support)
    bounds_out = u_partition %>% dplyr::arrange(leaf_id) %>%
        dplyr::select(-c("psi_hat", "leaf_id"))

    # param_out = u_star_cand(u_rpart, u_df, u_partition, D) # partition.R
    param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R
    # opt_part = param_out$optimal_part

    # ----------------------------------------------------------------------
    n_partitions = nrow(u_partition) # number of partitions

    psi_partition = param_out %>%
        dplyr::select(-c('leaf_id', 'psi_choice', 'logQ_cstar', 'n_obs'))

    bounds = psi_partition %>% dplyr::select(-c("psi_star"))

    log_vol_vec = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>%
        log() %>% rowSums()

    logzhat = (-psi_partition$psi_star + log_vol_vec) %>% log_sum_exp

    return(list(logz   = logzhat,
                bounds = bounds_out))
}
# end of hybml_const_call() function -------------------------------------------





#### hybml_const():
#### Wrapper function that calls the main logml_call() function -- we include
#### this here to catch potential errors during simulations so that rare/
#### problematic draws don't cause the entire simulation to stop running.
#### Replications that result in an error will return NA, which is something we
#### check for in the final simulation results
hybml_const <- function(u_df) {
    out <- tryCatch(
        {
            hybml_const_call(u_df)
            # The return value of `readLines()` is the actual value
            # that will be returned in case there is no condition
            # (e.g. warning or error).
            # You don't need to state the return value via `return()` as code
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("hybrid computation error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            # message("Here's the original warning message:")
            # message(cond)
            return(NULL)
        },
        finally={
        }
    )
    return(out)
} # end of hybml_const() function ----------------------------------------------





bridge_approx <- function(samples, log_density, prior, lb, ub, method = 'normal') {
    out <- tryCatch(
        {

            bridge_result <- bridgesampling::bridge_sampler(samples = samples,
                                                            log_posterior = log_density,
                                                            data = prior,
                                                            lb = lb, ub = ub,
                                                            silent = TRUE,
                                                            method = method)
            bridge_result$logml
            # The return value of `readLines()` is the actual value
            # that will be returned in case there is no condition
            # (e.g. warning or error).
            # You don't need to state the return value via `return()` as code
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("bridge error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            # message("Here's the original warning message:")
            message(cond)
            return(out)
        },
        finally={
        }
    )
    return(out)
} # end of bridge_approx() function --------------------------------------------

