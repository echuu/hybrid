
# partition.R


## functions in this file
## TODO: add function documentation/description
##
##     extractPartition()
##     u_star()
##     updatePartition()
##     node_partition()
##

#### extractPartition() --------------------------------------------------------
#
extractPartition = function(u_tree, param_support = NULL) {

    # pre-process the partition output
    rules = rpart.plot::rpart.rules(u_tree, digits = 6, roundint = FALSE)

    ## 1/23 -- updated following 3 lines
    rules_df_full = apply(rules, 2, as.character)
    # psi_hat_rules = round(as.numeric(rules_df_full[,1]), 4)
    rules_df = rules_df_full[,-c(1:2)]
    rules_str = data.frame(x = matrix(0, nrow(rules_df))) # n_partitions x 1

    partition_id = sort(unique(u_tree$where)) # row id of leaf node information

    # this is now in the order of the rules_df_full, rules_df dataframes
    psi_hat_id = cbind(leaf_id = partition_id,
                       psi_hat = u_tree$frame[partition_id, 5]) %>%
        data.frame() %>% dplyr::arrange(psi_hat)

    # psi_hat = round(u_tree$frame[partition_id,]$yval, 6)
    # psi_hat_sort = sort(psi_hat) # sort these to match the order of the rules

    n_params = length(u_tree$ordered)
    n_partitions = nrow(rules_str)


    for(r in 1:nrow(rules)) {
        rules_str[r,] = stringr::str_c(rules_df[r,][rules_df[r,] != ""],
                                       collapse = ' ')
    }


    # initialize storage for the partition intervals -- (n_params x 2)
    partition = data.frame(matrix(NA, n_partitions, n_params * 2))

    for (p in 1:n_params) {

        p_col = 2 * p - 1

        # form column names for {(u1_lb, u1_ub),...,(up_lb, up_ub)}
        lb = paste("u", p, "_lb", sep = "")
        ub = paste("u", p, "_ub", sep = "")
        names(partition)[p_col:(p_col + 1)] = c(lb, ub)

        # fill in the support for the p-th parameter
        partition[, p_col] = param_support[p, 1]       # param p lower bound
        partition[, p_col + 1] = param_support[p, 2]   # param p upper bound
    }


    # populate the partition with decision rules
    for (r in 1:n_partitions) {

        # (0)   split the string by & (each decision stored in matrix/vector)
        part_vec = stringr::str_split(rules_str[r,], "\\&", simplify = TRUE)

        # (1)   iterate over the partition/decisions
        for (p in 1:ncol(part_vec)) {


            ### fixed : node_partition() function is buggy
            processNode = node_partition(part_vec[p], param_support)

            # (1.1) identify the parameter corresponding to each partition
            param_id = processNode$id
            col_id = 2 * param_id - 1

            # (1.2) extract values defining the partition/decision
            bdry = processNode$interval

            # (1.3) store the values from (1.2) into correct column
            partition = updatePartition(partition, r, col_id, bdry)

        } # end of update step for each of the parameters

    } # end of loop storing the partition boundaries

    ## 1/23 -- updated the left-appended column
    partition_out = cbind(psi_hat_id, partition)



    # number of observations in each leaf node
    # part_obs_tbl = table(u_tree$where) %>% data.frame
    # names(part_obs_tbl) = c("leaf_id", "n_obs")

    #### (2) obtain predicted value for each of the observations
    # psi_hat_leaf = cbind(leaf_id = partition_id,
    #                     psi_hat = u_tree$frame[partition_id,]$yval) %>%
    #    data.frame()

    #partition_out = merge(psi_hat_leaf %>%
    #                          mutate(psi_hat = round(psi_hat, 4)),
    #                      partition, by = 'psi_hat')

    return(partition_out)

}
# end extractPartition() function ----------------------------------------------





#### plotPartition() -----------------------------------------------------------
# plot the 2-d partition with psi labels
# highest 3 value points are plotted to show areas of low probability (?)
# TODO: modify function to plot highest 3 value points rather than just the
# median; can plot the median value as a separate color
#
plotPartition = function(u_df, param_out) {

    plot(u_df[,1], u_df[,2], pch = 20, cex = 1, col = "cyan",
         xlab = 'u1', ylab = 'u2', main = '')
    rect(param_out$u1_lb,
         param_out$u2_lb,
         param_out$u1_ub,
         param_out$u2_ub)

    # add psi_hat labels for each partition
    text(x = param_out$u1_lb + (param_out$u1_ub - param_out$u1_lb) / 2,
         y = param_out$u2_lb + (param_out$u2_ub - param_out$u2_lb) / 2,
         labels = round(param_out$psi_hat, 5),
         cex = 0.8)

    # make the 'median' points red and large
    # points(x = param_out$u1_star, y = param_out$u2_star,
    #        col = 'red', pch = 19, cex = 1.2)

}
# end plotPartition() function -------------------------------------------------




#### rpart_mse() ---------------------------------------------------------------
# compute partition-wise MSE using the fitted values provided by the rpart
# regression tree; each mse can be extracted using a leaf_id
#
rpart_mse = function(rpart_obj, u_df_in) {

    # (1.1) determine which partition each observation is grouped in
    u_df = u_df_in %>% dplyr::mutate(leaf_id = rpart_obj$where)

    # (1.2) obtain the rows in rpart.object$frame associated with the leaf nodes
    # these rows contain the fitted value for nodes that fall w/in a
    # given partition
    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    # number of observations in each leaf node
    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    #### (2) obtain predicted value for each of the observations
    psi_hat_leaf = cbind(leaf_id = partition_id,
                         psi_hat = rpart_obj$frame[partition_id,]$yval) %>%
        data.frame()

    # append the fitted value for each row on the right as an additional column
    u_df = merge(u_df, psi_hat_leaf, "leaf_id")

    #### (3) compute squared residuals for each of the observations
    u_df = u_df %>% dplyr::mutate(dev_sq = (psi_u - psi_hat)^2)

    return(u_df)


} # end rpart_mse() function ---------------------------------------------------




#### u_star() ------------------------------------------------------------------
# updated 4/23 - choice of psi_star is the one that minimizes log(Q(c))
# output features: leaf_id, psi_choice (max, median, mean, hat), psi_star,
# logQ_cstar, n_obs (# of samples that partition)
#
u_star_max = function(rpart_obj, u_df, partition, n_params) {


    u_df = u_df %>% dplyr::mutate(leaf_id = rpart_obj$where)

    # extract fitted values from tree, merged with other psis later
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    n_partitions = length(partition_id)

    ### start updated code 4/27
    psi_all = u_df %>% dplyr::group_by(leaf_id) %>%
        dplyr::summarise(psi_max  = max(psi_u))

    psi_long = reshape2::melt(psi_all, id.vars = c("leaf_id"),
                              value.name = "psi_star",
                              variable.name = "psi_choice")

    psi_all_df = psi_long %>%
        dplyr::mutate(logQ_cstar = 0)


    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>%
            dplyr::mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                              sapply(psi_star, logQ, c_k = c_k),
                                              logQ_cstar))

    } # end of loop extracting representative points

    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_all_df, partition, by = 'leaf_id') %>%
        dplyr::select(-psi_hat)

    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

    return(partition_star)

}
# end of u_star() function -----------------------------------------------------





#### u_star() ------------------------------------------------------------------
# updated 4/23 - choice of psi_star is the one that minimizes log(Q(c))
# output features: leaf_id, psi_choice (max, median, mean, hat), psi_star,
# logQ_cstar, n_obs (# of samples that partition)
#
u_star_min = function(rpart_obj, u_df, partition, n_params) {


    u_df = u_df %>% dplyr::mutate(leaf_id = rpart_obj$where)

    # extract fitted values from tree, merged with other psis later
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    n_partitions = length(partition_id)

    ### start updated code 4/27
    psi_all = u_df %>% dplyr::group_by(leaf_id) %>%
        summarise(psi_min  = min(psi_u))

    psi_long = melt(psi_all, id.vars = c("leaf_id"), value.name = "psi_star",
                    variable.name = "psi_choice")

    psi_all_df = psi_long %>%
        dplyr::mutate(logQ_cstar = 0)


    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>%
            mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                       sapply(psi_star, logQ, c_k = c_k),
                                       logQ_cstar))

    } # end of loop extracting representative points

    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_all_df, partition, by = 'leaf_id') %>%
        dplyr::select(-psi_hat)

    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

    return(partition_star)

}
# end of u_star() function -----------------------------------------------------




#### u_star() ------------------------------------------------------------------
# updated 4/23 - choice of psi_star is the one that minimizes log(Q(c))
# output features: leaf_id, psi_choice (max, median, mean, hat), psi_star,
# logQ_cstar, n_obs (# of samples that partition)
#
u_star = function(rpart_obj, u_df, partition, n_params) {


    u_df = u_df %>% dplyr::mutate(leaf_id = rpart_obj$where)

    # extract fitted values from tree, merged with other psis later
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    n_partitions = length(partition_id)

    # compute max for each of the partitions
    # psi_all = u_df %>% dplyr::group_by(leaf_id) %>%
    #     summarise(psi_max  = max(psi_u),
    #               psi_med  = median(psi_u),
    #               psi_mean = mean(psi_u),
    #               psi_85   = quantile(psi_u, 0.85),
    #               psi_90   = quantile(psi_u, 0.90),
    #               psi_95   = quantile(psi_u, 0.95)) %>%
    #     merge(psi_hat_df, by = 'leaf_id')

    ### start updated code 4/27
    # psi_center = u_df %>% dplyr::group_by(leaf_id) %>%
    #     summarise(psi_med  = median(psi_u),
    #               psi_mean = mean(psi_u),
    #               psi_min  = min(psi_u)) %>%
    #     merge(psi_hat_df, by = 'leaf_id')

    # psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
    #    do(data.frame(t(quantile(.$psi_u, probs = seq(0.5, 1, 0.05)))))

    psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
        do(data.frame(t(quantile(.$psi_u, probs = seq(0.995, 1, 0.001)))))

    names(psi_quant) = c("leaf_id", paste('psi_', seq(95, 100, 1), sep = ''))

    # psi_all = merge(psi_center, psi_quant, by = 'leaf_id')
    psi_all = psi_quant
    ### end updated code 4/27

    psi_long = reshape2::melt(psi_all, id.vars = c("leaf_id"),
                              value.name = "psi_star",
                              variable.name = "psi_choice")

    psi_all_df = psi_long %>%
        dplyr::mutate(logQ_cstar = 0)


    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>%
            dplyr::mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                              sapply(psi_star, logQ, c_k = c_k),
                                              logQ_cstar))
    } # end of loop extracting representative points

    # for each partition (leaf_id), subset out rows for which log(Q(c)) is min
    psi_df = psi_all_df %>%
        dplyr::group_by(leaf_id) %>%
        dplyr::slice(which.min(logQ_cstar)) %>%  # extract rows that min log(Q(c))
        data.frame()

    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_df, partition, by = 'leaf_id') %>%
        dplyr::select(-psi_hat)

    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

    return(partition_star)

}
 # end of u_star() function ----------------------------------------------------


### 8/25 : this is the current version of u_star
u_star_cand = function(rpart_obj, u_df, partition, n_params) {

    u_df = u_df %>% mutate(leaf_id = rpart_obj$where)

    # extract fitted values from tree, merged with other psis later
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    n_partitions = length(partition_id)

    psi_center = u_df %>% dplyr::group_by(leaf_id) %>%
        summarise(psi_max  = max(psi_u))

    # psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
    #     do(data.frame(t(quantile(.$psi_u, probs = seq(0.85, 1, 0.05)))))
    psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
        do(data.frame(t(quantile(.$psi_u, probs = seq(0.9996, 0.9999, 0.00005)))))

    # psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
    #     do(data.frame(t(quantile(.$psi_u,
    #                              probs = seq(0.5, 0.99, length.out = 7)))))

    names(psi_quant) = c("leaf_id", paste('psi_', seq(1, 7, 1), sep = ''))

    psi_all = merge(psi_center, psi_quant, by = 'leaf_id')

    psi_long = melt(psi_all, id.vars = c("leaf_id"), value.name = "psi_star",
                    variable.name = "psi_choice")

    psi_all_df = psi_long %>%
        dplyr::mutate(logQ_cstar = 0)


    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>%
            mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                       sapply(psi_star, logJ, c_k = c_k),
                                       logQ_cstar))

    } # end of loop extracting representative points

    # for each partition (leaf_id), subset out rows for which log(Q(c)) is min
    psi_df = psi_all_df %>%
        group_by(leaf_id) %>%
        slice(which.min(logQ_cstar)) %>%  # extract rows that minimize log(Q(c))
        data.frame()

    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_df, partition, by = 'leaf_id') %>%
        dplyr::select(-psi_hat)

    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

    # this partition will be fed into any subsequent resampling functions;
    # omits extra information pertaining to loss functions, optimal psi choice
    optimal_part = partition_star %>%
        dplyr::select(-c('psi_choice', 'logQ_cstar'))

    return(list(partition_star = partition_star,
                candidate_psi  = psi_all,
                optimal_part   = optimal_part))

}
# end of u_star() function ----------------------------------------------------





#### u_star_hat() --------------------------------------------------------------
# updated 7/22 - choice of psi_star is the one that minimizes log(Q(c))
# output features: leaf_id, psi_choice (max, median, mean, hat), psi_star,
# logQ_cstar, n_obs (# of samples that partition)
#
u_star_hat = function(rpart_obj, u_df, partition, n_params) {

    u_df = u_df %>% mutate(leaf_id = rpart_obj$where)

    # extract fitted values from tree, merged with other psis later
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    n_partitions = length(partition_id)


    psi_long = melt(psi_hat_df, id.vars = c("leaf_id"), value.name = "psi_star",
                    variable.name = "psi_choice")

    psi_all_df = psi_long %>%
        dplyr::mutate(logQ_cstar = 0)


    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>%
            mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                       sapply(psi_star, logJ, c_k = c_k),
                                       logQ_cstar))

    } # end of loop extracting representative points

    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_all_df, partition, by = 'leaf_id') %>%
        dplyr::select(-psi_hat)

    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

    return(partition_star)
}





#### u_star() ------------------------------------------------------------------
# updated 4/23 - choice of psi_star is the one that minimizes log(Q(c))
# output features: leaf_id, psi_choice (max, median, mean, hat), psi_star,
# logQ_cstar, n_obs (# of samples that partition)
#
u_star_J = function(rpart_obj, u_df, partition, n_params) {


    u_df = u_df %>% mutate(leaf_id = rpart_obj$where)

    # extract fitted values from tree, merged with other psis later
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")

    n_partitions = length(partition_id)

    # compute max for each of the partitions
    # psi_all = u_df %>% dplyr::group_by(leaf_id) %>%
    #     summarise(psi_max  = max(psi_u),
    #               psi_med  = median(psi_u),
    #               psi_mean = mean(psi_u),
    #               psi_85   = quantile(psi_u, 0.85),
    #               psi_90   = quantile(psi_u, 0.90),
    #               psi_95   = quantile(psi_u, 0.95)) %>%
    #     merge(psi_hat_df, by = 'leaf_id')

    ### start updated code 4/27
    psi_center = u_df %>% dplyr::group_by(leaf_id) %>%
        summarise(psi_med  = median(psi_u),
                  psi_mean = mean(psi_u)) %>%
        merge(psi_hat_df, by = 'leaf_id')

    psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
        do(data.frame(t(quantile(.$psi_u, probs = seq(0.01, 0.15, 0.01)))))

    names(psi_quant) = c("leaf_id", paste('psi_', seq(1, 15, 1), sep = ''))

    psi_all = merge(psi_center, psi_quant, by = 'leaf_id')
    ### end updated code 4/27


    psi_long = melt(psi_all, id.vars = c("leaf_id"), value.name = "psi_star",
                    variable.name = "psi_choice")

    psi_all_df = psi_long %>%
        dplyr::mutate(logQ_cstar = 0)


    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>%
            mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                       sapply(psi_star, logJ, c_k = c_k),
                                       logQ_cstar))

    } # end of loop extracting representative points

    # for each partition (leaf_id), subset out rows for which log(Q(c)) is min
    psi_df = psi_all_df %>%
        group_by(leaf_id) %>%
        slice(which.min(logQ_cstar)) %>%  # extract rows that minimize log(Q(c))
        data.frame()

    # psi_df = psi_all_df %>%
    #     group_by(leaf_id) %>%
    #     dplyr::filter(psi_choice == 'psi_100') %>%  # extract rows that minimize log(Q(c))
    #     data.frame()


    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_df, partition, by = 'leaf_id') %>%
        dplyr::select(-psi_hat)

    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

    return(partition_star)

}
# end of u_star() function -----------------------------------------------------




#### updatePartition() ---------------------------------------------------------
# store the values from (1.2) into correct column of the final df
# partition: data frame that stores the current decisions/boundaries
# row_id   : the row corresponding to the decision to be updated
# col_id   : the TWO columns that need to be updated (parameter id)
# bdry     : the 2-d vector that has the numerical decision boundaries
updatePartition = function(partition, row_id, col_id, bdry) {

    partition[row_id, col_id:(col_id + 1)] = bdry

    return(partition)
}

# end updatePartition() function -----------------------------------------------





#### node_partition() ----------------------------------------------------------
# str_in needs only have ONE of the interval identifiers, i.e., str_split()
# already needs to have been called, splitting on "\\&"
node_partition = function(str_in, param_support) {

    LOWER = ">="
    UPPER = "<"

    # regex to subset out the numbers in the decision output
    INTERVAL = "(?>-)*[[:digit:]]+\\.*[[:digit:]]*"

    # subset out the numbers that appear in the decision string output
    decision = as.numeric(unlist(regmatches(str_in, gregexpr(INTERVAL, str_in,
                                                             perl = TRUE))))


    if (grepl(LOWER, str_in)) {          # decision: [x, ]

        # decision = as.numeric(str_split(str_in, LOWER, simplify = TRUE))
        # decision = as.numeric(unlist(regmatches(str_in,
        #                                        gregexpr(INTERVAL, str_in,
        #                                                 perl = TRUE))))

        param_id = decision[1]

        # TODO: upper bound should be the upper bound of the support
        interval = c(decision[2], param_support[param_id, 2])   # lb interval

    } else if (grepl(UPPER, str_in)) {   # decision: [ ,y]

        # decision = as.numeric(str_split(str_in, UPPER, simplify = TRUE))
        # decision = as.numeric(unlist(regmatches(str_in,
        #                                        gregexpr(INTERVAL, str_in,
        #                                                 perl = TRUE))))

        param_id = decision[1]

        # TODO: lower bound should be the lower bound of the support
        interval = c(param_support[param_id, 1], decision[2])  # ub interval

    } else {                             # decision: [x,y]

        # decision = as.numeric(unlist(regmatches(str_in,
        #                                        gregexpr(INTERVAL, str_in,
        #                                                 perl = TRUE))))

        param_id = decision[1]
        u_lb = decision[2]
        u_ub = decision[3]

        interval = c(u_lb, u_ub)
    }

    return(list(id = param_id, interval = interval))
}

# end node_partition() function ------------------------------------------------


