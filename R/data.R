###############################################################
###############################################################
######### Functions for generating hierarchical data ##########
###############################################################
###############################################################
# For generating labels and scores for a single node
GenerateNode <- function(parent_labels, cond_prob, 
                         null_shape1, null_shape2,
                         alt_shape1, alt_shape2,
                         test_n){
    # parent_labels: binary vector indicating if parent falls in alt class.
    #                same length as desired number of obs
    # cond_prob: P(node = 1 | parent = 1)
    # null/alt shapes: beta distribution parameters
    # test_n The last test_n entries of the vector are for the test data set
    # min_cases: for the node, the minimum number of observations of each
    #            type (null or alt). For example, if min_cases = 2, need
    #            at least 2 positive and 2 negative cases
    
    stopifnot(length(cond_prob) == 1)
    
    n <- length(parent_labels)
    train_n <- n - test_n
    shape1 <- c(null_shape1, alt_shape1)
    shape2 <- c(null_shape2, alt_shape2)
    probs <- c(0, cond_prob)
    
    labels <- rep(0, n)
    labels <- rbinom(n, 1, prob = probs[parent_labels + 1])
    
    #   # generate labels, ensuring that we have a
    #   # minimum number of cases of each type
    #   # to use this: add min_cases as an argument to the function    
    #   null_obs <- sum(labels[1:train_n])
    #   while(min(null_obs, n - null_obs) < min_cases){
    #     labels <- rbinom(n, 1, prob = probs[parent_labels + 1])
    #     null_obs <- sum(labels[1:train_n])
    #   }
    
    # generate scores based on these labels:
    scores <- rbeta(n, shape1 = shape1[labels + 1],
                    shape2 = shape2[labels + 1])
    
    return(list(labels = labels, scores = scores))
}

# For generating data from a tree or DAG structure
GenerateDAG <- function(parent_ind, node_info, 
                        n, test_n,
                        min_cases = 5, seed = NULL){
    # parent_ind A vector if the structure is a tree, and a list if
    #            it is a DAG. Each element contains the indices of 
    #            that node's parent nodes
    # node_info A named matrix with these six named columns:
    #           cond_prob, and the beta parameters--
    #           null_shape1, null_shape2, alt_shape1, alt_shape2
    # test_n The number of instances of the tree desired for testing. The last
    #        test_n columns of the returned matrices correspond to the test set
    # min_cases The minimum number of positives for each node in training
    # seed Random seed to use when generating data, for reproducibility
    #
    # Return two nodes x n matrices, for labels and scores, respectively.
    
    # Define some useful functions for creating nodes, findChildren and ParentsOf 
    findChildren <- function(finished, parent_ind){
        children <- which(parent_ind %in% finished)
        setdiff(children, finished)}
    if(is.list(parent_ind)) findChildren <- function(finished, parent_ind){
        # DAG: find children which have data for all of their parents:
        children <- which(sapply(parent_ind, function(x) all(x %in% finished)))
        setdiff(children, finished)}
    
    ParentsOf <- function(index, parent_ind) parent_ind[index]
    if(is.list(parent_ind)) ParentsOf <- function(index, parent_ind){
        parent_ind[[index]]
    }
    
    # Set random seed for reproducibility, if given
    if(!is.null(seed)) set.seed(seed)
    
    # Impose constraint on the minimum number of null or alt cases
    # for each node. If any node falls below this minimum, regenerate scores
    # While this requirement is not met, generate trees
    min_cases_observed <- 0 # Initialize while loop trivially
    while(min_cases_observed < min_cases){
        
        # Preallocate matrices to store labels and scores
        labels <- matrix(0, nrow = nrow(node_info), ncol = n)
        scores <- matrix(NA, nrow = nrow(node_info), ncol = n)
        parent_labels <- matrix(0, nrow = nrow(node_info), ncol = n)
        
        # Initialize tree generation loop:
        # current_nodes begins with all level 0 (root) nodes,
        # and the while loop continues to generate data until the last level
        current_nodes <- which(parent_ind == 0)
        parent_labels[current_nodes,] <- rep(1, n) # trivially true
        
        # finished is as long as number of nodes, TRUE if node's data has been created
        finished <- rep(FALSE, nrow(node_info))
        
        while(length(current_nodes) > 0){
            for(i in current_nodes){
                node_vals <- GenerateNode(parent_labels = parent_labels[i,], 
                                          cond_prob = node_info[i, "cond_prob"],
                                          null_shape1 = node_info[i, "null_shape1"],
                                          null_shape2 = node_info[i, "null_shape2"],
                                          alt_shape1 = node_info[i, "alt_shape1"],
                                          alt_shape2 = node_info[i, "alt_shape2"],
                                          test_n = test_n)
                labels[i,] <- node_vals$labels
                scores[i,] <- node_vals$scores
                finished[i] <- TRUE
            }
            
            # find child nodes, if any:
            children <- findChildren(which(finished), parent_ind)
            if(length(children) > 0){
                for(j in children){
                    # set parent label to 0 if not all parents are in the alt class
                    parent_labels[j,] <- apply(labels[ParentsOf(j, parent_ind),,drop=FALSE], 2, min)
                }
            }
            
            # continue the while loop:
            current_nodes <- children
        }
        
        pos_cases_per_node <- apply(labels[,1:n], 1, sum)
        # For checking if either pos or neg case are below threshold:
        min_cases_per_node <- pmin(pos_cases_per_node, n - pos_cases_per_node)
        
        min_cases_observed <- min(min_cases_per_node)
        min_node <- which.min(min_cases_per_node)
        
        if(min_cases_observed < min_cases) print(sprintf("Node %s less than %s cases",
                                                         min_node, min_cases))
    }
    
    return(list(labels = labels, scores = scores))
}

# For generating data from a non-hierarchical setting (just multilabel)
CreateMultiData <- function(n = 60000, test_n = 10000,
                            num_nodes = 3,
                            null_params, alt_params,
                            seed, outpath){
    # null_params and alt_params should be num_nodes x 2 matrices where
    # each row gives the neg and pos class beta parameters, respectively.
    set.seed(seed)
    prevalences <- rep(0, num_nodes)
    while(any(prevalences < 0.03)){ 
        # min prevalence set to be 0.03
        prevalences <- round(runif(num_nodes), 2)
    }
    
    # generate labels and scores, ensuring that test set actually has pos cases
    labels <- matrix(0, nrow = num_nodes, ncol = n)
    while(any(rowSums(labels[, (n-test_n+1):n]) <= 3)){
        # min positive cases is set to 3 for the test set
        labels <- t(sapply(prevalences, function(x) rbinom(n, 1, x)))
    }
    scores <- matrix(NA, nrow = num_nodes, ncol = n)
    for(i in 1:num_nodes){
        is_pos <- labels[i,] > 0
        scores[i, is_pos] <- rbeta(sum(is_pos),
                                   alt_params[i, 1], 
                                   alt_params[i, 2])
        scores[i, !is_pos] <- rbeta(sum(!is_pos), 
                                    null_params[i, 1],
                                    null_params[i, 2])
    }
    
    node_info <- cbind(parent_ind = rep(0, num_nodes),
                       cond_prob = prevalences,
                       null_shape1 = null_params[,1],
                       null_shape2 = null_params[,2],
                       alt_shape1 = alt_params[,1],
                       alt_shape2 = alt_params[,2])
    
    saveRDS(list(labels = labels, 
                 scores = scores,
                 seed = seed,
                 node_info = node_info),
            file = outpath)
}

# These are different functions for generating the data in batch
# for each simulation setting.
# They are written so you don't need to feed it any input except for
# the number of data sets you'd like to generate (the number of 
# times you'd like to repeat the simulation).
# The function is written to be reproducible, so the same random seeds
# are used to generate the data each time. 
#' Generate synthetic data based on Beta distributions.
#'
#' @param n Number of total objects, including training and testing.
#' @param test_n Number of testing objects
#' @param min_cases Minimal number of positive cases in each class.
#' @param unif_cond A logic value, whether using purely uniform conditional probablity.
#'                  A FALSE value accounts for the joint probablity to take into account the \code{min_cases}.
#' @param setting The simulation setting. There are 11 hard-coded setting.
#' @param seed Random seed.
#' 
#' @return A list with four elements.
#' \describe{
#'  \item{labels}{A K x n matrix containing labels for n objects across K classes.}
#'  \item{scores}{A K x n matrix containing simulated classifier scoresfor n objects across K classes.}
#'  \item{seed}{Random seed.}
#'  \item{node_info}{A K x 6 matrix containing the information of K classes: the first columns is the
#'                   the parent of the node, the second column is the conditional probablity \code{Pr(Y_i|Y_{pa(i)})},
#'                   the third and the columns are two parameters of the null Beta distribution, the fifth and the
#'                   sixth columns are two parameters of the alternative Beta distribution.}
#'  }
#'  
#' @export
#' 
CreateDataForSetting <- function(n = 60000, test_n = 10000, 
                                 min_cases = 150, unif_cond = TRUE,
                                 setting, seed){
    # Use the same random seed for each of the following steps:
    # 1. generating prevalences, 
    # 2. generating the data for the tree, given the prevalences
    # 3. generating the bags to use for bagged ClusHMC
    set.seed(seed)
    
    # Define a matrix with node parameters:
    if(setting == 1){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(NA, NA, NA),
                           null_shape1 = c(2, 5.5, 2),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(6, 6, 6),
                           alt_shape2 = c(2, 5.5, 2))
    }
    if(setting == 2){
        node_info <- cbind(parent_ind = c(0, 1, 1),
                           cond_prob = rep(NA, 3),
                           null_shape1 = c(2, 5.5, 2),
                           null_shape2 = rep(6, 3),
                           alt_shape1 = rep(6, 3),
                           alt_shape2 = c(2, 5.5, 2))
    }
    if(setting == 3){
        node_info <- cbind(parent_ind = c(0, 0, 0, 2, 3, 3, 4, 4, 
                                          4, 6, 6, 7, 7, 8, 8, 9, 
                                          9, 10, 10, 11, 12, 13, 
                                          13, 13, 14),
                           cond_prob = c(0.8, 0.8, 0.8, 0.3, 0.8, 0.1,
                                         0.5, 0.8, 0.1, 0.3, 0.5, 0.5, 
                                         0.8, 0.3, 0.1, 0.5, 0.8, 0.5, 
                                         0.8, 0.5, 0.05, 0.8, 0.3, 0.3, 0.1),
                           null_shape1 = c(5.5, 5.5, 5.5, 5.5, 4, 4, 2, 
                                           2, 2, 4, 2, 5.5, 5.5, 2, 5.5, 
                                           4, 5.5, 4, 4, 2, 5.5, 4, 4, 
                                           5.5, 4),
                           null_shape2 = rep(6, 25),
                           alt_shape1 = rep(6, 25),
                           alt_shape2 = c(5.5, 5.5, 5.5, 5.5, 4, 4, 2, 
                                          2, 2, 4, 2, 5.5, 5.5, 2, 5.5, 
                                          4, 5.5, 4, 4, 2, 5.5, 4, 4, 
                                          5.5, 4))
        
    }
    if(setting == 4){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(NA, NA, NA),
                           null_shape1 = c(2, 2, 5.5),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(6, 6, 6),
                           alt_shape2 = c(2, 2, 5.5))
    }
    if(setting == 5){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(NA, NA, NA),
                           null_shape1 = c(5.5, 5.5, 2),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(6, 6, 6),
                           alt_shape2 = c(5.5, 5.5, 2))
    }
    if(setting == 6){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(NA, NA, NA),
                           null_shape1 = c(2, 2, 2),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(2, 2, 2),
                           alt_shape2 = c(6, 6, 6))
    }
    if(setting == 7){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(NA, NA, NA),
                           null_shape1 = c(2, 2, 2),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(6, 6, 6),
                           alt_shape2 = c(2, 2, 2))
    }
    if(setting == 8){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(0.95, 0.95, 0.95),
                           null_shape1 = c(2, 2, 2),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(6, 6, 6),
                           alt_shape2 = c(2, 2, 2))
    }
    if(setting == 9){
        node_info <- cbind(parent_ind = c(0, 1, 2), 
                           cond_prob = c(0.95, 0.95, 0.95),
                           null_shape1 = c(2, 5.5, 2),
                           null_shape2 = c(6, 6, 6),
                           alt_shape1 = c(6, 6, 6),
                           alt_shape2 = c(2, 5.5, 2))
    }
    if(setting == 10){
        node_info <- cbind(parent_ind = c(0, 0, 0, 2, 3, 3, 4, 4, 
                                          4, 6, 6, 7, 7, 8, 8, 9, 
                                          9, 10, 10, 11, 12, 13, 
                                          13, 13, 14),
                           cond_prob = c(0.8, 0.8, 0.8, 0.3, 0.8, 0.1,
                                         0.5, 0.8, 0.1, 0.3, 0.5, 0.5, 
                                         0.8, 0.3, 0.1, 0.5, 0.8, 0.5, 
                                         0.8, 0.5, 0.05, 0.8, 0.3, 0.3, 0.1),
                           null_shape1 = rep(2, 25),
                           null_shape2 = rep(6, 25),
                           alt_shape1 = rep(6, 25),
                           alt_shape2 = rep(2, 25))
    }
    if(setting == 11){
        node_info <- cbind(parent_ind = c(0, 0, 0, 2, 3, 3, 4, 4, 
                                          4, 6, 6, 7, 7, 8, 8, 9, 
                                          9, 10, 10, 11, 12, 13, 
                                          13, 13, 14),
                           cond_prob = sample(c(0.8, 0.8, 0.8, 0.3, 0.8, 0.1,
                                                0.5, 0.8, 0.1, 0.3, 0.5, 0.5, 
                                                0.8, 0.3, 0.1, 0.5, 0.8, 0.5, 
                                                0.8, 0.5, 0.05, 0.8, 0.3, 0.3, 0.1)),
                           null_shape1 = c(5.5, 5.5, 5.5, 5.5, 4, 4, 2, 
                                           2, 2, 4, 2, 5.5, 5.5, 2, 5.5, 
                                           4, 5.5, 4, 4, 2, 5.5, 4, 4, 
                                           5.5, 4),
                           null_shape2 = rep(3.5, 25),
                           alt_shape1 = rep(3.5, 25),
                           alt_shape2 = c(5.5, 5.5, 5.5, 5.5, 4, 4, 2, 
                                          2, 2, 4, 2, 5.5, 5.5, 2, 5.5, 
                                          4, 5.5, 4, 4, 2, 5.5, 4, 4, 
                                          5.5, 4))
        
    }
    
    
    # Generate conditional probabilities so that the probability of 
    # obtaining more than min_cases pos. observations at the deepest level node
    # is sufficiently high (so as to not be changing the true distribution
    # of the tree generation process too much when we require having
    # at least min_cases observations.)
    computeJointProb <- function(cond_prob, parent_ind){
        while(any(parent_ind != 0)){
            idx <- which(parent_ind != 0)
            cond_prob[idx] <- cond_prob[idx]*cond_prob[parent_ind[idx]]
            parent_ind[idx] <- parent_ind[parent_ind[idx]]
        }
        return(cond_prob)
    }
    
    # Generate a set of conditional probabilities, but make sure none are
    # too small so that we can't realistically achieve the minimums we want
    accept_condprob <- FALSE
    if(unif_cond){
        while(!accept_condprob){
            cond_prob <- round(runif(n = nrow(node_info), max = 0.5), 2)
            p <- prod(cond_prob)
            accept_condprob <- p > min_cases/(n - test_n)*1.25 # x1.25 for wiggle room
        }
        node_info[,2] <- cond_prob
    } else{
        while(!accept_condprob){
            #       cond_prob <- sample(c(0.8, 0.5, 0.3, 0.15, 0.05), nrow(node_info),
            #                           prob = c(1.3, 1.75, 2, 1.5, 1), replace = TRUE)
            cond_prob <- runif(nrow(node_info))
            joint_probs <- computeJointProb(cond_prob, node_info[,1])
            p <- min(c(joint_probs, 1 - joint_probs)) # chance of pos or neg case
            accept_condprob <- p > min_cases/(n - test_n)*1.25 # x1.25 for wiggle room
        } 
        node_info[,2] <- cond_prob
    }
    
    sim_data <- GenerateDAG(parent_ind = node_info[,1], 
                            node_info = node_info[,2:6], 
                            n = n, test_n = test_n,
                            min_cases = min_cases, 
                            seed = seed)
    sim_data$seed <- seed
    sim_data$node_info <- node_info
    
    return(sim_data)
}