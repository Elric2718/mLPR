
#' Sort a vector of values using naive sorting or HierRank.
#' 
#' @param scores A K x n matrix 
#' @param graph_structure A list that indicates the parents of each node.
#' @param sorting_method The sorting method, including 'value' (naive sorting in descending order), 'HierRank'.
#'                       
#' @return A list.
#' \describe{
#'  \item{ordering}{The resulting ranking.}
#'  \item{values}{The transposed \code{scores}.}
#'  \item{block_ordering}{The ordered blocks in HierRank method (only for 'HierRank').}
#'  }
#'  
#' @export
Sorting <- function(scores, graph_structure, sorting_method = "value"){
    n_test <- nrow(scores)
    
    if(sorting_method == "value"){
        return(list(ordering = order(t(scores), decreasing = TRUE), values = t(scores)))
    }else if(sorting_method == "HierRank"){
        
        # for tree
        # parent_ind <- rep(unlist(graph_structure), times = n_test)
        # parent_ind <- parent_ind + rep((0:(n_test-1)) * n_class, each = n_class)
        # parent_ind[ parent_ind %% n_class == 0 ] <- 0
        
        # for DAG
        #parent_id <- replicateList(graph_structure, n_test)
        #cssa_pred <- cssa(parent_id, t(scores), parallel = FALSE)
        
        cssa_pred <- cssa(graph_structure, t(scores), parallel = TRUE)
        return(list(ordering = unlist(cssa_pred), values = t(scores), block_ordering = cssa_pred))
    }
}


#### #### #### #### #### #### #### ####
# coding the CSSA algorithm
#### #### #### #### #### #### #### ####

#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
cssa <- function(parent_ind, LPR, 
                 parallel = FALSE){
    # parallel cssa sorts each tree (per person) in parallel
    # and then collates into a final sorting vector
    # in other words, cssa is run separately for each column of LPR
    # parent_ind is a numeric vector whose length depends on whether
    #            you choose parallel or sequential running.
    #            if parallel, then parent_ind should be the same length as the
    #            number of nodes in one tree.
    #            if sequential, parent_ind is extended to be length
    #            ncol(LPR)*nrow(LPR) by repeating the original parent_ind ncol(LPR) 
    #            times and changing indices. Ex: c(0, 1, 2) for two trees becomes
    #            c(0, 1, 2, 0, 3, 4)
    # LPR is a matrix with same number of rows as nodes in the tree
    
    n <- length(parent_ind)
    tree <- ifelse(is.list(parent_ind), FALSE, TRUE)
    
    if(parallel){
        # parent_ind is a vector of length n
        # expect LPR to be an n x k matrix of LPR scores for k samples
        stopifnot("matrix" %in% class(LPR))
        stopifnot(nrow(LPR) == n)
        
        k <- ncol(LPR)
        
        if(tree){
            final_list <- foreach::foreach(i=1:k, .combine = c) %dopar% {
                i1 <- (i-1)*n + 1
                i2 <- i*n
                cssa_tree(parent_ind, LPR[,i], ind_names = i1:i2)
            }
        } else{
            final_list <- foreach::foreach(i=1:k, .combine = c) %dopar% {
                i1 <- (i-1)*n + 1
                i2 <- i*n
                cssa_dag_and(parent_ind, LPR[,i], ind_names = i1:i2)
            }
        }
        
        # sort the results in the list by avg LPR
        final_avgs <- sapply(final_list, function(x) mean(LPR[x]))
        final_list <- final_list[order(final_avgs, decreasing = TRUE)]
        
        #return(unlist(final_list))
        return(final_list)
    } else{
        stopifnot(length(LPR) == length(parent_ind))
        
        if(tree){
            output <- cssa_tree(parent_ind = parent_ind, LPR = as.vector(LPR))
            output <- unname(output)
        } else{ 
            output <- cssa_dag_and(parent_ind = parent_ind, LPR = as.vector(LPR))
            output <- unname(output)
        }
        
        final_avgs <- sapply(output, function(x) mean(LPR[x]))
        final_list <- output[order(final_avgs, decreasing = TRUE)]
        #return(unlist(final_list))
        return(final_list)
        
    }
    
}


cssa_tree <- function(parent_ind, LPR_vec, 
                      ind_names = 1:length(parent_ind), 
                      verbose=FALSE){
    # LPR is a vector of LPRs for each node
    # parent_ind is a vector giving the parent of each node,
    # it is assumed that nodes are numbered 1:n, and
    # a parent of 0 = that node is at root level
    # ind_names is an optional vector with the true indices of the
    # nodes, if they are coming from a larger graph and they are not
    # actually 1:n
    # more--return a more detailed list (TRUE) or just the final sort?
    
    stopifnot(is.vector(LPR_vec))
    stopifnot(length(parent_ind) == length(LPR_vec))
    
    # initialize while loop
    n <- length(LPR_vec)
    admitted <- rep(FALSE, n)     # logical: include supernode?
    nodelist <- split(1:n, 1:n) # to contain indices for supernodes
    avgs <- LPR_vec               # to contain averages across supernodes
    parent <- parent_ind          # parent[i] is the parent of node i
    i <- which.max(avgs)
    
    while(length(i) > 0){
        if(parent[i] == 0){ # admit the node if it is already at root level
            admitted[i] <- TRUE 
        } else{
            if(admitted[parent[i]] == TRUE){ # or if its parent has been admitted
                admitted[i] <- TRUE
            } else{ # otherwise, combine the node with its parent
                nodelist[[parent[i]]] <- c(nodelist[[parent[i]]], nodelist[[i]])
                avgs[parent[i]] <- mean(LPR_vec[nodelist[[parent[i]]]])
                parent[parent==i] <- parent[i]
            }
        }
        # replace node avg with NA so that while loop will
        # continue on to the next highest avg
        avgs[i] <- NA
        i <- which.max(avgs)
    }
    
    # keep only the nodes which have been admitted,
    # and use the real names of the indices, if provided
    finallist = lapply(nodelist[admitted], function(x) ind_names[x])
    
    # the supernodes are returned in arbitrary order, but that's okay because
    # the supernodes now have values that are T-decreasing
    # so it is impossible for a child node to be ordered before its parent
    return(finallist)
}

cssa_dag_and <- function(parent_ind, LPR_vec, 
                         ind_names = 1:length(parent_ind),
                         verbose=FALSE){
    # parent_ind is a list now
    
    stopifnot(is.vector(LPR_vec))
    stopifnot(length(parent_ind) == length(LPR_vec))
    
    # initialize before while loop
    n <- length(LPR_vec)
    admitted <- rep(FALSE, n)     # logical: include supernode?
    nodelist <- split(1:n, 1:n) # to contain indices for supernodes
    avgs <- LPR_vec               # to contain averages across supernodes
    parent <- parent_ind          # parent[i] is the parent of node i
    i <- which.max(avgs)
    
    while(length(i) > 0){
        if(all(parent[[i]] == 0)){ # admit the node if it is already at root level
            admitted[i] <- TRUE 
        } else{
            if(all(admitted[parent[[i]]]) == TRUE){ # or if all parents have been admitted
                admitted[i] <- TRUE
            } else{ # otherwise, combine node and parent with smallest avg LPR
                min_parent <- parent[[i]][which.min(avgs[parent[[i]]])]
                nodelist[[min_parent]] <- c(nodelist[[min_parent]], nodelist[[i]])
                avgs[min_parent] <- mean(LPR_vec[nodelist[[min_parent]]])
                parent <- lapply(parent, function(x) ifelse(x == i, min_parent, x))
            }
        }
        # replace node avg with NA so that while loop will
        # continue on to the next highest avg
        avgs[i] <- NA
        i <- which.max(avgs)
    }
    
    # keep only the nodes which have been admitted,
    # and use the real names of the indices, if provided
    finallist = lapply(nodelist[admitted], function(x) ind_names[x])
    return(finallist)  
}


