#############################################################################################
########################################## mLPR ##########################################
#############################################################################################
# require(dplyr)
# require(LPRclass)
# require(foreach)
# require(doParallel)
# require(gRain)
# require(locpol)


eps <- 1e-6
lower_clip = 0
upper_clip = 1


#' Estimate LPR
#'
#' This function estimates the mLPR on a data set.
#' The user needs to provide two data sets (one to use as training
#' and the other as a test set). 
#'
#' @param S_train Numeric matrix containing classifier scores for each training observation, 
#'                of shape=n_train x K, where n_train is the number of training objects and 
#'                K is the number of class. 
#' @param S_test  Numeric matrix containing classifier scores for each testing observation, 
#'                of shape=n_test x K, where n_test is the number of testing objects and 
#'                K is the number of class. 
#' @param Y_train Binary matrix containing labels for each training observation, of the shape
#'                as \code{S_train}.
#' @param graph_structure A list that indicates the parents of each node.
#' @param X_train Numeric matrix containing the features for training data. 
#' @param X_test Numeric matrix containing the features for testing data.
#' @param cdensity_method Methods of estimating the conditional probablity of \code{Pr(Y_i|Y_pa(i))},
#'                        including 'empirical', 'logistic' (requires \code{X_train} and \cdoe{X_test}), 
#'                        'svm' (requires \code{X_train} and \cdoe{X_test}).
#' @param pr_s_y_module A string indicating which method and model to estimate \code{Pr(S|Y) or Pr(Y|S)}. The method
#'                     and the model is concatenated by '-'. Methods include 'LPR', 'Gaussian', 'thresholding',
#'                     models include 'logistic', 'svm'.
#' @param pr_y_module A string indicating which method and model to estimate \code{Pr(Y)}. The method
#'                    and the model is concatenated by '-'. Methods include 'hierarchical', 'empirical',
#'                    models include 'logistic', 'svm'.
#' @param approx_mlpr Three versions of mLPR estimations: 'full', 'nbh', 'indpt'.
#' @param num_cores Number of cores for parallel computing.
#' @param saving_name The prefix name used to save intermediate objects.
#' 
#' @return A numetric matrix containing mLPR values for testing objects, of shape=n_test x K.
#' 
#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
#'
#' @export
#'
estimatemLPR <- function(S_train, 
                         S_test, 
                         Y_train,
                         graph_structure,
                         X_train = NULL,
                         X_test = NULL,
                         cdensity_method = "empirical",
                         pr_s_y_module = "LPR-",
                         pr_y_module= "hierarchical-", 
                         approx_mlpr = "full",
                         num_cores=1,
                         saving_name = "",
                         ...){
    doParallel::registerDoParallel(cores=num_cores)
    params <- list(...)
    n_train <- nrow(S_train)
    n_test <- nrow(S_test)
    n_class <- ncol(Y_train)
    
    pr_s_y_module <- strsplit(pr_s_y_module, '-')[[1]]
    pr_s_y_method <- pr_s_y_module[1]
    pr_s_y_fit_model <- pr_s_y_module[2]
    
    pr_y_module <- strsplit(pr_y_module, '-')[[1]]
    pr_y_method <- pr_y_module[1]
    pr_y_fit_model <- pr_y_module[2]
    
    
    if(saving_name == ""){
        saving_name = "~/mLPR_tmp/"
        unlink(saving_name, recursive=TRUE)
        dir.create(saving_name)
    }
    
    ###########################################
    # Step I: Compute the potential functions #
    ###########################################
    print("Step I: Compute the potential functions")
    print("(i) Compute p(y_i | y_{pa(i)})")
    ##### (i) Compute p(y_i | y_{pa(i)})
    if(!is.null(saving_name)){
        cdensity_path <- paste0(saving_name, "_cdensity_", ifelse(is.null(X_train), "null_", "cov_"), cdensity_method, ".RData")
    }else{
        cdensity_path <- "NULL"
    }
    if(!file.exists(cdensity_path)){
        cdensity <- CondDensity(Y_train, X_train, graph_structure, method = cdensity_method)
        # a function that takes in nodes (character) to investigate, Y values and X values.
        if(cdensity_path != "NULL"){
            save(cdensity, file = cdensity_path)
        }
    }else{
        print(paste0("Existed: ", cdensity_path, ". Loading."))
        load(cdensity_path)
    }
    
    
    
    # cdensity_test <- foreach(i = seq(n_test)) %dopar%{
    #   lapply(seq(n_class), function(node_id){
    #     if(is.matrix.csr(X_test)){
    #       cdensity(node_id, X_test[i,])
    #     }else{
    #       cdensity(node_id, X_test[i,,drop = FALSE])
    #     }
    #   })
    # }
    if(!is.null(saving_name)){
        cdensity_test_path <- paste0(saving_name, "_cdensity_test_", ifelse(is.null(X_train), "null_", "cov_"), cdensity_method, ".RData")
    }else{
        cdensity_test_path <- "NULL"
    }
    if(!file.exists(cdensity_test_path)){
        cdensity_test <- lapply(seq(n_class), function(node_id){
            cdensity(node_id, n_test, X_test)
        })
        if(cdensity_test_path != "NULL"){
            save(cdensity_test, file = cdensity_test_path)
        }
    }else{
        print(paste0("Existed: ", cdensity_test_path, ". Loading."))
        load(cdensity_test_path)
    }
    
    
    print("(ii) Compute p(y_i)")
    ##### (ii) Compute p(y_i)
    if(!is.null(saving_name)){
        pr_y_path <- paste0(saving_name, "_pr_y_", ifelse(is.null(X_train), "null_", "cov_"), cdensity_method, "_", pr_y_method, "_", pr_y_fit_model, ".RData")
    }else{
        pr_y_path <- "NULL"
    }
    if(!file.exists(pr_y_path)){
        pr_y <- LearnPr_y(graph_structure, n_test, cdensity_test, Y_train, X_train, method = pr_y_method, fit_model = pr_y_fit_model)
        # a function that takes nodes (character) to investigate and X values.
        if(pr_y_path != "NULL"){
            save(pr_y, file = pr_y_path)
        }
    }else{
        print(paste0("Existed: ", pr_y_path, ". Loading."))
        load(pr_y_path)
    }
    
    
    #pr_y_train <- pr_y(NULL, X_train)
    if(!is.null(saving_name)){
        pr_y_test_path <- paste0(saving_name, "_pr_y_test_", ifelse(is.null(X_train), "null_", "cov_"), cdensity_method, "_", pr_y_method, "_", pr_y_fit_model, ".RData")
    }else{
        pr_y_test_path <- "NULL"
    }
    if(!file.exists(pr_y_test_path)){
        pr_y_test <- pr_y(NULL, X_test) %>% matrix(nrow = n_test, ncol = n_class)
        if(pr_y_test_path != "NULL"){
            save(pr_y_test, file = pr_y_test_path)
        }
    }else{
        print(paste0("Existed: ", pr_y_test_path, ". Loading."))
        load(pr_y_test_path)
    }
    
    print("(iii) Compute p(s_i | y_i)")
    ##### (iii) Compute p(s_i | y_i)
    if(!is.null(saving_name)){
        pr_s_y_path <- paste0(saving_name, "_pr_s_y_", ifelse(is.null(X_train), "null_", "cov_"), pr_s_y_method, "_", pr_s_y_fit_model, ".RData")
    }else{
        pr_s_y_path <- "NULL"
    }
    if(!file.exists(pr_s_y_path)){
        if(pr_s_y_method == "thresholding"){
            
            threshold <- params$threshold_s_y
            Y_hat_train <- Threshold(S_train, threshold)
            #Y_hat_test <- Threshold(S_test, threshold)
            
            pr_s_y <- LearnPr_S_Y(Y_hat_train, Y_train, X_train, method = "binomial", fit_model = pr_s_y_fit_model) 
            
        }else if(pr_s_y_method == "Gaussian"){
            pr_s_y <- LearnPr_S_Y(S_train, Y_train, X_train, method = "gaussian", fit_model = pr_s_y_fit_model)
        }else if(pr_s_y_method == "LPR"){
            # Do not get a separate classifier.
        }
        # pr_s_y_train and pr_s_y_test are two functions that take in nodes (character) to investigate, Y values and X values.
        if(pr_s_y_method != "LPR"){
            if(pr_s_y_path != "NULL"){
                save(pr_s_y, file = pr_s_y_path)  
            }
        }
    }else{
        print(paste0("Existed: ", pr_s_y_path, ". Loading."))
        load(pr_s_y_path)
    }
    
    if(!is.null(saving_name)){
        lpr_test_path <- paste0(saving_name, "_lpr_test_", ifelse(is.null(X_train), "null_", "cov_"), pr_s_y_method, "_", pr_s_y_fit_model, ".RData")
    }else{
        lpr_test_path <- "NULL"
    }
    if(!file.exists(lpr_test_path)){
        if(pr_s_y_method == "thresholding"){
            threshold <- params$threshold_s_y
            #Y_hat_train <- Threshold(S_train, threshold)
            Y_hat_test <- Threshold(S_test, threshold) %>% matrix(nrow = n_test, ncol = n_class)
            
            #pr_s_y_train <- learner_s_y(NULL, Y_hat_train, X_train)
            pr_s_y_test <- pr_s_y(NULL, Y_hat_test, X_test) %>% matrix(nrow = n_test, ncol = n_class)
        }else if(pr_s_y_method == "Gaussian"){
            #pr_s_y_train <- pr_s_y(NULL, S_train, X_train)
            pr_s_y_test <- pr_s_y(NULL, S_test, X_test) %>% matrix(nrow = n_test, ncol = n_class)
        }else if(pr_s_y_method == "LPR"){
            try_rlt <- try({
                lpr_test <- foreach::foreach(j=seq(n_class), .combine = cbind) %dopar% {
                    estimateLPR(scores = S_train[,j], 
                                labels = Y_train[,j], 
                                newdata = S_test[, j], 
                                num_bags = 0,
                                method = "spline")
                } %>% matrix(nrow = n_test, ncol = n_class)
            }, silent = TRUE)
            if(class(try_rlt) == "try-error"){
                lpr_test <- foreach::foreach(j=seq(n_class), .combine = cbind) %dopar% {
                    estimateLPR(scores = S_train[,j], 
                                labels = Y_train[,j], 
                                newdata = S_test[, j], 
                                num_bags = 0,
                                method = "kernel")
                } %>% matrix(nrow = n_test, ncol = n_class)
            }
        }
        
        if(pr_s_y_method != "LPR"){
            # lpr_train <- sapply(seq(n_class), function(j){
            #     pr_joint_s_y_train <- pr_s_y_train[[j]] * cbind(pr_y_train[, j], 1 - pr_y_train[, j])
            #     pr_joint_s_y_train[, 1]/rowMeans(pr_joint_s_y_train)
            # }) 
            lpr_test <- sapply(seq(n_class), function(j){
                pr_joint_s_y_test <- pr_s_y_test[[j]] * cbind(pr_y_test[, j], 1 - pr_y_test[, j])
                pr_joint_s_y_test[, 1]/(rowSums(pr_joint_s_y_test) + eps)
            }) %>% matrix(nrow = n_test, ncol = n_class)
        }
        if(lpr_test_path != "NULL"){
            save(lpr_test, file = lpr_test_path)
        }
    }else{
        print(paste0("Existed: ", lpr_test_path, ". Loading."))
        load(lpr_test_path)
    }
    
    ##########################################
    # Step II: Compute p(y_i| s_1, ..., s_N) #
    ##########################################
    print("Step II: Compute p(y_i| s_1, ..., s_N)")    
    if(approx_mlpr == "indpt"){
        #mlpr_train <- lpr_train
        mlpr_test <- lpr_test
    }else if(approx_mlpr == "neighborhood"){
        yn <- c("yes", "no")
        ch_set <- lapply(seq(n_class), function(node_id){
            #pattern <- paste0("^", node_id, ";|;", node_id, ";|;", node_id, "$")
            sapply(seq(n_class), function(j){
                if(node_id %in% graph_structure[[j]]){
                    return(j)
                }else{
                    return(NULL)
                }
            }) %>% unlist()
        })
        mlpr_test <- foreach::foreach(i=seq(n_test), .combine = rbind) %dopar% {
            sapply(seq(n_class), function(node_id){
                #print(node_id)
                pa_id <- graph_structure[[node_id]]#as.numeric(strsplit(graph_structure[node_id], split = ";")[[1]])
                ch_id <- ch_set[[node_id]]
                
                ### (a) prepare the potential functions
                # pr(j); j in pa(i)
                if(!all(pa_id == 0)){
                    pr_pa <- lapply(pa_id, function(j){
                        #pr_pa <- lapply(seq(pa_id), function(j){
                        vpar <- paste0("~ Node", j)
                        values <- c(lpr_test[i, j], 1 - lpr_test[i, j])
                        
                        eval(parse(text = paste0("gRain::cptable(", vpar,
                                                 ", values = c(", paste0(values, collapse = ","),
                                                 "), levels = yn, normalize = FALSE, smooth = eps)")))
                    })
                }else{
                    pr_pa <- NULL
                }
                
                # pr(i | pa(i))
                if(all(pa_id == 0)){
                    vpar <- paste0("~ Node", node_id)
                    values <- c(lpr_test[i, node_id], 1 - lpr_test[i, node_id])
                }else{
                    vpar <- paste0("~ Node", node_id, " | ", paste0("Node", pa_id, collapse = " + "))
                    values <-  c(rbind(cdensity_test[[node_id]][["joint"]][i, ], 1 - cdensity_test[[node_id]][["joint"]][i, ]) *
                                     c(lpr_test[i, node_id], 1 - lpr_test[i, node_id]) *
                                     1/c(Clip(pr_y_test[i, node_id] + eps, lower_clip, upper_clip), Clip(1 - pr_y_test[i, node_id] + eps, lower_clip, upper_clip)))
                    # values <-  c(rbind(cdensity_test[[node_id]][["joint"]][i, ], 1 - cdensity_test[[node_id]][["joint"]][i, ]) * 
                    #                c(lpr_test[i, node_id], 1 - lpr_test[i, node_id]))
                    #c(lpr_test[i, node_id], 1 - lpr_test[i, node_id]) * 
                    #1/c(pr_y_test[i, node_id], 1 - pr_y_test[i, node_id])
                }
                eval(parse(text = paste0("pr_i_pa <- gRain::cptable(", vpar,
                                         ", values = c(", paste0(values, collapse = ","),
                                         "), levels = yn, normalize = FALSE, smooth = eps)")))
                
                # pr(j | i); j in ch(i)
                if(length(ch_id) > 0){
                    pr_ch_i <- lapply(ch_id, function(j){
                        try_rlt <- try({
                            loc <- which(node_id == graph_structure[[j]])
                            
                            vpar <- paste0("~ Node", j, " | Node", node_id)
                            values <- c(lpr_test[i, j], 1 - lpr_test[i, j]) *
                                1/c(Clip(pr_y_test[i, j] + eps, lower_clip, upper_clip), Clip(1 - pr_y_test[i, j] + eps, lower_clip, upper_clip)) *
                                c(rbind(cdensity_test[[j]][["single"]][i, c(2 * loc - 1, 2 * loc)],
                                        1 - cdensity_test[[j]][["single"]][i, c(2 * loc - 1, 2 * loc)]))
                            # values <- c(lpr_test[i, j], 1 - lpr_test[i, j]) *
                            #   #1/c(pr_y_test[i, j] + eps, 1 - pr_y_test[i, j] + eps) *
                            #   c(rbind(cdensity_test[[j]][["single"]][i, c(2 * loc - 1, 2 * loc)],
                            #           1 - cdensity_test[[j]][["single"]][i, c(2 * loc - 1, 2 * loc)]))
                            
                            eval(parse(text = paste0("gRain::cptable(", vpar,
                                                     ", values = c(", paste0(values, collapse = ","),
                                                     "), levels = yn, normalize = FALSE, smooth = eps)")))
                        })
                        if(class(try_rlt) == 'try-error'){
                            browser()
                        }else{
                            try_rlt
                        }
                    })
                }else{
                    pr_ch_i <- NULL
                }
                # (b) run the junction tree algorithm                
                plist <- gRain::compileCPT(c(pr_pa, list(pr_i_pa), pr_ch_i))
                gin <- gRain::grain(plist)
                
                return(gRain::querygrain(gin, nodes = paste0("Node", node_id), type = "marginal")[[1]][1])
            })
        }
        
    }else if(approx_mlpr == "full"){
        yn <- c("yes", "no")
        mlpr_test <- foreach::foreach(i=seq(n_test), .combine = rbind) %dopar% {
            # (a) prepare the potential functions
            potential_list <- lapply(seq(n_class), function(node_id){
                pa_id <- graph_structure[[node_id]] #as.numeric(strsplit(graph_structure[node_id], split = ";")[[1]])
                if(all(pa_id == 0)){
                    vpar <- paste0("~ Node", node_id)
                    values <- c(lpr_test[i, node_id], 1 - lpr_test[i, node_id])
                }else{
                    vpar <- paste0("~ Node", node_id, " | ", paste0("Node", pa_id, collapse = " + "))
                    values <- c(rbind(cdensity_test[[node_id]][["joint"]][i, ], 1 - cdensity_test[[node_id]][["joint"]][i, ]) *
                                    c(lpr_test[i, node_id], 1 - lpr_test[i, node_id]) *
                                    1/c(Clip(pr_y_test[i, node_id] + eps, lower_clip, upper_clip), Clip(1 - pr_y_test[i, node_id] + eps, lower_clip, upper_clip)))
                    # values <- rbind(rep(lpr_test[i, node_id], length(cdensity_test[[node_id]][["joint"]][i, ])), 
                    #                 1 - rep(lpr_test[i, node_id], length(cdensity_test[[node_id]][["joint"]][i, ])))
                    # values <- c(rbind(cdensity_test[[node_id]][["joint"]][i, ], 1 - cdensity_test[[node_id]][["joint"]][i, ]) *
                    #               c(lpr_test[i, node_id], 1 - lpr_test[i, node_id]))
                    #cdensity_test[[i]][[node_id]][["joint"]] * 
                    #c(rbind(lpr_test[i, node_id], 1 - lpr_test[i, node_id])) * 
                    #1/c(rbind(pr_y_test[i, node_id], 1 - pr_y_test[i, node_id]))
                }
                eval(parse(text = paste0("tmp <- gRain::cptable(", vpar,
                                         ", values = c(", paste0(values, collapse = ","),
                                         "), levels = yn, normalize = FALSE, smooth = eps)")))
                return(tmp)
            })
            # (b) run the junction tree algorithm
            plist <- gRain::compileCPT(potential_list)
            gin <- gRain::grain(plist)
            
            mlpr_list <- gRain::querygrain(gin, nodes = paste0("Node", seq(n_class)), type = "marginal")
            return(sapply(seq(n_class), function(node_id) mlpr_list[[paste0("Node", node_id)]][[1]][1]))
        }
    }
    
    return(mlpr_test)
}

#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
CondDensity <- function(Y_train, X_train, graph_structure, method){
    n_nodes <- length(graph_structure)
    pop_marginal_pi <- colMeans(Y_train)
    
    cond_rlt <- foreach::foreach(node_id = seq(n_nodes)) %dopar% {
        
        pa_id <- graph_structure[[node_id]] #as.numeric(strsplit(graph_structure[node_id], split = ";")[[1]])
        
        if(all(pa_id == 0)){
            if(is.null(X_train) | method == "empirical"){
                return(pop_marginal_pi[node_id])
            }else if(all(Y_train[, node_id] == 1) | all(Y_train[, node_id] == 0)){
                return(ExtremeClassifier(X_train, Y_train[, node_id]))
            }else if(method == "logistic"){
                return(glm(y~., family = binomial, data = data.frame(y = Y_train[, node_id], x = X_train)))
            }else if(method == "svm"){
                return(LinearSVM(X_train, Y_train[, node_id]))
            }
        }else{
            n_pa <- length(pa_id)
            # pr(i| pa(i))
            tmp_rlt_joint <- lapply(seq(2^(n_pa), 1), function(j){
                y_val <- as.integer(intToBits(j - 1))[seq(n_pa)]
                valid_sums <- sapply(seq(n_pa), function(tmp_id){
                    Y_train[, pa_id[tmp_id]] == y_val[tmp_id]
                }) %>% rowSums()
                
                if(sum(valid_sums == n_pa) <= 1){
                    if(is.null(X_train) | method == "empirical"){
                        return(pop_marginal_pi[node_id])
                    }else{
                        return(NA)
                    }
                }
                
                if(is.null(X_train) | method == "empirical"){
                    return(mean(Y_train[valid_sums == n_pa, node_id]))
                }else if(all(Y_train[valid_sums == n_pa, node_id, drop = FALSE] == 1) | all(Y_train[valid_sums == n_pa, node_id, drop = FALSE] == 0)){
                    return(ExtremeClassifier(X_train[valid_sums == n_pa,], Y_train[valid_sums == n_pa, node_id, drop = FALSE]))
                }else if(method == "logistic"){
                    return(glm(y~., family = binomial, data = data.frame(y = Y_train[valid_sums == n_pa, node_id, drop = FALSE], x = X_train[valid_sums == n_pa,,drop=FALSE ])))
                }else if(method == "svm"){
                    if(is.matrix.csr(X_train)){
                        return(LinearSVM(X_train[valid_sums == n_pa,], Y_train[valid_sums == n_pa, node_id, drop = FALSE]))
                    }else{
                        return(LinearSVM(X_train[valid_sums == n_pa,,drop=FALSE ], Y_train[valid_sums == n_pa, node_id, drop = FALSE]))
                    }
                }
            })
            
            # pr(i | j); j in pa(i); used for the neighborhood approximation
            tmp_rlt_single <- lapply(pa_id, function(j){
                
                if(mean(Y_train[, j]) == 0 | mean(Y_train[, j]) == 1){
                    if(is.null(X_train) | method == "empirical"){
                        return(c(pop_marginal_pi[node_id], pop_marginal_pi[node_id])) #might not be correct
                    }else{
                        return(NA)
                    }
                }
                
                
                if(is.null(X_train) | method == "empirical"){
                    return(c(mean(Y_train[Y_train[, j] == 1, node_id]),
                             mean(Y_train[Y_train[, j] == 0, node_id])))
                }else{
                    tmp_list <- list()
                    for(tmp_val in c(1, 0)){
                        if(all(Y_train[Y_train[, j] == tmp_val, node_id] == 1) | all(Y_train[Y_train[, j] == tmp_val, node_id] == 0)){
                            tmp_list <- c(tmp_list, list(ExtremeClassifier(X_train[Y_train[, j] == tmp_val,], Y_train[Y_train[, j] == tmp_val, node_id])))
                        }else if(method == "logistic"){
                            tmp_list <- c(tmp_list, list(glm(y~., family = binomial, data = data.frame(y = Y_train[Y_train[, j] == tmp_val, node_id], x = X_train[Y_train[, j] == tmp_val,]))))
                        }else if(method == "svm"){
                            tmp_list <- c(tmp_list, list(LinearSVM(X_train[Y_train[, j] == tmp_val,], Y_train[Y_train[, j] == tmp_val, node_id])))
                        }
                    }
                    return(tmp_list)
                } 
            })
            
            #if(node_id == 4){browser()}
            
            #if(sum(lapply(tmp_rlt_joint, function(l) is.na(l)) %>% unlist) > 0 | sum(lapply(tmp_rlt_single, function(l) is.na(l)) %>% unlist)) browser()
            
            if(method == "empirical"){
                return(list("joint" = unlist(tmp_rlt_joint), "single" = unlist(tmp_rlt_single)))
            }else if(method %in% c("logistic", "svm")){
                return(list("joint" = tmp_rlt_joint, "single" = unlist(tmp_rlt_single, recursive = FALSE)))
            }
            
        }
    }
    
    
    pred_cdensity <- function(node_id, n_X, X = NULL){
        #pa_id <- as.numeric(strsplit(graph_structure[node_id], split = ";")[[1]])
        n_pa <- length(graph_structure[[node_id]])
        object_list <- cond_rlt[[node_id]]
        
        # if(all(pa_id == 0)){
        #     object_list <- list(cond_rlt[[node_id]])
        # }else{
        #     object_list <- lapply(seq(2^(n_pa)), function(j){
        #         y_val <- as.integer(intToBits(j - 1))[seq(n_pa)]
        #         obj_idx <- sum(y_val * 2^(seq(length(y_val) - 1)))
        #         cond_rlt[[node_id]][obj_idx]
        #     })
        # }
        
        if(all(graph_structure[[node_id]] == 0)){
            rlt <- list()
            if(is.null(X_train) | method == "empirical"){
                rlt[["joint"]] <- matrix(rep(object_list, n_X), nrow = n_X)  
            }else if(class(object_list) == "ExtremeClassifierClass"){
                rlt[["joint"]] <- predict.ExtremeClassifierClass(object_list, newdata = X) %>% matrix(nrow = n_X)
            }else if(method == "logistic"){
                rlt[["joint"]] <- predict.glm(object_list, newdata = data.frame(x = X), type = "response") %>% matrix(nrow = n_X)
            }else if(method == "svm"){
                rlt[["joint"]] <- attributes(predict(object_list, newdata = X, probability = TRUE))$probabilities[, "1"] %>% matrix(nrow = n_X)
            }
        }else{
            rlt <- list()
            try_rlt <- try({
                rlt[["joint"]] <- sapply(object_list$joint, function(object){
                    if(!is.na(object)){
                        if(is.null(X_train) | method == "empirical"){
                            return(rep(object, n_X))
                        }else if(class(object) == "ExtremeClassifierClass"){
                            return(predict.ExtremeClassifierClass(object, newdata = X)) 
                        }else if(method == "logistic"){
                            return(predict.glm(object, newdata = data.frame(x = X), type = "response"))
                        }else if(method == "svm"){
                            return(attributes(predict(object, newdata = X, probability = TRUE))$probabilities[, "1"])
                        }
                    }else{
                        return(rep(pop_marginal_pi[node_id], n_X)) # can be improved to mean response given pa == 1
                    }
                }) %>% matrix(nrow = n_X)
                rlt[["single"]] <- sapply(object_list$single, function(object){
                    if(!is.na(object)){
                        if(is.null(X_train) | method == "empirical"){
                            return(rep(object, n_X))
                        }else if(class(object) == "ExtremeClassifierClass"){
                            return(predict.ExtremeClassifierClass(object, newdata = X)) 
                        }else if(method == "logistic"){
                            return(predict.glm(object, newdata = data.frame(x = X), type = "response"))
                        }else if(method == "svm"){
                            return(attributes(predict(object, newdata = X, probability = TRUE))$probabilities[, "1"])
                        }
                    }else{
                        return(rep(pop_marginal_pi[node_id], n_X)) # can be improved to mean response given pa == 1
                    }
                }) %>% matrix(nrow = n_X)
            })
            if(class(try_rlt) == "try-error") browser()
        }
        return(rlt)
    }
    
    return(pred_cdensity)
}


Threshold <- function(S, threshold){
    return(t(as.numeric(t(S) - threshold >= 0)))
}

Clip <- function(val, lower = 0.1, upper = 0.9){
    val * as.numeric(val >= lower & val <= upper) + lower * as.numeric(val < lower) + upper * as.numeric(val > upper)
}

#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
LearnPr_S_Y <- function(S_train, Y_train, X_train, method = "gaussian", fit_model = "logistic"){
    n_class <- ncol(S_train)
    
    rlt <- foreach::foreach(node_id = seq(n_class)) %dopar%{
        tmp <- lapply(c(1, 0), function(val){
            if(is.null(X_train)){
                if(method == "binomial"){
                    mean(S_train[Y_train[, node_id] == val, node_id])
                }else if(method == "gaussian"){
                    c(mean(S_train[Y_train[, node_id] == val, node_id]), sd(S_train[Y_train[, node_id] == val, node_id]))
                }
            }else{
                if(method == "binomial"){
                    if(all(S_train[Y_train[, node_id] == val, node_id] == 1) | all(S_train[Y_train[, node_id] == val, node_id] == 0)){
                        ExtremeClassifier(X_train, S_train[Y_train[, node_id] == val, node_id])
                    }else if(fit_model == "logistic"){
                        glm(y~., family = binomial,
                            data = dataframe(y = S_train[Y_train[, node_id] == val, node_id],
                                             x = X_train[Y_train[, node_id] == val, node_id]))
                    }else if(fit_model == "svm"){
                        LinearSVM(X_train[Y_train[, node_id] == val, node_id], S_train[Y_train[, node_id] == val, node_id])
                    }
                }else if(method == "gaussian"){
                    c(mean(S_train[Y_train[, node_id] == val, node_id]), sd(S_train[Y_train[, node_id] == val, node_id]))
                }
            }
        })
        names(tmp) <- c("1", "0")
        return(tmp)
    }
    
    pred_pr_s_y <- function(nodes, S_test, X_test){
        if(is.null(nodes)){
            nodes <- seq(n_class)
        }
        
        lapply(nodes, function(node_id){
            if(is.null(X_test)){
                if(method == "binomial"){
                    cbind(abs(1 - S_test[,node_id] - rlt[[node_id]][["1"]]),
                          abs(1 - S_test[,node_id] - rlt[[node_id]][["0"]]))
                }else if(method == "gaussian"){
                    cbind(pnorm(S_test[,node_id], mean = rlt[[node_id]][["1"]][1], sd = rlt[[node_id]][["1"]][2]),
                          pnorm(S_test[,node_id], mean = rlt[[node_id]][["0"]][1], sd = rlt[[node_id]][["0"]][2]))
                }
            }else{
                if(method == "binomial"){
                    sapply(c(1, 0), function(val){
                        if(class(rlt[[node_id]][[as.character(val)]]) == "ExtremeClassifier"){
                            abs(1 - S_test[, node_id] - predict.ExtremeClassifierClass(rlt[[node_id]][[as.character(val)]], X_test)) 
                        }else if(fit_model == "logistic"){
                            abs(1 - S_test[,node_id] - predict.glm(rlt[[node_id]][[as.character(val)]], newdata = data.frame(x = X_test), type = "response"))
                        }else if(fit_model == "svm"){
                            abs(1 - S_test[,node_id] - attributes(predict(rlt[[node_id]][[as.character(val)]], newdata = X_test, probability = TRUE))$probabilities[, "1"])
                        }
                    })
                }else if(method == "gaussian"){
                    cbind(pnorm(S_test[,node_id], mean = rlt[[node_id]][["1"]][1], sd = rlt[[node_id]][["1"]][2]),
                          pnorm(S_test[,node_id], mean = rlt[[node_id]][["0"]][1], sd = rlt[[node_id]][["0"]][2]))
                }
            }
        })
    }
    
    return(pred_pr_s_y)
}

#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
LearnPr_y <- function(graph_structure, n_test, cdensity_test, Y_train, X_train, method = "hierarchical", fit_model = "logistic"){
    n_class <- ncol(Y_train)
    #n_test <- length(cdensity_test[[1]])
    
    #all_ancestors <- apply(AllAncestorsList_cpp(lapply(graph_structure, function(pa) pa - 1)), 2, function(x) which(x == 1))
    all_ancestors <- AllAncestorList(graph_structure)
    if(!is.list(all_ancestors)){
        all_ancestors <- lapply(seq(n_class), function(j) all_ancestors[, j])
    }
    rlt <- foreach::foreach(node_id = seq(n_class)) %dopar%{
        if(method == "learning"){
            if(is.null(X_train)){
                return(mean(Y_train[, node_id]))
            }else if(all(Y_train[, node_id] == 1) | all(Y_train[, node_id] == 0)){
                return(ExtremeClassifier(X_train, Y_train[, node_id]))
            }else{
                if(fit_model == "logistic"){
                    return(glm(y~., family = binomial, data = data.frame(x = X_train, y = Y_train[, node_id])))
                }else if(fit_model == "svm"){
                    return(LinearSVM(X_train, Y_train[, node_id]))
                }
            }
        }else if(method == "hierarchical"){
            sapply(seq(n_test), function(i){
                sapply(all_ancestors[[node_id]], function(j) cdensity_test[[j]][["joint"]][i, 1]) %>% prod() #be careful
            }) 
        }
    }
    
    pred_pr_y <- function(nodes, X_test){
        if(is.null(nodes)){
            nodes <- seq(n_class)
        }
        
        if(method == "empirical"){
            if(is.null(X_test)){
                return(sapply(nodes, function(j) rep(rlt[j], n_test)))
            }else{
                return(sapply(nodes, function(j){
                    if(class(rlt[[j]]) == "ExtremeClassifier"){
                        predict.ExtremeClassifierClass(rlt[[j]], X_test)
                    }else if(fit_model == "logistic"){
                        predict.glm(rlt[[j]], newdata = data.frame(x = X_test), type = "response")
                    }else if(fit_model == "svm"){
                        attributes(predict(rlt[[j]], newdata = X_test, probability = TRUE))$probabilities[, "1"]
                    }
                }))
            }
        }else if(method == "hierarchical"){
            return(sapply(nodes, function(j) rlt[[j]]))
        }
    }
    return(pred_pr_y)
}


AllAncestorList <- function(parent_id){
    lapply(seq(length(parent_id)), function(index_node){
        ancestors <- index_node
        queue <- parent_id[[index_node]]
        while(any(queue != 0)&&(length(queue)>0)){
            pointer <- queue[1]
            queue <- queue[-1]
            if(pointer != 0){
                ancestors <- unique(c(ancestors, pointer))
                new_mems <- parent_id[[pointer]]
                queue <- c(queue, new_mems)
            }else{
                queue <- c(queue, pointer)
            }
        }
        sort(ancestors)
    })
}

replicateList <- function(parents, num){
    final_list <- parents
    times <- 1
    while(times < num){
        to_add <- unlist(parents) + times*length(parents)
        to_add[unlist(parents) == 0] <- 0
        to_add <- relist(to_add, parents)
        final_list <- append(final_list, to_add)
        times <- times + 1
    }
    return(final_list)
}

LinearSVM <- function(train_X, train_Y){
    param_set <- expand.grid(cost = exp(seq(-5, 5)),
                             gamma = c(1),
                             stringsAsFactors = FALSE)
    
    cv_rlt <- sapply(seq(nrow(param_set)), function(param_idx){
        cost = param_set$cost[param_idx]
        gamma = param_set$gamma[param_idx]
        
        S_cv <- e1071::svm(x = train_X,
                           y = as.factor(train_Y),
                           #scale = TRUE,
                           type = "C-classification",
                           kernel = "linear",
                           cost = cost,
                           #gamma = gamma,
                           #cross = 5,
                           class.weights = 100/table(train_Y),
                           tolerance = 0.001,
                           probability = TRUE)
        sum(diag(table(S_cv$fitted, train_Y)))/length(train_Y)
    })
    
    best_param_idx <- which.max(cv_rlt)
    best_cost <- param_set$cost[best_param_idx]
    best_gamma <- param_set$gamma[best_param_idx]
    
    svm_train <- e1071::svm(x = train_X,
                            y = as.factor(train_Y),
                            #scale = TRUE,
                            type = "C-classification",
                            kernel = "linear",
                            cost = best_cost,
                            #gamma = best_gamma,
                            #cross = 5,
                            class.weights = 100/table(train_Y),
                            probability = TRUE)
    
    return(svm_train)
    
}



ExtremeClassifier <- function(X, Y, ...) {
    extreme_y = ifelse(all(Y == 1), 1, 0)
    model = structure(list(x = X, y = extreme_y), 
                      class = "ExtremeClassifierClass") 
    return(model)
}

# create a method for function print for class myClassifierClass
predict.ExtremeClassifierClass = function(modelObject, newdata, ...) {
    n_newdata = nrow(newdata)
    return(rep(modelObject$y, n_newdata))
}


