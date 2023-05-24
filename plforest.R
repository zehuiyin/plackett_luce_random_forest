library(tidyverse)
library(PlackettLuce) # for pltree function
library(rsample) # for taking samples
library(wrapr) # for sort vector
library(foreach) # for parallel computing
library(doParallel)
library(caret) # for predict function

###############################################################################
###############################################################################
# Create bagging or random forest based on pltree from package PlackettLuce
plforest <- function(rankings, covariates,
                     bootstrapping_n = 10, k_fold_CV = 5,
                     grouped_rankings = FALSE,
                     ...) {
  # bootstrapping covariate dataframe
  training_boots <- bootstraps(covariates, 
                               times = bootstrapping_n)
  
  # create a table to store the best tree for each bootstrapping iteration
  best_tree_table <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(best_tree_table) <- c("n", "k", "accurate_rate", "tree_name")
  
  # loop over each bootstrapping iteration
  for (n in 1:bootstrapping_n) {
    # subset the training dataset to n-th bootstrapping iteration
    subsample <- covariates[training_boots$splits[[n]]$in_id,]
    subsample_rank <- rankings[training_boots$splits[[n]]$in_id,]
    
    # perfrom k-fold CV
    fold <- vfold_cv(subsample, v = k_fold_CV)
    
    # create a table to store the accurate rate of k-fold tree iterations
    accurate_rate_table <- as.data.frame(matrix(nrow = 0, ncol = 4))
    colnames(accurate_rate_table) <- c("n", "k", "accurate_rate", "tree_name")
    
    # loop over each n-fold CV iteration
    for (k in 1:k_fold_CV) {
      # subset the training dataset to k-th k-fold CV iteration
      subsample_fold <- subsample[fold$splits[[k]]$in_id,]
      subsample_rank_fold <- subsample_rank[fold$splits[[k]]$in_id,]
      
      # convert rankings to group rankings format if needed
      if (grouped_rankings == FALSE) {
        subsample_rank_fold <- group(subsample_rank_fold, 
                                     index = 1:length(subsample_rank_fold))
      }
      
      # perform Plackett-Luce trees
      # ... for additional arguments pass to mob_control function
      subsample_rank_fold_tree <- pltree(formula = subsample_rank_fold ~.,
                                         data = subsample_fold, ...)
      
      # subset the validating dataset to k-th k-fold CV iteration
      subsample_fold_test <- subsample[-fold$splits[[k]]$in_id,]
      
      # generate prediction based on the fitted tree
      subsample_fold_predict <- 
        predict(subsample_rank_fold_tree, subsample_fold_test)
      
      # create a table to hold predicted rankings
      subsample_fold_predict_rank <- as.data.frame(
        matrix(nrow = nrow(subsample_fold_test), 
               ncol = dim(as.rankings(rankings))[2]))
      
      # add results to the predicted rankings table
      for (i in 1:nrow(subsample_fold_predict_rank)) {
        sorted_colnames <- subsample_fold_predict[i,] |>
          sort(decreasing = T) |>
          as.data.frame() |>
          rownames()
        for (c in 1:ncol(subsample_fold_predict_rank)) {
          subsample_fold_predict_rank[i, c] <- sorted_colnames[c]
        }
      }
      
      # convert predicted rankings table to rankings format
      subsample_fold_predict_rank <- 
        as.rankings(subsample_fold_predict_rank, input = "orderings")
      
      # subset to create true rankings for validating subsample
      subsample_fold_true_rank <-
        subsample_rank[-fold$splits[[k]]$in_id,]
      
      # check for whether prediction match the true rankings
      checker <-
        subsample_fold_predict_rank == subsample_fold_true_rank
      
      # create a counter to store the number of matching
      correct_case <- 0
      
      # generate the correct case count
      for (i in 1:nrow(subsample_fold_predict_rank)) {
        if (all(checker[i,])) {
          correct_case <- correct_case + 1
        }
      }
      
      # calculate the accurate rate for this fold iteration
      accurate_rate <- correct_case/nrow(subsample_fold_predict_rank)*100
      
      # store the iteration results and the accurate_rate
      assign(paste0(n, "_", k, "_tree"), subsample_rank_fold_tree)
      additional_row <- cbind(n, 
                              k, 
                              accurate_rate, 
                              paste0(n, "_", k, "_tree")) |> as.data.frame()
      additional_row |> colnames() <-
        c("n", "k", "accurate_rate", "tree_name")
      accurate_rate_table <- rbind(accurate_rate_table,
                                   additional_row)
    }
    # find the tree name with the highest accuracy rate
    best_tree_name <- 
      accurate_rate_table[which(accurate_rate_table$accurate_rate ==
                                  max(accurate_rate_table$accurate_rate)),
                          "tree_name"]
    
    # handle the case when two or more trees have tied accuracy
    # return only the first one for simplicity
    best_tree_name <- best_tree_name[1]
    
    # store the iteration results and the accurate_rate
    additional_row <- cbind(n, k, 
                            max(accurate_rate_table$accurate_rate), 
                            best_tree_name) |> as.data.frame()
    additional_row |> colnames() <-
      c("n", "k", "accurate_rate", "tree_name")
    best_tree_table <- rbind(best_tree_table,
                             additional_row)
    
    # print a counter to show percentage of completion for each bootstrapping
    cat('\b\b\b\b\b\b\b\b\b\b\b\b\b',
        n/bootstrapping_n*100,
        "% Completed", sep = "")
  }
  # combine outputs into a large list for returning
  output <- list(best_tree_table, mget(best_tree_table$tree_name))
  
  return(output)
}

###############################################################################
###############################################################################
# Predict worth based on new covariates according to the plforest model
predict_plforest <- function(plforest_output, covariates) {
  # record the bootstrapping number of the output
  n <- dim(plforest_output[[1]])[1]
  
  # record the number of ranked levels
  model_predict_n <- predict(plforest_output[[2]][1], covariates)
  model_predict_n <- as.data.frame(model_predict_n[[1]])
  n_level <- ncol(model_predict_n)
  
  # create a holder to store the bagging predicted results
  bagged_predict <- as.data.frame(matrix(0,
                                         nrow = nrow(covariates), 
                                         ncol = n_level))
  
  # store the level names to the holder
  colnames(bagged_predict) <- colnames(model_predict_n)
  
  # convert the holder into matrix for summation
  bagged_predict <- as.matrix(bagged_predict)
  
  # loop over each bootstrapping iteration's best model prediction
  # and summing the predicted worth up
  for (i in 1:n) {
    model_predict_n <- predict(plforest_output[[2]][i], covariates)
    model_predict_n <- model_predict_n[[1]]
    bagged_predict <- bagged_predict + 
      as.matrix(model_predict_n)
  }
  
  # create a holder to store the bagging predicted results in rank order
  bagged_predict_rank <- as.data.frame(matrix(nrow = nrow(bagged_predict),
                                              ncol = ncol(bagged_predict)))
  colnames(bagged_predict_rank) <- 1:ncol(bagged_predict_rank)
  
  # loop over each prediction to convert worth to ranking
  for (i in 1:nrow(bagged_predict)) {
    rank <- as.numeric(bagged_predict[i,])
    names(rank) <- colnames(bagged_predict)
    bagged_predict_rank[i, ] <-
      names(sort(rank, decreasing = T))
  }
  
  # convert rank ordering in rankings format
  bagged_predict_rank <- as.rankings(bagged_predict_rank, input = "ordering")
  
  # standardize the predicted worth matrix (enforce each row sum to 1)
  bagged_predict <- bagged_predict/n
  
  # combine the predicted worth matrix and predicted rankings together into 
  # a list for return
  output <- list(bagged_predict, bagged_predict_rank)
  
  return(output)
}

###############################################################################
###############################################################################
# A parallel computing implementation of plforest function
plforest_mc <- function(rankings, covariates,
                        bootstrapping_n = 10, k_fold_CV = 5,
                        grouped_rankings = FALSE, cores = detectCores()-1,
                        ...) {
  # bootstrapping covariate dataframe
  training_boots <- bootstraps(covariates, 
                               times = bootstrapping_n)
  
  # create a table to store the best tree for each bootstrapping iteration
  best_tree_table <- as.data.frame(matrix(nrow = bootstrapping_n,
                                          ncol = 4))
  colnames(best_tree_table) <- c("n", "k", "accurate_rate", "tree_name")
  
  # create or empty the progress txt file to denote progress
  close(file("./progress.txt", open="w"))
  
  # register cluster
  cl <- makeCluster(cores) # not to overload computer
  registerDoParallel(cl)
  
  # loop over each bootstrapping iteration
  # with parallel computing
  best_tree_list <- 
    foreach(n=1:bootstrapping_n, .combine = rbind,
            .packages = c("tidyverse", "PlackettLuce",
                          "rsample", "wrapr")) %dopar% 
    {
      # subset the training dataset to n-th bootstrapping iteration
      subsample <- covariates[training_boots$splits[[n]]$in_id,]
      subsample_rank <- rankings[training_boots$splits[[n]]$in_id,]
                            
      # perfrom k-fold CV
      fold <- vfold_cv(subsample, v = k_fold_CV)
                            
      # create a table to store the accurate rate of k-fold tree iterations
      accurate_rate_table <- as.data.frame(matrix(nrow = 0,
                                                  ncol = 4))
      colnames(accurate_rate_table) <- c("n", "k", "accurate_rate", "tree_name")
                            
      # loop over each n-fold CV iteration
      for (k in 1:k_fold_CV) {
        # subset the training dataset to k-th k-fold CV iteration
        subsample_fold <- subsample[fold$splits[[k]]$in_id,]
        subsample_rank_fold <- subsample_rank[fold$splits[[k]]$in_id,]
                              
        # convert rankings to group rankings format if needed
        if (grouped_rankings == FALSE) {
          subsample_rank_fold <- group(subsample_rank_fold, 
                                       index = 1:length(subsample_rank_fold))
          }
                              
        # subset the validating dataset to k-th k-fold CV iteration
        subsample_fold_test <- subsample[-fold$splits[[k]]$in_id,]
                              
        # loop over all the factor class variables
        # add missing levels in the training dataset from the testing dataset
        for (x in 1:ncol(subsample_fold)) {
          if (class(subsample_fold[,x]) == "factor") {
            if (all(levels(subsample_fold[,x]) == 
                    levels(subsample_fold_test[,x]))) {
              next
            } else {
              all_levels <- union(levels(subsample_fold[,x]),
                                  levels(subsample_fold_test[,x]))
              levels(subsample_fold[,x]) <- all_levels
              levels(subsample_fold_test[,x]) <- all_levels
            }
          } else {
            next
          }
        }
        
        # perform Plackett-Luce trees
        # ... for additional arguments pass to mob_control function
        subsample_rank_fold_tree <- pltree(formula = subsample_rank_fold ~.,
                                           data = subsample_fold, ...)
        
        # generate prediction based on the fitted tree
        subsample_fold_predict <- 
          predict(subsample_rank_fold_tree, subsample_fold_test)
                              
        # create a table to hold predicted rankings
        subsample_fold_predict_rank <- as.data.frame(
          matrix(nrow = nrow(subsample_fold_test), 
                 ncol = dim(as.rankings(rankings))[2]))
                              
        # add results to the predicted rankings table
        for (i in 1:nrow(subsample_fold_predict_rank)) {
          sorted_colnames <- subsample_fold_predict[i,] |>
            sort(decreasing = T) |>
            as.data.frame() |>
            rownames()
          for (c in 1:ncol(subsample_fold_predict_rank)) {
            subsample_fold_predict_rank[i, c] <- sorted_colnames[c]
            }
          }
                              
        # convert predicted rankings table to rankings format
        subsample_fold_predict_rank <- 
          as.rankings(subsample_fold_predict_rank, input = "orderings")
                              
        # subset to create true rankings for validating subsample
        subsample_fold_true_rank <-
          subsample_rank[-fold$splits[[k]]$in_id,]
                              
        # check for whether prediction match the true rankings
        checker <-
          subsample_fold_predict_rank == subsample_fold_true_rank
                              
        # create a counter to store the number of matching
        correct_case <- 0
                              
        # generate the correct case count
        for (i in 1:nrow(subsample_fold_predict_rank)) {
          if (all(checker[i,])) {
            correct_case <- correct_case + 1
            }
          }
                              
        # calculate the accurate rate for this fold iteration
        accurate_rate <- correct_case/nrow(subsample_fold_predict_rank)*100
                              
        # store the iteration results and the accurate_rate
        assign(paste0(n, "_", k, "_tree"), subsample_rank_fold_tree)
        additional_row <- cbind(n, 
                                k, 
                                accurate_rate, 
                                paste0(n, "_", k, "_tree")) |> 
          as.data.frame()
        additional_row |> colnames() <-
          c("n", "k", "accurate_rate", "tree_name")
        accurate_rate_table <- rbind(accurate_rate_table,
                                     additional_row)
        }
      # find the tree name with the highest accuracy rate
      best_tree_name <- accurate_rate_table[
        which(accurate_rate_table$accurate_rate ==
                max(accurate_rate_table$accurate_rate)),
                                            "tree_name"]
                            
      # handle the case when two or more trees have tied accuracy
      # return only the first one for simplicity
      best_tree_name <- best_tree_name[1]
                            
      # write a counter check progress
      # add a dot to the progress txt file in the directory
      write(".", file = "./progress.txt", append = T, sep = "")
                            
      # store the iteration results and the accurate_rate
      additional_row <- cbind(n, k, 
                              max(accurate_rate_table$accurate_rate), 
                              best_tree_name) |> as.data.frame()
      additional_row |> colnames() <-
        c("n", "k", "accurate_rate", "tree_name")
                            
      # list together the iteration results for returning
      list(additional_row, get(additional_row$tree_name))
      }
  
  # stop parallel cluster
  stopCluster(cl)
  
  # record the list results into the table for return
  for (i in 1:nrow(best_tree_table)) {
    best_tree_table[i,] <- best_tree_list[[i,1]]
  }
  
  # combine outputs into a large list for returning
  output <- list(best_tree_table, best_tree_list[,2])
  
  return(output)
}

###############################################################################
###############################################################################
# Partition large list data into smaller files for data transfer
split.file <- function(db, rows, basename) {
  n = nrow(db)
  m = n %/% rows
  for (k in seq_len(m)) {
    db.sub <- db[seq(1 + (k-1)*rows, k*rows), , drop = F]
    saveRDS(db.sub, file = sprintf("%s%.5d.rds", basename, k),
            compress = "xz", ascii = F)
  }
  if (m * rows < n) {
    db.sub <- db[seq(1 + m*rows, n), , drop = F]
    saveRDS(db.sub, file = sprintf("%s%.5d.rds", basename, m+1),
            compress = "xz", ascii = F)
    m <- m + 1
  }
  m
}

###############################################################################
###############################################################################
# Join the partitioned files together into 1
join.files <- function(basename) {
  files <- sort(list.files(pattern = sprintf("%s[0-9]{5}\\.rds", basename)))
  do.call("rbind", lapply(files, readRDS))
}
