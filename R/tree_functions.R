# Creating a stump for a tree
stump <- function(x_train,
                  x_test,
                  B_train_arr,
                  B_test_arr){

  # Creating the base node
  node <- list()
  node[["node0"]] <- list(
    # Creating the node number
    node_number = 0,
    isRoot = TRUE,
    # Creating a vector with the tranining index
    train_index = 1:nrow(x_train),
    test_index = 1:nrow(x_test),
    depth_node = 0,
    node_var = NA,
    node_cutpoint_index = NA,
    left = NA,
    right = NA,
    parent_node = NA,
    terminal = TRUE,
    gamma = 0,
    B_train_arr = B_train_arr,
    B_test_arr  = B_test_arr
   )

  # Returning the node
  return(node)

}

# Get all the terminal nodes
get_terminals <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Get nog terminal nodes
get_nogs <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Getting the maximum node index number
get_max_node <- function(tree){

  # Return the name of the termianl nodes
  return(max(unlist(lapply(tree, function(x){x$node_number}),use.names =  TRUE)))
}


# A function to calculate the loglikelihood
nodeLogLike <- function(curr_part_res,
                        index_node,
                        B_train_arr,
                        tau_beta_vec,
                        P){

  # Matrix aux basis sum likelihood
  n_leaf <- length(index_node)
  basis_aux <- matrix(0,n_leaf,n_leaf)
  for(i in 1:dim(B_train_arr)[3]){
    basis_subset <- B_train_arr[index_node,,i]
    basis_aux <- basis_aux + (1/tau_beta_vec[i])*basis_subset%*%solve(P,t(basis_subset))
  }

  return(mvnfast::dmvn(X = curr_part_res[index_node],mu = rep(0,length(index_node)),
                sigma = diag(1/tau,nrow = n_leaf)+tau_mu+basis_aux,log = TRUE))

}

# Grow a tree
grow <- function(tree,
                 x_train,
                 curr_part_res){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]

  # Sample a split var
  p_var <- sample(1:ncol(tree[[1]]$B_train_arr),size = 1)

  # Selecting an available cutpoint from this terminal node
  valid_range_grow <- range(x_train[g_node$train_index,p_var])

  # Getting which cutpoints are valid and sample onde index
  sample_cutpoint <- sample(which(xcut_m[,p_var]>valid_range_grow[1] & xcut_m[,p_var]<valid_range_grow[2]),
         size = 1)

  # Getting the left & right index
  left_index  <- all_var_splits[[p_var]][[sample_cutpoint]]$left_train[all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% g_node$train_index]
  right_index  <- all_var_splits[[p_var]][[sample_cutpoint]]$right_train[all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% g_node$train_index]

  left_test_index  <- all_var_splits[[p_var]][[sample_cutpoint]]$left_test[all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% g_node$test_index]
  right_test_index  <- all_var_splits[[p_var]][[sample_cutpoint]]$right_test[all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% g_node$test_index]



  # Verifying that the correct number was used
  if((length(left_index)+length(right_index))!=length(g_node$train_index)){
    stop("Something went wrong here --- train grown index doest match")
  }

  if((length(left_test_index)+length(right_test_index))!=length(g_node$test_index)){
    stop("Something went wrong here --- test grown index doest match")
  }

  # Calculating loglikelihood for the grown node, the left and the right node

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           B_train_arr = B_train_arr,
                           index_node = g_node$train_index,
                           tau_beta_vec = tau_beta_vec,
                           P = P)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               B_train_arr = B_train_arr,
                               index_node = left_index,
                               tau_beta_vec = tau_beta_vec,
                               P = P)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               B_train_arr = B_train_arr,
                               index_node = right_index,
                               tau_beta_vec = tau_beta_vec,
                               P = P)

  # Calculating the prior
  prior_loglike <- log(alpha*(1+g_node$depth_node)^(-beta)) + # Prior of the grown node becoming nonterminal
                   2*log(1-alpha*(1+g_node$depth_node+1)^(-beta)) - # plus the prior of the two following nodes being terminal
                   log(1-alpha*(1+g_node$depth_node)^(-beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){
    left_node <- list(node_number = max_index+1,
                      isRoot = FALSE,
                      train_index = left_index,
                      test_index = left_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      rigth = NA,
                      parent_node = g_node_name,
                      terminal = TRUE,
                      gamma = 0)

    right_node <- list(node_number = max_index+2,
                      isRoot = FALSE,
                      train_index = right_index,
                      test_index = right_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      rigth = NA,
                      parent_node = g_node_name,
                      terminal = TRUE,
                      gamma = 0)

    # Modifying the current node
    tree[[g_node_name]]$left = paste0("node",max_index+1)
    tree[[g_node_name]]$right = paste0("node",max_index+2)
    tree[[g_node_name]]$terminal = FALSE

    tree[[paste0("node",max_index+1)]] <- left_node
    tree[[paste0("node",max_index+2)]] <- right_node


  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)
}


# Get all the NOGs (it parent of terminal children)
get_nog <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}


