# Creating a stump for a tree
stump <- function(data){

  # Creating the base node
  node <- list()
  node[["node0"]] <- list(
    # Creating the node number
    node_number = 0,
    isRoot = TRUE,
    # Creating a vector with the tranining index
    train_index = 1:nrow(data$x_train),
    test_index = 1:nrow(data$x_test),
    depth_node = 0,
    node_var = NA,
    node_cutpoint_index = NA,
    left = NA,
    right = NA,
    parent_node = NA,
    terminal = TRUE,
    gamma = 0,
    betas_vec = matrix(0,nrow = ncol(data$B_train_arr),ncol = data$n_tree)
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
  non_terminal <- names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)]

  # In case there are non nonterminal nondes
  if(length(non_terminal)==0){
    return(non_terminal)
  }

  bool_nog <- vector("logical",length = length(non_terminal))
  for(i in 1:length(bool_nog)){
    # Checking if both children are terminal
    if( tree[[tree[[non_terminal[i]]]$left]]$terminal & tree[[tree[[non_terminal[i]]]$right]]$terminal) {
      bool_nog[i] <- TRUE
    }
  }

  return(  non_terminal[bool_nog])
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
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0
  while(valid_terminal_node){
      # Convinience while to avoid terminal nodes of 2
      # Sample a split var
      p_var <- sample(1:ncol(data$B_train_arr),size = 1)

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

      if( (length(left_index) > 1) & (length(right_index)>1)){
          # Getting out of the while
          break
      } else {

          # Adding one to the counter
          valid_count = valid_count + 1

          # Stop trying to search for a valid cutpoint
          if(valid_count > 2) {
            valid_terminal_node = FALSE
            return(tree)
          }
      }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }
  # Calculating loglikelihood for the grown node, the left and the right node

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           B_train_arr = data$B_train_arr,
                           index_node = g_node$train_index,
                           tau_beta_vec = data$tau_beta_vec,
                           P = data$P)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               B_train_arr = data$B_train_arr,
                               index_node = left_index,
                               tau_beta_vec = data$tau_beta_vec,
                               P = data$P)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               B_train_arr = data$B_train_arr,
                               index_node = right_index,
                               tau_beta_vec = data$tau_beta_vec,
                               P = data$P)

  # Calculating the prior
  prior_loglike <- log(alpha*(1+g_node$depth_node)^(-beta)) + # Prior of the grown node becoming nonterminal
                   2*log(1-alpha*(1+g_node$depth_node+1)^(-beta)) - # plus the prior of the two following nodes being terminal
                   log(1-alpha*(1+g_node$depth_node)^(-beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(0<acceptance){
    left_node <- list(node_number = max_index+1,
                      isRoot = FALSE,
                      train_index = left_index,
                      test_index = left_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      right = NA,
                      parent_node = g_node_name,
                      terminal = TRUE,
                      gamma = 0,
                      betas_vec = g_node$betas_vec)

    names(tree[[1]])

    right_node <- list(node_number = max_index+2,
                      isRoot = FALSE,
                      train_index = right_index,
                      test_index = right_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      right = NA,
                      parent_node = g_node_name,
                      terminal = TRUE,
                      gamma = 0,
                      betas_vec = g_node$betas_vec)

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


# Pruning a tree
prune <- function(tree,
                 curr_part_res,
                 data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  children_left_index <- tree[[p_node$left]]$train_index
  children_right_index <- tree[[p_node$right]]$train_index

  # Calculating loglikelihood for the grown node, the left and the right node

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           B_train_arr = data$B_train_arr,
                           index_node = p_node$train_index,
                           tau_beta_vec = data$tau_beta_vec,
                           P = data$P)


  p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               B_train_arr = data$B_train_arr,
                               index_node = children_left_index,
                               tau_beta_vec = data$tau_beta_vec,
                               P = data$P)

  p_right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                B_train_arr = data$B_train_arr,
                                index_node = children_right_index,
                                tau_beta_vec = data$tau_beta_vec,
                                P = data$P)

  # Calculating the prior
  prior_loglike <- log(1-alpha*(1+p_node$depth_node)^(-beta)) - # Prior of the new terminal node
    log(alpha*(1+p_node$depth_node)^(-beta)) - # Prior of the grown node becoming nonterminal
    2*log(1-alpha*(1+p_node$depth_node+1)^(-beta))  # plus the prior of the two following nodes being terminal
    # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(p_loglike-p_left_loglike-p_right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]] <- NULL
    tree[[p_node$right]] <- NULL

    # Modifying back the pruned node
    tree[[p_node_name]]$left <- NA
    tree[[p_node_name]]$right <- NA
    tree[[p_node_name]]$terminal <- TRUE

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}


# Change a tree
change <- function(tree,
                 curr_part_res,
                 data){

  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0


  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$B_train_arr),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(x_train[c_node$train_index,p_var])

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

    # Avoiding having terminal nodes with just one observation
    if( (length(left_index) > 1) & (length(right_index)>1)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }
  # Calculating loglikelihood for the new changed nodes and the old ones

  c_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                           B_train_arr = data$B_train_arr,
                           index_node = tree[[c_node$left]]$train_index,
                           tau_beta_vec = data$tau_beta_vec,
                           P = data$P)


  c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                               B_train_arr = data$B_train_arr,
                               index_node = tree[[c_node$right]]$train_index,
                               tau_beta_vec = data$tau_beta_vec,
                               P = data$P)

  new_c_loglike_left <-  nodeLogLike(curr_part_res = curr_part_res,
                                B_train_arr = data$B_train_arr,
                                index_node = left_index,
                                tau_beta_vec = data$tau_beta_vec,
                                P = data$P)

  new_c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                     B_train_arr = data$B_train_arr,
                                     index_node = right_index,
                                     tau_beta_vec = data$tau_beta_vec,
                                     P = data$P)


  # Calculating the acceptance probability
  acceptance <- exp(new_c_loglike_left+new_c_loglike_right-c_loglike_left-c_loglike_right)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    # Updating the left and the right node
    # === Left =====
    tree[[c_node$left]]$node_var <- p_var
    tree[[c_node$left]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$left]]$train_index <- left_index
    tree[[c_node$left]]$test_index <- left_test_index

    #==== Right ====
    tree[[c_node$right]]$node_var <- p_var
    tree[[c_node$right]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$right]]$train_index <- right_index
    tree[[c_node$right]]$test_index <- right_test_index

  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)

}


# ==========================
# Updating the gamma
# ==========================
updateGamma <- function(tree,
                        curr_part_res,
                        data){

  # Getting the terminals
  t_nodes_names <- get_terminals(tree)
  basis_dim <- dim(data$B_train_arr)[3]
  knots_dim <- ncol(data$B_train_arr)


  # Creating each element for each basis
  basis_sum_list <- vector("list",basis_dim)


  for(i in 1:length(t_nodes_names)){

    cu_t <- tree[[t_nodes_names[i]]]

    # Updating the sum of the residuals
    r_sum = curr_part_res[cu_t$train_index]

    basis_sum <- numeric(basis_dim)

    for(j in 1:basis_dim){
      basis_sum[j] <- crossprod(cu_t$betas_vec[,j,drop = FALSE],colSums(data$B_test_arr[cu_t$train_index,,j]))
    }

    s_gamma_inv <- 1/(length(cu_t$train_index)+data$tau_gamma/data$tau)

    # Computing mean and sd from the intercep
    gamma_mean <- s_gamma_inv*(r_sum+sum(basis_sum))
    gamma_sd <- sqrt(s_gamma_inv/(data$tau))

    tree[[t_nodes_names[i]]]$gamma = stats::rnorm(n = 1,mean = gamma_mean,sd = gamma_sd)

  }

  return(tree)
}

# ============
# Update Betas
# ============

updateBetas <- function(tree,
                        curr_part_res,
                        data){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)
  basis_dim <- dim(data$B_train_arr)[3]
  knots_dim <- ncol(data$B_train_arr)


  # Creating each element for each basis
  basis_sum_list <- vector("list",basis_dim)


  for(i in 1:length(t_nodes_names)){

    cu_t <- tree[[t_nodes_names[i]]]

    for(j in 1:basis_dim){

      Gamma_beta_tau_chol <-  chol(crossprod(data$B_train_arr[cu_t$train,,j])+data$tau_beta_vec[j]/data$tau*data$P)
      Gamma_beta_inv <- chol2inv(Gamma_beta_tau_chol)

      sum_aux <- matrix(0,nrow = length(cu_t$train_index),ncol = 1)
      # Sum j exception
      for(k in (1:basis_dim)[-j]){
          sum_aux <- sum_aux + (data$B_train_arr[cu_t$train_index,,k]%*%cu_t$betas_vec[,k,drop = FALSE])
      }

      # Calculating the mean to be sampled
      beta_mean <- Gamma_beta_inv%*%crossprod(data$B_train_arr[cu_t$train_index,,j],curr_part_res[cu_t$train_index])-crossprod(data$B_train_arr[cu_t$train_index,,j],(cu_t$gamma+sum_aux))

      # Check this line again if there's any bug on the cholesky decomposition
      tree[[t_nodes_names[i]]]$betas_vec[,j] <- mvnfast::rmvn(n = 1,mu = beta_mean,sigma = (1/(data$tau))*Gamma_beta_inv,isChol = FALSE)

      # #If we want to use the cholesky decomposition, doesn;t seem any faster though
      # test1 <- mvnfast::rmvn(n = 100,mu = beta_mean,sigma = sqrt(data$tau)*Gamma_beta_tau_chol,isChol = TRUE),
    }

  }

}

# =================
# Update \tau_betas
# =================
update_tau_betas <- function(forest,
                             data){

  # Getting the tau_b size
  basis_size <- dim(data$B_test_arr)[3]
  knots_size <- ncol(data$B_train_arr)
  tau_beta_vec_aux <- numeric(basis_size)

  # Same default as the paper;
  nu <- 2

  tau_b_shape <- numeric(basis_size)
  tau_b_rate <- numeric(basis_size)

  # Iterating over all basis
  for(k in 1:basis_size){


      # Iterating over all trees
      for(i in 1:length(forest)){

        # Getting terminal nodes
        t_nodes_names <- get_terminals(forest[[i]])
        n_t_nodes <- length(t_nodes_names)

        # Iterating over the terminal nodes
        for(j in 1:length(t_nodes_names)){

          tau_b_shape[k] <- tau_b_shape[k] + n_t_nodes
          tau_b_rate[k] <- tau_b_rate[k] + crossprod(forest[[i]][[t_nodes_names[j]]]$betas_vec[,j,drop = FALSE],(data$P%*%forest[[i]][[t_nodes_names[j]]]$betas_vec[,j,drop = FALSE]))

        }


      }

    tau_beta_vec_aux[k] <- rgamma(n = 1,
                                   shape = 0.5*knots_size*tau_b_shape[k] + 0.5*nu,
                                   rate = 0.5*tau_b_rate[k] + 0.5*nu*data$delta_vec[k])
  }

  return(tau_beta_vec_aux)

}

# ===================
# Updating the \delta
# ===================

update_delta <- function(data){

  delta_aux <- numeric(dim(data$B_train_arr)[3])
  nu <- 2
  a_delta <- d_delta <- 0.001

  for(i in 1:dim(data$B_train_arr)[3]){

    delta_aux[i] <- stats::rgamma(n = 1,shape =  0.5*nu + a_delta,
                               rate = 0.5*nu*data$tau_beta_vec[i]+d_delta)
  }


  # Returning the delta sampled vector
  return(delta_aux)
}

