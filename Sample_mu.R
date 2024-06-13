sample_mu <- function(tree, tmax){

  # we start with a tree and prior 
  #
  # ex: 0 1 2 0 0 
  #
  # first event has 0 (background) as its parent
  # 
  # 2nd event has 1 as its parent
  # 3rd event has 2 as its parent
  # 4th event has 0 as its parent
  # 5th event has 0 as its parent  
  # 
  # ------------------------------------------------------
  # inputs:
  #   - tree (vector): vector of N elements, where the ith element is the
  #                    value of its parent
  #   - tmax (number): tmax
  #                      
  # outputs:
  #   - mu_new (number): posterior sample of mu 
  # ------------------------------------------------------
    
  # M = number of immigrants
  M <- sum(tree == 0)
  
  # sample from the posterior: gamma(2M, 2T)
  mu_new <- rgamma(1, 2 * M, 2 * tmax)
  
  return(mu_new)
}



update_phi <- function(thetas, alphas, betas, total_clusters){
  
  # update the (alpha, beta) parameters for clusters that are not observed
  # 
  # new_beta ~ gamma(alpha_0, beta_0)
  # new_alpha ~ unif[0, new_beta]
  # 
  # input:
  #     - thetas (N-dim vector): membership vector for each observation
  #     - alphas (N-dim vector): alpha value for each observation
  #     - betas (N-dim vector): beta value for each observation
  #     - total_clusters (integer): number of total clusters in the DPMM
  #     
  # output:
  #     - params_distinct (total_clusters x 3 dataframe)
  #       - col1: thetas in ascending order
  #       - col2: alphas
  #       - col3: betas
  
  params <- data.frame(thetas, alphas, betas)
  
  params_distinct <- distinct(params) %>% 
    arrange(thetas) 
  
  # generate new parameters for clusters that are not observed
  
  missing_thetas <- setdiff(1:total_clusters, params_distinct$thetas)
  
  # for each cluster thats not observed, generate a new alpha, beta pair

  for(j in missing_thetas){
    alpha_0 <- 2
    beta_0 <- 3
    new_beta <- rgamma(1, alpha_0, beta_0)
    new_alpha <- runif(1) * new_beta
    new_row <- c(j, new_alpha, new_beta)
    
    params_distinct <- params_distinct <- rbind(params_distinct, new_row)
  }
  
  params_distinct <- params_distinct %>% arrange(thetas)
  
  return(params_distinct)
    
}



sample_memberships <- function(tree, 
                               events, thetas, alphas, betas,
                               total_clusters, alpha){
  

  # implementing Neal algorithm (2000)
  # 
  # 
  # we have a sampled tree
  # we have our events, and their memberships
  # 
  # for each event, resample its membership based on the likelihood described in 
  # neal algorithm 8 (2000)
  # 
  # ------------------------------------------------------------------------
  # input:
  #   - tree (N-dim vector): each element describes its parent
  #   - events (N-dim vector): vector of N observed events timestamps
  #   - thetas (N-dim vector): membership vector for each observation
  #   - alphas (N-dim vector): alpha value for each observation
  #   - betas (N-dim vector): beta value for each observation
  #   - total_clusters (number): number of total clusters pre-specified in the DPMM
  #   - alpha (number): discount parameter when creating new "tables"
  # 
  # output:
  #   - param_df (N x 3 dataframe): dataframe with N rows that describes the new
  #                                 membership and therefore alpha and beta values
  #                                 for each observation
  #                                 
  #                                 if you have no kids, you wont switch
  # -----------------------------------------------------------------------
  
  N <- length(events)

  
  # for each event (AKA parent)
  
  for(i in 1:N){
    
    # update k and m in case we have occupied a new table in the past iteration
    
    # k = number of occupied "tables"
    # m = number of unoccupied "tables"
    
    k <- length(unique(thetas))
    m <- total_clusters - k
    
    # update alphas and betas for clusters that are empty
    
    params_distinct <- update_phi(thetas, alphas, betas, total_clusters) %>% 
      arrange(thetas)
    
    alpha_dict <- params_distinct$alphas
    beta_dict  <- params_distinct$betas
    theta_dict <- params_distinct$thetas
    
    # obtain the times of its children 
    
    # if there are no children, SKIP? OR DO SOMETHING? IDK
    
    t_parent <- events[i]
    children_ID <- which(tree == i)
    
    
    # if there are no children, do nothing, else, proceed
    
    if(length(children_ID > 0)){ 
    
      t_children <- events[children_ID]
      
      t_children_normalized <- t_children - t_parent
      
      # calculate F(y_i, phi_c) for each cluster c
      
      loglik <- rep(0, total_clusters)
      
      for(j in 1:total_clusters){
        
        loglik[j] <- decay_exp_ll(alpha_dict[j], 
                                  beta_dict[j], 
                                  t_children_normalized)
        
      }
      
      # calculate the n_-i/(n-1+alpha) fraction
      
      n_frac <- rep(0, total_clusters)
      
      thetas_no_i <- thetas[-c(i)]
      
      for(j in 1:total_clusters){
        
        # if we are working with a cluster with members
        if(j %in% thetas_no_i){
          
          n_i <- sum(thetas_no_i == j)
          n_frac[j] <- (n_i/(N - 1 + alpha))
          
        } else{ # if we are working with a cluster with no members
          if(m == 0){
            n_frac[j] <- 0
          } else{
            n_frac[j] <- (alpha / m) / (N - 1 + alpha)
          }
          
        }
        
        
      }
      
      # multiply the alpha_frac with F(y_i, phi_c) and normalize
      
      # here i do an exponential trick when working with log likelihoods
      
      log_probs <- loglik + log(n_frac)
      
      max_log <- max(log_probs)
      
      sum_exp <- sum(exp(log_probs - max_log))
      
      log_sum_exp_value <- max_log + log(sum_exp)
      
      probs <- exp(log_probs - log_sum_exp_value)
      
      # new_draw_vec_exp <- exp(new_draw_vec)
      # new_draw_vec_probs <- new_draw_vec_exp / sum(new_draw_vec_exp)
      
      
      # sample a new membership value
      
      new_membership <- sample(1:total_clusters, 1, prob = probs)
      thetas[i] <- new_membership
      alphas[i] <- alpha_dict[new_membership]
      betas[i] <- beta_dict[new_membership]
    }
  }
  
  param_df <- data.frame(thetas, alphas, betas)
  
  return(param_df)

}



params_MH <- function(alpha, beta, times){
 
  
  # for each cluster, we want to update its alpha and beta parameters
  # in a MH fashion
  # 
  # we suggest the following proposal:  
  #   1) propose a new beta where we take a normally distributed step in the 
  #      log space to ensure that beta remains positive
  #   2) propose a new alpha that is unif[0, 1] * new beta. This is to ensure
  #      the kernel doesnt explode
  #   3) the likelihood of jumping back is equal. 
  # 
  #
  # --------------------------------------------------------------
  # input:
  #   - alpha (number): alpha parameter for a particular cluster
  #   - beta (number): beta parameter for a particular cluster
  #   - times (vector): times of all the offspring of parents in the particular
  #                     cluster standardized by their parents
  #                     
  # output:
  #   - alpha, beta (vector of length 2): accepted or original alpha and beta 
  #                                       values for MH algorithm 
  # --------------------------------------------------------------
   
  # proposal likelihoods
  
  step <- rnorm(1, 0, 0.1)
  # lik_propose <- dnorm(step)
  propose_beta <- exp(log(beta) + step)
  propose_alpha <- runif(1) * propose_beta
  # 
  # lik_reverse <- dnorm(step)
  
  # data likelihoods
  
  ll_num <- decay_exp_ll(propose_alpha, propose_beta, times)
  
  ll_den <- decay_exp_ll(alpha, beta, times)
  
  # MH part 
  
  log_a_star <- min(0, ll_num - ll_den)
  
  log_u <- log(runif(1))
  
  if(log_u < log_a_star){
    return(c(propose_alpha, propose_beta))
  } else{
    return(c(alpha, beta))
  }
  
}



neal_8 <- function(events,
                   thetas,
                   alphas,
                   betas,
                   mu, tmax, alpha,
                   tree,
                   total_clusters){
  

  # with our data, events
  # and our current values of the latent parameters (thetas, alphas, betas)
  # and our current sampled mu
  # and our current sampled tree
  # 
  # perform neal algorihtm 8 (2000) to
  #   - resample all the memberships
  #   - update all memberships alpha and beta parameters
  # 
  # 
  # ---------------------------------------------------------------------
  # inputs:
  #   
  #   - events (N-dim vector): vector of N observed events timestamps
  #   - thetas (N-dim vector): membership vector for each observation
  #   - alphas (N-dim vector): alpha value for each observation
  #   - betas (N-dim vector): beta value for each observation
  # 
  #   - mu (number): background intensity 
  #   - tmax (number): time limit of observation window
  #   - alpha (number): discount parameter when creating new "tables"
  #   - tree (N-dim vector): each element describes its parent
  #   - total_clusters (number): number of total clusters pre-specified in the DPMM
  # 
  # outputs:
  #   - params (N x 4 dataframe): 
  #     - 1st column = events
  #     - 2nd column = thetas
  #     - 3rd column = alphas
  #     - 4th column = betas
  # ----------------------------------------------------------------

  
  params <- data.frame(events, thetas, alphas, betas)
  # RM <- rate_matrix(events, mu, thetas, alphas, betas)
  # tree <- sample_tree(RM, length(events))
  
  
  # PART 1 - RESAMPLE MEMBERSHIP - UPDATE THETAS, ALPHAS, BETAS
  
  params <- sample_memberships(tree, 
                               events, thetas, alphas, betas,
                               total_clusters, alpha)
  
  thetas <- params$thetas
  alphas <- params$alphas
  betas <- params$betas
  
  params_distinct <- update_phi(thetas, alphas, betas, total_clusters) %>% 
    arrange(thetas)
  
  alpha_dict <- params_distinct$alphas
  
  beta_dict <- params_distinct$betas
  
  theta_dict <- params_distinct$thetas
  
  # PHASE 2 - update params_distinct for clusters with parents
  
  cluster_ID <- unique(thetas)
  
  
  for(i in cluster_ID){ # for cluster ID i
    
    # times of all parents in cluster i
    
    parents <- which(thetas == i)
    
    # times of all children whose parents are in cluster i
    children <- c()
    
    for(j in parents){ # for the jth parent in cluster i
      children_j <- which(tree == j) # obtain their children
      children_times_j <- events[children_j]
      parent_j <- events[j]
      time_j <- children_times_j - parent_j
      children <- c(children, time_j)
    }
    
    # update
    
    param_update <- params_MH(alpha_dict[i], beta_dict[i], children)
      
    # update alpha and beta
    
    params_distinct[i,2] <- param_update[1]
    params_distinct[i,3] <- param_update[2] 
  }
  
  # PHASE 3 - store alphas and betas in the MCMC matrix
  
  distinct_alphas <- params_distinct$alphas
  distinct_betas <- params_distinct$betas
  alphas <- distinct_alphas[thetas]
  betas  <- distinct_betas[thetas]
  
  params <- data.frame(events, thetas, alphas, betas)
  
  return(params)

}


