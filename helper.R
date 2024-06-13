decay_exp <- function(alpha, beta, t){

  #
  # caluculate the density of observing a time for a exponential decay kernel
  #
  # input:
  #     - alpha (number): alpha parameter
  #     - beta (number): beta parameter
  #     - t (number or vector): observed time(s) of event(s)
  #   
  # output:
  #     - likelihood (number or vectors): likelihood of observing each time
  
  likelihood <- alpha * exp(- beta * t)
  
  return(likelihood)
}





decay_exp_ll <- function(alpha, beta, times){

  # with given alpha and beta parameters, we find the log likelihood of 
  # observing children at the times

  # --------------------------------------------------------------------  
  # input:
  #     - alpha (number): alpha parameter
  #     - beta (number): beta parameter
  #     - times (vector): vector of observed times
  #   
  # output:
  #     - ll (number) log likelihood of observing the vector of times
  #
  #----------------------------------------------------------------------

  N <- length(times)
  
  ll <- N * log(alpha) - beta * sum(times)
  
  return(ll)
}





rate_matrix <- function(events, mu, theta, alpha, beta){
  
  # calculate the matrix that describes the likelihood of the ith event being 
  # the parent of the jth event.
  # 
  # An event can only be the parent of another event if it came before it. 
  #
  # ------------------------------------------------------------------
  # input:
  #   - events (N-dim vector): vector of N observed events timestamps
  #   - mu (number): background intensity value
  #   - theta (N-dim vector): membership vector for each observation
  #   - alpha (N-dim vector): alpha value for each observation
  #   - beta (N-dim vector): beta value for each observation
  # 
  # output:
  #   - RM5 ((N+1) x N matrix): [i,j]th entry specifies the (cumulative)
  #                             probablity that event
  #                             (i-1) is the parent of event j.
  # 
  #                             EX: [1,4]th entry = event 0 is 
  #                                 the parent of event 4
  # ---------------------------------------------------------------------
                            
                    
  N <- length(events)
  
  RM <- matrix(0, N+1, N+1)
  
  # fill out background first
  
  for(i in 2:(N+1)){
    RM[1,i] <- mu
  }
  
  # fill out the rest
  
  # for the jth row and ith column
  # jth event causes ith event
  
  for(j in 2:N){
    for(i in (j+1):(N+1)){ 
      time_diff <- events[i-1] - events[j-1]
      
      # extract membership parameter and cluster parameters
      theta_i <- theta[j-1]
      alpha_i <- alpha[j-1]
      beta_i <- beta[j-1]
      RM[j,i] <- decay_exp(alpha_i, beta_i, time_diff)
    }
  }
  
  # delete first column because 0th event can't have parents
  
  RM <- RM[,2:(N+1)]
  
  # normalize each column 
  col_sums <- colSums(RM)
  RM2 <- sweep(RM, 2, col_sums, FUN = "/")
  
  # flip the matrix vertically
  RM3 <- RM2[(N+1):1,]
  
  # cumsum vertically
  RM4 <- apply(RM3, 2, cumsum)
  
  # flip the matrix vertically again
  # this way, the matrix seems upper-diagonal with values that "trickle up"
  # and add to 1 at the very top
  
  RM5 <- RM4[(N+1):1,]
  
  return(RM5)
}




sample_parent <- function(RM, i){

# Sampling a parent for the ith event
# i = the event we want to sample a parent for.
# 
# Possible candidates are events 1 to (i-1) or the background intensity.
# 
# 
# --------------------------------------------------------------------------
# input:
#   - RM ((N+1) x N matrix):  rate matrix we created with probabilities
#                             that trickle up and sum to 1
#   - i (number): event ID corresponding to the order in which the event
#                 occurred. ID = 4 means its the 4th event that occurred
# output:
#   - parent (number): sampled parent of the ith event (0 to i-1 are possible)
# ---------------------------------------------------------------------------


  
  
  RM_slice = RM[,i]
  u <- runif(1)
  
  parent <- tail(which(RM_slice > u), 1) - 1
  
  return(parent)
}





sample_tree <- function(events, mu, thetas, alphas, betas){
  
# -----------------------------------------------------------------------
# Now, sample an entire tree
# 
# input:
#   - events (N-dim vector): vector of N observed events timestamps
#   - mu (number): background intensity value
#   - theta (N-dim vector): membership vector for each observation
#   - alpha (N-dim vector): alpha value for each observation
#   - beta (N-dim vector): beta value for each observation
# 
# output:
#   - parents (N-dim vector): vector of parents for each event
# ------------------------------------------------------------------------
  
  
  N <- length(events)
  
  RM <- rate_matrix(events, mu, thetas, alphas, betas)
  
  parents <- rep(0, N)
  
  for(i in 1:N){
    parents[i] <- sample_parent(RM, i)
  }
  
  return(parents)
}

hawkes_likelihood_recursive <- function(times, a, b, mu){
  
# with values of a, b, and mu, calculate the likelihood of a hawkes process with
# a single exponential trigger function
# 
# ----------------------------------------------------------------------------
# input:
#   - times (vector): vector of events
#   - a (number): alpha parameter in an exponential trigger kernel
#   - b (number): beta parameter in an exponential trigger kernel
#   - mu (number): background intensity
# 
# output:
#   - logl (number): log likelihood of the hawkes process given the data (times)
# ------------------------------------------------------------------------------
  
  k <- length(times)
  
  
  # get A(i)'s 
  vec_a <- rep(0,k)
  
  for(i in 2:k){
    vec_a[i] <- exp(-b * (times[i] - times[i-1])) * (1 + vec_a[i-1])
  }
  
  # final observation
  tk <- times[k]
  
  logl <- - mu * tk
  
  for(i in 1:k){
    
    
    term1 <- log(mu + a * vec_a[i])
    
    
    term2 <- (a/b) * (exp(-b * (tk - times[i])) - 1)
    
    logl <- logl + term1 + term2
    
  }
  
  return(logl)
}
