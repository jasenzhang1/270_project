
immigrant_children_simulation <- function(tmax, a, b, mu){
  
  # simulate a hawkes process with an exponential trigger function
  # with the immigrant children process
  # 
  # input:
  #   - tmax (number): time limit of the simulation
  #   - a (number): alpha parameter in the exponential trigger function
  #   - b (number): beta parameter in the exponential trigger function
  #   - mu (number): baseline intensity value
  # 
  # output:
  #   - events (vector of numbers): timepoint of the simulated events
  
  t <- 0
  
  # first simulate the immigrants
  backgroundEvents <- vector()
  while(t < tmax){
    
    t <- t + rexp(n=1,rate=mu)
    
    if(t < tmax){
      backgroundEvents <- c(backgroundEvents, t)
    }
    
  }
  
  events     <- backgroundEvents
  currentGen <- backgroundEvents
  
  while(length(currentGen) > 0) {
    
    nextGen <- vector()
    for(k in 1:length(currentGen)) {
      
      # the number of children is pois(lambda = a/b)
      nChilds <- rpois(n=1,lambda=a/b)
      
      # if there is at least one child, simulate the time the event
      # occurs and store the time in the nextGen vector
      
      if (nChilds > 0) {
        childs <- currentGen[k] + rexp(n = nChilds, rate = b)
        nextGen <- c(nextGen,childs)
      }
      
    }
    
    # once all the events in currentGen have been simulated for children
    # store nextGen as currentGen
    
    currentGen <- nextGen
    events <- c(events,currentGen)
    
  }
  
  # when there are no further children, sort and cut off the children
  # whose times are past tmax
  
  events <- sort(events)
  events <- events[events < tmax]
  
  return(events)
}

make_stepfunction <- function(events, tmax){
  
  # create a stepfunction graph from a dataframe of events over time
  # 
  # input:
  # - events (vector): time at which each event occurs
  # - tmax (number): tmax
  # 
  # output:
  # - g (graph): a graph showing the stepfunction described by the dataframe
  
  N <- length(events)
  
  df <- data.frame(events, 1:N)
  
  
  colnames(df) <- c('x', 'cdf')
  new_steps <- tibble(df)
  
  
  df <- new_steps %>%
    mutate(type = "cdf") %>%
    bind_rows(new_steps %>%
                mutate(type = "prior",
                       cdf = lag(cdf))) %>%
    drop_na() %>%
    arrange(x, desc(type))
  # 
  #df <- rbind(c(0,0,'prior'), df)
  #df <- rbind(c(0,0,'cdf'), df)
  
  g <- ggplot(df) + 
    geom_point(aes(x, cdf, fill = type),
               shape = 21) +
    scale_fill_manual(values = c("black", "white")) +
    geom_segment(aes(x = lag(x), y = lag(cdf),
                     xend = x, yend = cdf,
                     lty = type)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    xlim(0, tmax) +
    theme(legend.position="none")
  
  return(g)
}
