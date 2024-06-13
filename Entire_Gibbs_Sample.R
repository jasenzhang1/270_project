library(tidyverse)
library(gridExtra)
library(coda)
library(reshape2)

source("helper.R")
source("Sample_mu.R")
source("generate_data.R")

set.seed(270)

# simulation parameters

tmax <- 200
a1 <- 5
b1 <- 10
mu1 <- 0.05

a2 <- 1
b2 <- 4
mu2 <- 0.1

# simulated data

events1 <- immigrant_children_simulation(tmax, a1, b1, mu1)
events2 <- immigrant_children_simulation(tmax, a2, b2, mu2)

# save simulated data

g_events1 <- make_stepfunction(events1)
g_events2 <- make_stepfunction(events2)

g_all_events <- grid.arrange(g_events1, g_events2, nrow = 2, ncol = 1)

ggsave(file="sim_data.pdf", width = 8, height = 6, g_all_events)



events <- c(events1, events2)
events <- sort(events)

n_events <- length(events)


# assumed fixed constants
alpha <- 1
total_clusters <- 5

# starting values
thetas  <- rep(c(1,2),500)
thetas <- thetas[1:n_events]

alphas <- rep(c(1.5, 2), 500)
alphas <- alphas[1:n_events]

betas <- rep(c(2, 5), 500)
betas <- betas[1:n_events]

mu <- 0.2

# MCMC parameters
niter <- 1000
burnin <- 100

theta_MCMC <- matrix(0, niter, total_clusters)
alpha_MCMC <- matrix(0, niter, total_clusters)
beta_MCMC <- matrix(0, niter, total_clusters)
mu_MCMC <- matrix(0, niter, 1)



for(i in 1:niter){

  # STEP 1 - SAMPLE THE TREE
  
  tree <- sample_tree(events, mu, thetas, alphas, betas)

  # STEP 2 - SAMPLE MU
  
  mu <- sample_mu(tree, tmax)
  
  # STEP 3 - SAMPLE NEW MEMBERSHIPS
  
  
  params <- neal_8(events,
                   thetas,
                   alphas,
                   betas,
                   mu, tmax, alpha,
                   tree,
                   total_clusters)
  
  
  thetas <- params$thetas
  alphas <- params$alphas
  betas  <- params$betas
  
  # STEP 3.5 - SORT THEM
  
  params_distinct <- update_phi(thetas, alphas, betas, total_clusters) %>% 
    arrange(thetas)
  
  alpha_dict_store <- params_distinct$alphas
  beta_dict_store  <- params_distinct$betas
  theta_dict_store <- params_distinct$thetas
  
  alpha_dict_store[! (1:total_clusters %in% unique(thetas))] <- NA
  beta_dict_store[! (1:total_clusters %in% unique(thetas))]  <- NA
  theta_dict_store[! (1:total_clusters %in% unique(thetas))]  <- NA
  theta_dict_store[! is.na(theta_dict_store) ] <- 1
  # STEP 4 - STORE VALUES
  
  mu_MCMC[i,1] <- mu
  theta_MCMC[i,] <- theta_dict_store
  alpha_MCMC[i,] <- alpha_dict_store
  beta_MCMC[i,] <- beta_dict_store
  
  if(i %% 100 == 0){
    print(i)
  }
  
}

## Burnin

theta_MCMC <- theta_MCMC[(burnin+1):niter,]
alpha_MCMC <- alpha_MCMC[(burnin+1):niter,]
beta_MCMC <- beta_MCMC[(burnin+1):niter,]
mu_MCMC <- mu_MCMC[(burnin+1):niter,]

apply(theta_MCMC, 2,function(x){sum(x, na.rm = T)}) / 
  rep((niter - burnin ), total_clusters)

g_mu <- ggplot() + geom_line(aes(x = (burnin+1):niter, y = mu_MCMC))
g_mu

g_alpha <- ggplot() + geom_line(aes(x = (burnin+1):niter, y = log(alpha_MCMC[,1])))
g_alpha

g_beta <- ggplot() + geom_line(aes(x = (burnin+1):niter, y = log(beta_MCMC[,1])))
g_beta

## ACF

# acf(mu_MCMC)
# acf(alpha_MCMC[,1])
# 
# mcmc.obj1 <- as.mcmc(mu_MCMC)
# mcmc.obj2 <- as.mcmc(alpha_MCMC)
# mcmc.obj3 <- as.mcmc(beta_MCMC)
# 
# ess_mu <- mean(coda::effectiveSize(mcmc.obj1)) 
# ess_alpha <- mean(coda::effectiveSize(mcmc.obj2)) 
# ess_beta <- mean(coda::effectiveSize(mcmc.obj3)) 
# 
# print(ess_mu)
# print(ess_alpha)
# print(ess_beta)

## histograms

g_mu_hist <- ggplot() + geom_density(aes(x = mu_MCMC)) + 
  geom_vline(xintercept = 0.15, linetype = 'dashed', alpha = 0.4, color = 'black') 
g_mu_hist

alpha_df <- data.frame(alpha_MCMC)
alpha_df$ID <- (burnin+1):niter 
alpha_df_tall <- melt(alpha_df, id.vars = "ID")

beta_df <- data.frame(beta_MCMC)
beta_df$ID <- (burnin+1):niter 
beta_df_tall <- melt(beta_df, id.vars = "ID")

g_alpha_total <- ggplot(data = alpha_df_tall) + 
  geom_density(aes(x = value, color = variable)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', alpha = 0.4, color = 'red') + 
  geom_vline(xintercept = 5, linetype = 'dashed', alpha = 0.4, color = 'blue') + 
  xlim(0, 15) + 
  labs(color = "alpha")
g_alpha_total

g_beta_total <- ggplot(data = beta_df_tall) + 
  geom_density(aes(x = value, color = variable)) + 
  geom_vline(xintercept = 4, linetype = 'dashed', alpha = 0.4, color = 'red') + 
  geom_vline(xintercept = 10, linetype = 'dashed', alpha = 0.4, color = 'blue') +   
  xlim(0,15) + 
  labs(color = "beta")
g_beta_total

g_all_clusters <- grid.arrange(g_mu_hist, g_alpha_total, g_beta_total, nrow = 3, ncol = 1)


ggsave(file="alpha_beta.pdf", width = 8, height = 9, g_all_clusters)

alpha_beta_graphs <- list()

for(i in 1:total_clusters){
  g <- ggplot() + 
    geom_density(aes(x = alpha_MCMC[,i], color = 'alpha')) + 
    geom_density(aes(x = beta_MCMC[,i], color = 'beta'))
  
  index <- i
  
  g
  alpha_beta_graphs[[index]] <- g
}


grid.arrange(alpha_beta_graphs[[1]],
             alpha_beta_graphs[[2]],
             alpha_beta_graphs[[3]],
             alpha_beta_graphs[[4]],
             alpha_beta_graphs[[5]], nrow = 2, ncol = 3)

 
apply(alpha_MCMC, 2, function(x){mean(x, na.rm = T)})
apply(beta_MCMC, 2, function(x){mean(x, na.rm = T)})

