
# libraries
library(here)
library(ggplot2)
library(gridExtra)
library(mvtnorm)


# Maximum likelihood estimator of Beta in the Poisson regression model

df <- read.table("eBayNumberOfBidderData.dat", header = TRUE)

df_woc <- df[, -2] # data frame without Const

poisson_reg <- glm(formula = nBids ~ ., data = df_woc, family = poisson(link = "log"))

print(summary(poisson_reg))
# Essential covariates are VerifyID, Sealed, LogBook, MinBidShare


X <- as.matrix(df[, -1])
y <- as.matrix(df[, 1])
features_dim <- dim(X)[2]

# prior mean 'mu' and prior covariance matrix 'prior_covmat'
prior_covmat <- 100 * solve(t(X)%*%X)
mu <- rep(0, features_dim) 

# initial Beta parameters
betas_init <- rep(0,features_dim)

logposterior <- function(betas, X, y, mu, sigma){
  
  # transform matrix beta into single column vector
  betas = as.vector(betas)
  
  # log-likelihood for poisson regression model
  loglik <- sum((X %*% betas) * y) - sum(exp(X %*% betas))
  
  # log-normal prior
  logprior <- dmvnorm(betas, mean = mu, sigma = sigma, log = TRUE)
  
  return(loglik + logprior)
}

# perform optimization
optim_result <- optim(par = betas_init,  # initial values for parameters
                      fn = logposterior, # function to be minimized
                      gr = NULL, 
                      X, y, mu, prior_covmat, # arguments of function to be minimized
                      method = c("BFGS"), 
                      control = list(fnscale = -1), # minimize negative logposterior
                      hessian = TRUE)

# calculated posterior mode aka Beta parameters
posterior_betas <- optim_result$par
names(posterior_betas) <- colnames(X)

# calculated posterior covariance matrix by using hessian matrix from optim function
posterior_covmat <- -solve(optim_result$hessian)


# Function for Random Walk Metropolis Algorithm
random_walk_metropolis <- function(c, logposterior){ 
  
  c = c # tuning parameter
  N <- 1000 # number of iterations for algorithm
  
  # matrix for storing all accepted draws 
  thetas <- matrix(0, nrow = N, ncol = features_dim) # accepted draws are stored here
  # set counter for rejection rate
  rejection_rate <- 0  
  # store alpha values for comparison and evaluation of average acceptance probability
  alphas <- c()
  
  # run Random Walk Metropolis algorithm in for loop for N draws
  for (i in 2:N){
    
    # simulate sample proposal from proposed distribution
    sample_proposal <- rmvnorm(1, mean = thetas[i - 1, ], sigma = c * posterior_covmat)
    
    # evaluate probability for sample proposal
    p_sample_proposal <- logposterior(betas = sample_proposal, X, y, mu, prior_covmat)
    
    # evaluate probability for the last accepted draw
    p_last_accepted_draw <- logposterior(betas = thetas[i - 1, ], X, y, mu, prior_covmat)
    
    # calculate acceptance probability
    alpha <- min(1, exp(p_sample_proposal - p_last_accepted_draw)) # use exp because of logs
    alphas[i] <- round(alpha, 2) # this is just for myself, no need to store in reality
    
    # generate random number from uniform distribution
    U <- runif(1)
    # compare it with acceptance probability and if lower
    if (U < alpha){
      # add sample proposal to matrix of accepted draws
      thetas[i, ] <- sample_proposal
      
    }else{
      # if not accepted add to current i iteration of thetas a previous value i - 1
      thetas[i, ] <- thetas[i - 1, ]
      rejection_rate <- rejection_rate + 1
    }
    
  } 
  cat("Average acceptance rate = ", 
      1 - rejection_rate/N) 
  
  # or mean(alphas) where alpha is acceptance probability
  print(mean(alphas, na.rm = TRUE))
  # return accepted draws
  return(thetas)
}

N <- 1000
thetas <- random_walk_metropolis(c = 1, logposterior = logposterior)

# Graphical representation of convergence
results <- cbind(iterations = c(1:N), as.data.frame(thetas))
colnames(results)[-1] <- colnames(X)

pl <- ggplot(data = results) +
  geom_line(aes(x = iterations, y = VerifyID, color = "P1")) +
  geom_line(aes(x = iterations, y = Sealed, color = "P2")) +
  geom_line(aes(x = iterations, y = LogBook, color = "P3")) +
  geom_line(aes(x = iterations, y = MinBidShare, color = "P4")) +
  scale_color_manual(name = 'Parameters:', 
                     values = c(P1 = "darkblue", P2 = "orange",
                                P3 = "darkgrey", P4 = "black")) +
  labs(x = "Iteration #", y = "Parameter Value") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.text = element_text(size = 12),
        legend.title = element_text(face = 'bold'))

plot(pl)  

ggsave('mcmc.svg', width = 8, height = 6)



