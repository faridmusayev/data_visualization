

# libraries
library(ggplot2)
library(gridExtra)
library(grid)

# Given 
y <- readRDS("Precipitation.rds") # daily records of weather
y <- log(y + 10^(-30))
n <- length(y)
sample_mean <- mean(y)
ndraws <- n 

# Prior input parameters
mean_0 <- 1         # mean for prior mean
tau_0 <- 1          # variance for prior mean
variance_0 <- 1     # variance for prior variance
v_0 <- 1            # degrees of freedom for prior variance
v_n <- v_0 + n      # degrees of freedom for full conditional posterior variance

# posterior initial parameters
init_posterior_mean <- 1


# full conditional posterior variance simulation
posterior_variance <- function(ndraws = 1, posterior_mean = 1){
  x <- rchisq(n = ndraws, v_n)
  s_2 <- (v_0 * variance_0 + sum((y - posterior_mean)**2))/v_n
  posterior_variance <- v_n * s_2/ x
  return(posterior_variance)
}

# full conditional posterior mean simulation
posterior_mean <- function(ndraws = 1, posterior_variance){
  tau_n <- (1/(n/posterior_variance + 1/tau_0**2))**0.5
  w <- (n/posterior_variance) / ((n/posterior_variance) + (1/tau_0**2))
  mean_n <- w * sample_mean + (1 - w) * mean_0
  posterior_mean <- rnorm(1, mean_n, tau_n**2)
  return(posterior_mean)
}

# Gibbs sampling process
sim_vars <- c()
sim_means <- c()
simulated_mean <- 1 # which is initial posterior mean
for (i in 1:ndraws){
  simulated_variance <- posterior_variance(posterior_mean = simulated_mean)
  simulated_mean <- posterior_mean(posterior_variance = simulated_variance)
  sim_vars[i] <- simulated_variance
  sim_means[i] <- simulated_mean
}


# Graphical Representation
df <- data.frame(sim_vars = sim_vars,
                 sim_means = sim_means, 
                 avg_sim_vars = cumsum(sim_vars)/seq(1, ndraws), 
                 avg_sim_means = cumsum(sim_means)/seq(1, ndraws),
                 ndraws = 1:ndraws)

pl <- ggplot(data = df) + 
  geom_line(aes(x = ndraws, y = sim_vars), size = 0.4, color = "orange") + 
  labs(title = "Simulated Variances at each iteration", 
       x = "Draw #", y = "Simulated variance") +
  theme_minimal() +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = "plain"))


pl2 <- ggplot(data = df) + 
  geom_line(aes(x = ndraws, y = sim_means), size = 0.4, color = "darkblue") + 
  labs(title = "Simulated Means at each iteration", 
       x = "Draw #", y = "Simulated mean") +
  theme_minimal() +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = "plain"))

pl3 <- ggplot(data = df) + 
  geom_line(aes(x = ndraws, y = avg_sim_vars), size = 0.4, color = "orange") + 
  labs(title = "Average of simulated variances at each iteration", 
       x = "Draw #", y = "Simulated variance") +
  theme_minimal() +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = "plain"))

pl4 <- ggplot(data = df) + 
  geom_line(aes(x = ndraws, y = avg_sim_means), size = 0.4, color = "darkblue") + 
  labs(title = "Average of simulated means at each iteration", 
       x = "Draw #", y = "Simulated mean") +
  theme_minimal() +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = "plain"))


grid.arrange(grobs = list(pl, pl2, pl3, pl4), 
             top = textGrob("Gibbs Sampling Convergence",
                            gp = gpar(fontsize = 24)))



g <- arrangeGrob(pl, pl2, pl3, pl4) 

ggsave('gibbs_sampling.svg', g, width = 12, height = 8)
