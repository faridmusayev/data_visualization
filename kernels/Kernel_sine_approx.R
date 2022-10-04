
# libraries
library(ggplot2)

# Kernel sin() function approximation
# Using different "h" width parameters

# Code for sampling the training data.
n <- 500
x <- runif(n, 0, 20)
y <- sin(x)
set.seed(12345)

# Gaussian Kernel function
k <- function(value, h){
  # if the distance is the smallest, then kernel values are the highest
  gaussian_kernel <- exp(-(value)**2/(2*h**2))
  return(gaussian_kernel)
}

# loop through different h distances 
h = c(0.1, 1, 10)

# data frame for storing all y results for each h distance
results <- data.frame(h_01 = numeric(n), 
                  h_1 = numeric(n), 
                  h_10 = numeric(n))

y_preds <- c()
for (i in 1: length(h)){
  
  for (j in 1:n){
    diff = x - x[j]
    kernel_weights <- k(value = diff, h = h[i])
    y_preds[j] <- sum(kernel_weights * y)/sum(kernel_weights)
  } 
  results[,i] <- y_preds
}


df <- data.frame(x, y, results)

# MSE values for y predictions with different h
mse_h_01 <- mean((df$y - df$h_01)**2)
cat("h = 0.1,  MSE = ", mse_h_01, "\n") 

mse_h_1 <- mean((df$y - df$h_1)**2)
cat("h = 1,  MSE = ", mse_h_1, "\n") 

mse_h_10 <- mean((df$y - df$h_10)**2)
cat("h = 10,  MSE = ", mse_h_10, "\n") 

# Plot using h with smallest MSE
pl <- ggplot(data = df) +
      geom_point(aes(x = x , y = y, color = "1"), size = 1.2) +
      geom_point(aes(x = x, y = h_01, color = "2"), size = 1.2) +
      scale_color_manual(name = NULL,
                         labels = c("sin(x)", "Approximated with Gaussian Kernel"), 
                         values = c("darkblue","orange")) +
      labs(x = 'x', y = 'sin(x)') + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.text = element_text(size = 12)) +
      guides(colour = guide_legend(override.aes = list(size = 2)))


plot(pl)


ggsave('kernel_sin.svg', pl, width = 8, height = 4)




