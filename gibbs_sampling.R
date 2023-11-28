
##w  <- 1.999
##xv <- seq(-1, 1, by=0.01) * 1/sqrt(1-w^2/4)  # a range of x1-values, where the term below the root is non-negative (compare Lecture 4)
##plot(xv, xv, type="n", xlab=expression(x[1]), ylab=expression(x[2]), las=1)
### ellipse
##lines(xv, -(w/2)*xv-sqrt(1-(1-w^2/4)*xv^2), lwd=2, col=10)
##lines(xv, -(w/2)*xv+sqrt(1-(1-w^2/4)*xv^2), lwd=2, col=10)

# load necessary library
library(ggplot2)

# Function for the target bivariate distribution
target_bivariate_function <- function(x1, x2, w) {
  indicator <- (x2^2 + w*x1*x2 + x1^2 < 1)
  return(indicator)
}

w <- 1.999
x1 <- seq(-1, 1, length.out = 200) * 1/sqrt(1-w^2/4) # a range of x1-values, where the term below the root is non-negative (compare Lecture 4)
x2 <- seq(-1, 1, length.out = 200) * 1/sqrt(1-w^2/4) # a range of x1-values, where the term below the root is non-negative (compare Lecture 4)
df <- expand.grid(x1 = x1, x2 = x2)
df$inside <- target_bivariate_function(df$x1, df$x2, w)
ggplot(df, aes(x1, x2)) +
  geom_tile(aes(fill = inside), alpha = 0.3) +
  ggtitle("The boundaries of the region where X has a uniform distribution") + coord_fixed()


# Gibbs Sampling Function
gibbs_sampling <- function(n, w) {
  samples <- matrix(0, n, 2)
  x1 <- 0
  x2 <- 0
  set.seed(123)

  for (i in 1:n) {
    # Sample x1 given x2
    x1 <- runif(1, -(w/2)*x2-sqrt(1-(1-w^2/4)*x2^2), -(w/2)*x2+sqrt(1-(1-w^2/4)*x2^2))
    # Sample x2 given x1
    x2 <- runif(1, -(w/2)*x1-sqrt(1-(1-w^2/4)*x1^2), -(w/2)*x1+sqrt(1-(1-w^2/4)*x1^2))
    samples[i,] <- c(x1, x2)
  }
  return(samples)
}

# Set parameters
n <- 1000
w <- 1.999

# Run Gibbs sampling
samples <- gibbs_sampling(n, w)

# Plot the boundaries
ggplot(df, aes(x1, x2)) +
  geom_tile(aes(fill = inside), alpha = 0.3) +
  geom_point(data = as.data.frame(samples), aes(x = V1, y = V2), color = "red", alpha = 0.5) +
  ggtitle("Gibbs Sampling for X") + coord_fixed()


# Determine P(X1 > 0)
prob_x1_positive <- mean(samples[,1] > 0)
#prob_x1_positive <- sum(samples[,1] > 0) / nrow(samples)
cat("P(X1 > 0):", prob_x1_positive, "\n")


# U Boundaries
u_boundary_function <- function(u1, u2, w) {
  u_indicator <- ((2+w)*u2^2 + (2-w)*u1^2 < 4)
  return(u_indicator)
}

w <- 1.999
u1 <- seq(-2, 2, length.out = 200) * sqrt(4/(2-w))  # a range of u1-values
u2 <- seq(-2, 2, length.out = 200) * sqrt(4/(2+w))  # a range of 21-values
u_df <- expand.grid(u1 = u1, u2 = u2)
u_df$density <- u_boundary_function(u_df$u1, u_df$u2, w)
ggplot(u_df, aes(u1, u2)) +
  geom_tile(aes(fill = density), alpha = 0.3) + xlim(-70,70) + ylim(-2,2) +
  ggtitle("The boundaries of U")

# Gibbs Sampling function for U
u_function <- function(n, w) {
  set.seed(123)
  u <- matrix(nrow=n, ncol=2)
  u[1,] <- c(0, 0)  # starting value

  for (i in 2:n) {
    u[i,1] <- runif(1, -sqrt((4-(2+w)*u[i-1,2]^2)/(2-w)), sqrt((4-(2+w)*u[i-1,2]^2)/(2-w)))
    u[i,2] <- runif(1, -sqrt((4-(2-w)*u[i,1]^2)/(2+w)), sqrt((4-(2-w)*u[i,1]^2)/(2+w)))
  }
  return(u)
}

# Set parameters
n <- 1000
w <- 1.999

# Run u function
u_samples <- u_function(n, w)

# Plot gibbs sampling for U
ggplot(u_df, aes(u1, u2)) +
  geom_tile(aes(fill = density), alpha = 0.3) + xlim(-70,70) + ylim(-2,2) +
  geom_point(data = as.data.frame(u_samples), aes(x = V1, y = V2), color = "blue", alpha = 0.5) +
  ggtitle("Gibbs Sampling for U")


# Determine P(X1 > 0) = P((U2+U1)/2 > 0)
prob_x1_positive_u <- mean((u_samples[,2] + u_samples[,1])/2 > 0)
cat("P(X1 > 0) = P((U2+U1)/2 > 0):", prob_x1_positive_u, "\n")

