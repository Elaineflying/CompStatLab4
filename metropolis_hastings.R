# Function for the target distribution
target_distribution <- function(x) {
  return(x^5 * exp(-x))
}

# Metropolis hastings method when proposal distribution is log-normal function
metropolis_hastings_log_normal <- function(n) {
  set.seed(123)
  chain <- rep(0,n)
  chain[1] <- 1
  acceptance_rate <- 0

  for (i in 2:n) {
    proposal <- rlnorm(1, meanlog = log(chain[i-1]), sdlog = 1)
    ratio <- (target_distribution(proposal) * dlnorm(chain[i-1], meanlog = log(proposal), sdlog = 1)) / (target_distribution(chain[i-1]) * dlnorm(proposal, meanlog = log(chain[i-1]), sdlog = 1))
    if (runif(1) < min(1,ratio)) {
      chain[i] <- proposal
      acceptance_rate <- acceptance_rate + 1
    } else {
      chain[i] <- chain[i-1]
    }
  }
  acceptance_rate <- acceptance_rate / n
  return(list(chain=chain,acceptance_rate=acceptance_rate))
}

samples1 <- metropolis_hastings_log_normal(n=10000)
cat("The acceptance rate is :", samples1$acceptance_rate,"\n")

plot(samples1$chain, type = "l", main = "Metropolis-Hastings Chain (log-normal)", xlab = "Iterations", ylab = "Sample Value")
hist(samples1$chain, main = "Histogram of the sample (log-normal)", xlab = "Sample Value", breaks = 50,  col = "lightblue")


samples1$chain[20]

# Burn-in period analysis
burn_in_period <- 11
burn_in_chain <- samples1$chain[(burn_in_period + 1):10000]

# Plot the trace plot after burn-in
plot(burn_in_chain, type = 'l', xlab = 'Iterations', ylab = 'Sample Value', main = 'Metropolis-Hastings Chain (After Burn-In)')

# Autocorrelation plot
acf(burn_in_chain, lag.max = 50, main = 'Autocorrelation Plot (After Burn-In)')

# Histogram of the sample after burn-in
hist(burn_in_chain, main = 'Histogram of Samples (After Burn-In)', xlab = 'Sample Value', col = 'lightblue')


# Metropolis hastings method when proposal distribution is chi-square function
metropolis_hastings_chi_square <- function(n) {
  set.seed(123)
  chain <- rep(0,n)
  chain[1] <- 1
  acceptance_rate <- 0

  for (i in 2:n) {
    proposal <- rchisq(1, df = floor(chain[i-1]) + 1)
    ratio <- (target_distribution(proposal) * dchisq(chain[i-1], df = floor(proposal) + 1)) / (target_distribution(chain[i-1]) * dchisq(proposal, df = floor(chain[i-1]) + 1))
    if (runif(1) < min(1,ratio)) {
      chain[i] <- proposal
      acceptance_rate <- acceptance_rate + 1
    } else {
      chain[i] <- chain[i-1]
    }
  }
  acceptance_rate <- acceptance_rate / n
  return(list(chain=chain,acceptance_rate=acceptance_rate))
}

samples2 <- metropolis_hastings_chi_square(n=10000)
cat("When the proposal function is chi-squre function, the acceptance rate is :", samples2$acceptance_rate,"\n")

plot(samples2$chain, type = "l", main = "Metropolis-Hastings Chain (chi-square)", xlab = "Iterations", ylab = "Sample Value")
hist(samples2$chain, main = "Histogram of the sample (chi-suare)", xlab = "Sample Value", breaks = 50,  col = "lightblue")


# Metropolis hastings method when proposal distribution is normal function
metropolis_hastings_normal <- function(n) {
  set.seed(123)
  chain <- rep(0,n)
  chain[1] <- 1
  acceptance_rate <- 0

  for (i in 2:n) {
    proposal <- rnorm(1, mean = chain[i-1], sd = 1)
    # proposal normal function is symmetric, therefore no need to calculate g(x(t)|x(*)) and g(x(*)|x(t))
    ratio <- target_distribution(proposal) / target_distribution(chain[i-1])
    if (runif(1) < min(1,ratio)) {
      chain[i] <- proposal
      acceptance_rate <- acceptance_rate + 1
    } else {
      chain[i] <- chain[i-1]
    }
  }
  acceptance_rate <- acceptance_rate / n
  return(list(chain=chain,acceptance_rate=acceptance_rate))
}

samples3 <- metropolis_hastings_normal(n=10000)
cat("When the proposal function is normal function, the acceptance rate is :", samples3$acceptance_rate,"\n")

plot(samples3$chain, type = "l", main = "Metropolis-Hastings Chain (normal)", xlab = "Iterations", ylab = "Sample Value")
hist(samples3$chain, main = "Histogram of the sample (normal)", xlab = "Sample Value", breaks = 50,  col = "lightblue")

results_df <- data.frame(proposal_function=c("log-normal", "chi-square", "normal"),
                         acceptance_rate=c(samples1$acceptance_rate, samples2$acceptance_rate, samples3$acceptance_rate),
                         mean=c(mean(samples1$chain), mean(samples2$chain), mean(samples3$chain)))

print(results_df)




# Define the probability density function of the gamma distribution
gamma_distribution <- function(x) {
  constant <- 6  # Given constant of proportionality
  return(constant * x^5 * exp(-x))
}

# Use the integrate function to compute the integral
result <- integrate(target_distribution, lower = 0, upper = Inf)

# Display the result
cat("Actual value of the integral:", result$value, "\n")
