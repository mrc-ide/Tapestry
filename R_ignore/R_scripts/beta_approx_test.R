
# Beta_approx_test.R
#
# Author: Bob Verity
# Date: 2022-10-14
#
# Purpose:
# Test Beta approximation to Beta-binomial from DEploidIBD paper
# Conclusion: it's pretty good, although gets a bit janky for very low allele
# counts
#
# ------------------------------------------------------------------

# define counts
a <- 5
r <- 10
n <- a + r

# range of true WSAF to explore
q <- seq(0, 1, l = 101)

# variance parameter
c <- 100

# get true Beta-binomial distribution
alpha <- q*c
beta <- (1 - q)*c
y1 <- extraDistr::dbbinom(a, n, alpha = alpha, beta = beta)

# manually compute mean and variance of distribution
m <- sum(q*y1) / sum(y1)
m2 <- sum(q^2*y1) / sum(y1)
v <- m2 - m^2

# get Beta parameters from mean and variance, i.e. match on first two moments
u <- m^2*(1 - m) / v - m
v <- u / m - u

# get fitted Beta distribution
y2 <- dbeta(q, u, v)

# plot and compare
plot(q, y1 / sum(y1))
lines(q, y2 / sum(y2))



