
# bob_deploy.R
#
# Author: Bob Verity
# Date: 2022-10-06
#
# Purpose:
# Test package functions.

# RStudio shortcuts:
#    cmd+shift+L     : load package from local version
#    cmd+shift+D     : document (NB, option must be activated in Build Tools)
#    cmd+shift+E     : check
#    cmd+shift+T     : test

# Useful commands:
# devtools::install()  # install package
# pkgdown::build_site() # build all pages of pkgdown website
# pkgdown::build_article('blocks')  # build single vignette
# check('.', args = '--no-examples')  # run checks without examples
# covr::report()    # interactive coverage report
# devtools::build_vignettes()
# ------------------------------------------------------------------

#library(Rcpp)
library(dplyr)
library(ggplot2)

set.seed(1)

# ------------------------------------------------------------------

# make up some data
L <- 1e3
z <- sample(c(0.3, 0.7), L, replace = TRUE)
wsaf <- rbeta(L, shape1 = z*30, shape2 = (1 - z)*30)
a <- round(wsaf*100)
r <- round((1 - wsaf)*100)
p <- rep(0.5, L)

# plot what the WSAF looks like
plot(a / (a + r), ylim = c(0, 1))

# set vector of thermodynamic powers
beta <- seq(0, 1, 0.05)^3
#beta <- 1

# run MCMC
my_mcmc <- run_mcmc(a = a, r = r, p = p,
                    burnin = 1e2,
                    samples = 1e3,
                    beta = beta)

plot(my_mcmc$diagnostics$MC_accept_burnin, ylim = c(0, 1))
plot(my_mcmc$diagnostics$MC_accept_sampling, ylim = c(0, 1))

# have a look at output
head(my_mcmc$draws)
tail(my_mcmc$draws)

# plot posterior mu draws
my_mcmc$draws %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = mu_1), size = 0.5) +
  geom_point(aes(x = iteration, y = mu_2), size = 0.5, col = "red") +
  ylim(c(0, 1)) +
  ylab("mu")

# plot posterior sigma draws
my_mcmc$draws %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = sigma))

# plot posterior w draws
my_mcmc$draws %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = w)) +
  ylim(c(0, 1))
