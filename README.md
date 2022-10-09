 <!-- badges: start -->
  [![checks_develop](https://github.com/mrc-ide/Tapestry/workflows/checks_develop/badge.svg)](https://github.com/mrc-ide/Tapestry/actions)
  [![checks_master](https://github.com/mrc-ide/Tapestry/workflows/checks_master/badge.svg)](https://github.com/mrc-ide/Tapestry/actions)
  <!-- badges: end -->

# Tapestry

Malaria infections often contain multiple genotypes, and when sequenced together
these produce a complex signal that is a mixture of the individual genotypes.
Building on the framework of earlier programs like DEploid and DEploidIBD,
Tapestry attempts to pull these individual genotypes apart by exploiting allele
frequency imbalances within a sample, while simultaneously estimating segments
of identity by descent (IBD) between sequences. Unlike previous programs,
Tapestry uses advanced MCMC methods to ensure that results are robust even for
high complexity of infection (COI).
