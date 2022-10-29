
# make_IBD_transition_lookup.R
#
# Author: Bob Verity
# Date: 2022-10-28
#
# Purpose:
# If K is the COI of a sample then the number of possible IBD configurations is
# given by the Kth Bell number, which grows very quickly. These configurations are
# equivalent to all the ways that K elements can be partitioned into up to K
# groups. We want to enumerate all these partitions, then make a hidden Markov
# model (HMM) over IBD states. The transition matrix of this HMM should give the
# rate (or the probability) of transitioning from any IBD state to any other.
#
# Partitions can be enumerated using the Restrictive Growth Function (RGF)
# method, in which each partition is represented as a numeric sequence. Elements
# that have the same numeric value are in the same group. For example, all
# possible IBD staets for K=5 can be enumerated as follows:
# {1,1,1}
# {1,1,2}
# {1,2,1}
# {1,2,2}
# {1,2,3}
# We have functions for moving from any given RGF sequence to it's index and
# vice versa.
#
# To reduce the number of parameters in the model, we assume a single coalescent
# rate parameter between each pair of strains. This leads to K-choose-2
# parameters, which grows much more slowly than the number of partitions. Now we
# need to work out how this smaller number of parameters maps to the larger
# transition matrix. We can do this by enumerating all partitions, then working
# out which coalescent events change this partition to another partition. The
# end result is a dataframe giving the index of every starting partition, the
# end partition that it moves to, and which strains need to coalesce for this to
# happen. These results are saved to file to be loaded into Tapestry.
#
# ------------------------------------------------------------------

library(dplyr)

# --------------------------------
# returns the number of RGF sequences that:
# 1) are length l
# 2) assume the largest value seen before the start of this sequence is j.
# When j=1 this is simply the number of possible partitions of (l+1) objects
# into distinct groups, which is given by the (l+1)th Bell number, but when j>1
# we get a larger number of possible partitions. For example, when l=2 and j=2
# we get the following sequences(l, j):
# {1,1}
# {1,2}
# {1,3}
# {2,1}
# {2,2}
# {2,3}
# {3,1}
# {3,2}
# {3,3}
# {3,4}
# hence our function should return 10. This is solved recursively by noting the
# structure in the above set of sequences, which consists of 3 sequences
# starting with 1, then 3 sequences starting with 2, then 4 sequences starting
# with 3. To generalise, we see all sequences(l-1, j), starting with the numbers
# 1:j, and additionally we see all sequences(l-1, j+1) starting with the number
# j+1. This creates our recursion in the number of sequences.
seq_count <- function(l, j) {
  if (l == 0) {
    return(1)
  }
  ret <- j * seq_count(l-1, j) + seq_count(l-1, j+1)
  return(ret)
}

# --------------------------------
# convert an RGF sequence into its unique integer index. This is done by looping
# through values and counting the number of seqences that are ruled out at each
# step. For example, imagine all sequences of length 3:
# {1,1,1}
# {1,1,2}
# {1,2,1}
# {1,2,2}
# {1,2,3}
# and that are target sequence is z=c(1,2,1). When we hit the value z[2]=2, this
# rules out the sequences {1,1,1} and {1,1,2}, so the index must be at least 3.
# Continuing in this way we eventually arrive at the unique index of a RGF
# sequence
get_partition_no <- function(z) {
  n <- length(z)
  if (n == 1) {
    return(1)
  }
  ret <- 1
  for (i in 2:n) {
    # get greatest value in z up to this point
    j <- max(z[1:(i-1)])
    
    # for any value of z[i] from 1:j there is the same number of sub-sequences
    # in the remaining (n-i) values. This number is given by:
    subseq <- seq_count(n-i, j)
    
    # we can rule out (z[i]-1) blocks of this size
    ret <- ret + (z[i] - 1) * subseq
  }
  return(ret)
}

# --------------------------------
# get RGF sequence of length n at index position x. This uses a similar logic to
# get_partition_no(), but now filling in the sequence z such that the index
# never exceeds x as we add new values. Continuing in this way we eventually
# arrive at the unique RGF sequence corresponding to an index
get_partition_from_no <- function(n, x) {
  z <- rep(1, n)
  x_lower <- 0
  for (i in 2:n) {
    m <- max(z[1:(i-1)])
    v <- x_lower + (0:m) * seq_count(n-i, m)
    w <- max(which(v < x))
    z[i] <- w
    x_lower <- v[w]
  }
  return(z)
}

# ------------------------------------------------------------------

# length of partition sequence (equal to COI)
n <- 8

# create matrix of all pairwise comparisons for indexing
m_pairs <- cbind(rep(1:(n-1), times = (n-1):1),
                 unlist(mapply(function(x) x:n, 2:n)))

# create matrix of all partitions
partition_mat <- t(mapply(function(x) {
  get_partition_from_no(n,x)
}, 1:seq_count(n-1,1)))


# populate results. First value in each row gives the index of the starting
# partition, second value gives the index of the partition that it transitions
# to. Third and fourth values gives the two strains that must coalesce for this
# to happen. Fifth value gives the same information as third and fourth
# together, but now in linear sequence over all n-choose-2 values
ret_list <- list()
for (i in 1:nrow(partition_mat)) {
  x <- partition_mat[i,]
  w <- which(x[m_pairs[,1]] != x[m_pairs[,2]])
  if (any(w)) {
    ret_list[[i]] <- matrix(NA, length(w), 5)
    for (j in seq_along(w)) {
      this_pair <- m_pairs[w[j],]
      y <- x
      y[y == y[this_pair[1]]] <- y[this_pair[2]]
      y <- match(y, unique(y))
      ret_list[[i]][j,] <- c(i, get_partition_no(y), this_pair, w[j])
    }
  }
}

# convert to data.frame
ret_df <- do.call(rbind, ret_list) %>%
  as.data.frame() %>%
  set_names(c("S_from", "S_to", "coal_from", "coal_to", "param"))

head(ret_df)

# explore how many values have to be updated given a specific parameter update.
# This is (perhaps unsurprisingly) the full dimension of the data.frame divided
# by the number of parameters, which is given by choose(n, 2)
ret_df %>%
  dim()

ret_df %>%
  filter(param == 1) %>%
  dim()

# optionally save to file by running the above and saving into a single list
# transition_lookup_list <- list()
# transition_lookup_list[[n]] <- ret_df
# saveRDS(transition_lookup_list, file = "inst/extdata/transition_lookup_list.rds")
