
#------------------------------------------------
#' @title Check that Tapestry package has loaded successfully
#'
#' @description Simple function to check that Tapestry package has loaded 
#'   successfully. Prints "Tapestry loaded successfully!" if so.
#'
#' @export

check_tapestry_loaded <- function() {
  message("Tapestry loaded successfully!")
}

#------------------------------------------------
#' @title Test function to run an example MCMC
#'
#' @description This should be replaced by a more carefully thought out
#'   structure, probably involving defining a parameters dataframe outside of
#'   this function that can be loaded in.
#'
#' @param a vector of alternative allele counts.
#' @param r vector of reference allele counts.
#' @param p vector of allele frequencies per locus.
#' @param burnin the number of burn-in iterations.
#' @param samples the number of sampling iterations.
#' @param beta vector of thermodynamic powers. Final value in the vector should
#'   always be 1.
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100\%. This avoids large amounts of
#'   output being printed to markdown files.
#' @param silent whether to suppress all console output.
#'
#' @importFrom stats setNames
#' @importFrom utils txtProgressBar
#' @export

run_mcmc <- function(a, r, p,
                     burnin = 1e2,
                     samples = 1e3,
                     beta = 1,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  # check inputs
  assert_vector_pos_int(a)
  assert_vector_pos_int(r)
  assert_vector_bounded(p)
  assert_same_length_multiple(a, r, p)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_vector_bounded(beta)
  #assert_eq(beta[length(beta)], 1)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # calculate WSAFs
  x <- a / (a + r)
  
  # make a list of model parameters
  args_params <- list(error_0 = 0.01)
  
  # make a list of MCMC parameters
  args_MCMC <- list(burnin = burnin,
                    samples = samples,
                    beta = beta,
                    pb_markdown = pb_markdown,
                    silent = silent)
  
  # make progress bars
  pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args_progress <- list(pb_burnin = pb_burnin,
                        pb_samples = pb_samples)
  
  # R functions to pass to C++
  args_functions <- list(update_progress = update_progress)
  
  # ---------- run MCMC ----------
  
  # run efficient C++ function
  output_raw <- run_mcmc_cpp(x, args_params, args_MCMC, args_progress, args_functions)
  
  #return(output_raw)
  
  # ---------- process output ----------
  
  # get parameters draws from burn-in phase into data.frame
  df_burnin <- do.call(rbind, output_raw$mu_burnin) %>%
    as.data.frame() %>%
    setNames(sprintf("mu_%s", 1:2)) %>%
    dplyr::mutate(phase = "burnin",
                  iteration = 1:burnin,
                  .before = 1) %>%
    dplyr::mutate(sigma = output_raw$sigma_burnin)
  
  # equivalent for sampling phase
  df_sampling <- do.call(rbind, output_raw$mu_sampling) %>%
    as.data.frame() %>%
    setNames(sprintf("mu_%s", 1:2)) %>%
    dplyr::mutate(phase = "sampling",
                  iteration = burnin + 1:samples,
                  .before = 1) %>%
    dplyr::mutate(sigma = output_raw$sigma_sampling)
  
  # return
  ret <- list(draws = rbind(df_burnin, df_sampling),
              diagnostics = list(MC_accept_burnin = output_raw$MC_accept_burnin / burnin,
                                 MC_accept_sampling = output_raw$MC_accept_sampling / samples))
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}
