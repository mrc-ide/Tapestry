#------------------------------------------------
#' @title Define empty tapestry_project object
#'
#' @description Define empty tapestry_project object
#'
#' @details TODO
#'
#' @export
#' @examples
#' # TODO
#'
tapestry_project <- function() {
  
  # create some empty data frames for storing results
  GTI_logevidence_model <- data.frame(set = numeric(),
                                      name = character(),
                                      mean = numeric(),
                                      SE = numeric())
  class(GTI_logevidence_model) <- "maverick_GTI_logevidence_model"
  
  GTI_posterior_model <- data.frame(set = numeric(),
                                    name = character(),
                                    Q2.5 = numeric(),
                                    Q50 = numeric(),
                                    Q97.5 = numeric())
  class(GTI_posterior_model) <- "maverick_GTI_posterior_model"
  
  # initialise project with default values
  ret <- list(data = NULL,
              data_processed = NULL,
              parameter_sets = NULL,
              active_set = 0,
              output = list(single_set = list(),
                            all_sets = list(GTI_logevidence_model = GTI_logevidence_model,
                                            GTI_posterior_model = GTI_posterior_model)
              )
  )
  
  # create class and return invisibly
  class(ret) <- "tapestry_project"
  invisible(ret)
}

#------------------------------------------------
#' @title Custom print function for class tapestry_project
#'   
#' @description Custom print function for class tapestry_project, printing a 
#'   summary of the key elements (also equivalent to \code{summary(x)}). To do 
#'   an ordinary \code{print()} of all elements of the project, use the 
#'   \code{print_full()} function.
#'   
#' @param x object of class \code{tapestry_project}
#' @param ... other arguments (ignored)
#'   
#' @export

print.tapestry_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Ordinary print function for class tapestry_project
#'
#' @description Calling \code{print()} on an object of class tapestry_project
#'   results in custom output, therefore this function stands in for the base
#'   \code{print()} function, and is equivalent to running
#'   \code{print(unclass(x))}.
#'
#' @param x object of class \code{tapestry_project}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_class(x, "tapestry_project")
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class tapestry_project
#'   
#' @description Overload summary function for class tapestry_project
#'   
#' @param object object of class \code{tapestry_project}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.tapestry_project <- function(object, ...) {
  
  # print data summary
  cat("TODO - some summary of project here\n")
  
}

#------------------------------------------------
#' @title Determine if object is of class tapestry_project
#'
#' @description Determine if object is of class tapestry_project.
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO
#'
is.tapestry_project <- function(x) {
  inherits(x, "tapestry_project")
}
