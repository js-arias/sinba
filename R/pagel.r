#' @export
#' @title
#' Maximum Likelihood Estimation of the Pagel's Model Under the Sinba Model
#'
#' @description
#' `fit_pagel()` searches for the maximum likelihood estimate
#' of a given model for two traits
#' using the Pagel (1994) model.
#' It is a wrapper of `fit_fixed_births()`
#' with the births at the root.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
#' @param root_method Method for root calculation at the root.
#'   By default it use the root prior.
#'   If set as "FitzJohn" it will use the FitzJohn et al. (2009)
#'   method,
#'   in which ancestral states are weighted by its own likelihood.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_pagel <- function(
    tree, data, model = NULL,
    root = NULL, root_method = "prior", opts = NULL) {
  root_id <- length(tree$tip.label) + 1
  births <- list()
  for (i in 1:2) {
    b <- list(
      node = root_id,
      age = 0
    )
    births[[i]] <- b
  }
  obj <- fit_fixed_births(tree, data, births, model, root, root_method, opts)
  return(obj)
}
