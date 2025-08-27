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
#' @param model Model of evolution for the traits.
#'   By default it uses the independent model ("IND").
#'   The standard model for correlated traits
#'   can be referred as "CORR", "DEP", "xy", or "ARD".
#'   In the "ER" model both traits have equal rates
#'   in any direction;
#'   the "ER2" model also has equal rates,
#'   but rates are different for each trait;
#'   in the "ERs" model the rates of state transitions are equal
#'   but can be different depending on the state.
#'   If the "SYM" model changes between states are equal.
#'   There a two full dependant models,
#'   "x" for a model in which trait x depends on y;
#'   and "y" in which trait y depends on x.
#'   The "coll" model collapse (i.e., removes)
#'   entries for unobserved traits.
#'   In the "sCORR" model,
#'   rates are correlated by the state of the other trait.
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
    tree, data, model = "IND",
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
  obj$model <- sprintf("Pagel + %s", model)
  return(obj)
}
