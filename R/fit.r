#' @export
#' @title
#' Maximum Likelihood Estimation of the Sinba Model
#'
#' @description
#' `fit_sinba()` searches for the maximum likelihood estimate
#' using the Sinba model
#' for one or two traits.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   A second,
#'   and optionally a third column,
#'   should contain the data,
#'   coded as 0 and 1.
#' @param root Root prior probabilities.
fit_sinba <- function(tree, data, root = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)
  cond <- init_conditionals(t, data)
  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_sinba: invalid size for `root` vector.")
  }
  root <- root / sum(root)
  youngest <- youngest_birth_event(t, cond)
  print(youngest)

  # nloptr function
  optFn <- function(p) {
    if (any(p < 0)) {
      return(Inf)
    }
    if (any(p[3:length(p)] > 1000)) {
      return(Inf)
    }
  }
}
