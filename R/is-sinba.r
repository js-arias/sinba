#' @export
#' @title
#' Return true if the data set is a full Sinba dataset
#'
#' @description
#' `is_full_sinba()` checks if a dataset is a full Sinba dataset,
#' that is,
#' if for each trait at least one state born
#' at a node different from the root,
#' and,
#' if there are two traits,
#' they form a nested set.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
is_full_sinba <- function(tree, data, model = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop("fit_sinba: `model` must be an object of class \"sinba_model\".")
  }
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
  }

  t <- phylo_to_sinba(tree)

  if (model$traits == 1) {
    et <- encode_traits(t, data, 1)

    root_states <- c("0", "1")
    for (r in root_states) {
      youngest <- youngest_birth_event(t, et, r)
      if (youngest[[1]] != t$root_id) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  et <- encode_traits(t, data, 2)
  root_states <- c("00", "01", "10", "11")
  for (r in root_states) {
    youngest <- youngest_birth_event(t, et, r)
    y1 <- youngest[[1]]
    if (y1 == t$root_id) {
      next()
    }
    y2 <- youngest[[2]]
    if (y2 == t$root_id) {
      next
    }

    # the nodes are equal
    if (y1 == y2) {
      return(TRUE)
    }
    # the nodes are on the same path
    if (is_parent(t, y1, y2) || is_parent(t, y2, y1)) {
      return(TRUE)
    }

    # the nodes are on different paths
    # but share a common ancestor
    # that is not the root
    x <- most_recent_common_ancestor(t, y1, y2)
    if (x != t$root_id) {
      return(TRUE)
    }
  }

  return(FALSE)
}
