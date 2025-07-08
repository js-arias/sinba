#' @import stats
#' @import nloptr
#'
#' @export
#' @title Maximum Likelihood Estimation of the Pagel's model
#'
#' @description
#' `fit_pagel()` searches for the maximum likelihood estimate
#' of a given model for two traits
#' using the Pagel (1994) model.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param x A vector of phenotypic values for a binary trait
#'   for the tips in `tree`.
#' @param y A vector of phenotypic values for a binary trait
#'   for the tips in `tree`.
#' @param model Model of evolution for the traits.
#'   By default it uses the independent model.
#'   Other valid values are "ARD"
#'   or "xy" for a fully correlated model;
#'   "ER" for a model in which both traits have equal rates
#'   in any direction;
#'   "ER2" for an equal rates model,
#'   but rates are different for each trait;
#'   "SYM" for the symmetric model
#'   in which changes between states are equal;
#'   "x" for a model in which trait x depends on y;
#'   and "y" in which trait y depends on x.
#' @param fixed Use a transition matrix with fixed values.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_pagel <- function(tree, x, y, model = "IND", fixed = NULL, opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_pagel: `tree` must be an object of class \"phylo\".")
  }

  if (!is.factor(x)) {
    x <- as.factor(x)
  }
  x_levels <- levels(x)
  if (!is.factor(y)) {
    y <- as.factor(y)
  }
  y_levels <- levels(y)
  if (length(x_levels) != 2 || length(y_levels) != 2) {
    stop("fit_pagel: `x` and `y` must be binary traits.")
  }
  xy <- stats::setNames(
    factor(paste(x, y, sep = "|"),
      levels = sapply(x_levels, paste, y_levels, sep = "|")
    ),
    names(x)
  )
  xy_levels <- levels(xy)

  if (!is.null(fixed)) {
    if (nrow(fixed) != ncol(fixed)) {
      stop("fit_pagel: `fixed` must be an square matrix.")
    }
    if (nrow(fixed) != 4) {
      stop("fit_pagel: `fixed` must be an square matrix.")
    }
    rownames(fixed) <- xy_levels
    colnames(fixed) <- xy_levels

    xt <- tree_to_cpp(tree)
    cond <- init_tree_conditionals(tree, xy)

    l <- pagel_like(tree, fixed, xt, cond)

    pi <- rep(0.25, 4)
    names(pi) <- xy_levels
    obj <- list(
      logLik = l,
      Q = normalize_Q(fixed),
      states = xy_levels,
      pi = pi,
      root.prior = "flat",
      data = xy,
      tree = tree
    )
    class(obj) <- "fixed_pagel"
    return(obj)
  }

  mQ <- model_matrix(model)
  k <- max(mQ)
  rownames(mQ) <- xy_levels
  colnames(mQ) <- xy_levels

  if (is.null(opts)) {
    v <- 1e-04
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      # set the upper bound using the number of replicates
      xtol_abs = rep(v, max(mQ)),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  # the upper bound is a change per branch
  upper <- length(tree$edge) / sum(tree$edge.length)

  # closure for the nloptr function
  like_func <- function(t, d, m) {
    xt <- tree_to_cpp(t)
    cond <- init_tree_conditionals(t, d)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p > upper)) {
        return(Inf)
      }

      Q <- from_model_to_Q(m, p)
      l <- pagel_like(t, Q, xt, cond)
      return(-l)
    })
  }
  fn <- like_func(tree, xy, mQ)

  best <- list()

  # for the equal rates model use Brent
  if (max(mQ) == 1) {
    par <- stats::runif(1, max = upper)
    best <- stats::optim(par, fn,
      method = "Brent", lower = 0, upper = upper
    )
  } else {
    # initial guess using an equal rates model
    erm <- model_matrix("ER")
    erf <- like_func(tree, xy, erm)
    ep <- stats::runif(1, max = upper)
    best <- stats::optim(ep, erf,
      method = "Brent", lower = 0, upper = upper
    )

    par <- rep(best$par[1], max(mQ))
    res <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    best$par <- res$solution
    best$value <- res$objective
  }
  q <- from_model_to_Q(mQ, best$par)
  q <- normalize_Q(q)
  rownames(q) <- xy_levels
  colnames(q) <- xy_levels
  pi <- rep(0.25, 4)
  names(pi) <- xy_levels
  obj <- list(
    logLik = -best$value,
    model = model,
    k = k,
    rates = best$par,
    index.matrix = mQ,
    Q = q,
    states = xy_levels,
    pi = pi,
    root.prior = "flat",
    data = xy,
    tree = tree
  )
  class(obj) <- "fit_pagel"
  return(obj)
}

# pagel_like calculates the likelihood
# of a Pagel's model.
pagel_like <- function(t, Q, xt, cond) {
  # make sure that Q matrix is valid
  Q[1, 4] <- 0
  Q[2, 3] <- 0
  Q[3, 2] <- 0
  Q[4, 1] <- 0
  Q <- normalize_Q(Q)

  root_id <- length(t$tip.label) + 1

  l <- full_conditionals(xt$parent, xt$nodes, xt$branch, cond, Q)

  mx <- max(l[root_id, ])
  like <- log(sum(exp(l[root_id, ] - mx))) + mx + log(0.25)
  return(like)
}
