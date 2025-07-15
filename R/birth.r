#' @import stats
#' @import nloptr

#' @export
#' @title
#' Maximum Likelihood Estimation of the Sinba Model Parameters Including Births
#'
#' @description
#' `fit_sinba_births()` searches for the maximum likelihood estimate
#' of a given model for two traits,
#' assuming a single birth for each trait.
#' Different from `fit_sinba()`
#' in this function the birth events
#' are also taken as parameters to be optimized.
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
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba_births <- function(tree, x, y, model = "IND", opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba_births: `tree` must be an object of class \"phylo\".")
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
    stop("fit_sinba_birth: `x` and `y` must be binary traits.")
  }
  xy <- stats::setNames(
    factor(paste(x, y, sep = "|"),
      levels = sapply(x_levels, paste, y_levels, sep = "|")
    ),
    names(x)
  )
  xy_levels <- levels(xy)

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
  ages <- node_lengths(tree)
  youngest <- youngest_event(tree, x, y)

  # closure for the nloptr function
  like_func <- function(t, d, dx, dy, m, ye, a) {
    root <- sinba_root(t, dx, dy)

    xt <- tree_to_cpp(t)
    cond <- init_tree_conditionals(t, d)

    # event probability is scaled with the total tree length
    # scaled by the maximum tip distance.
    ev_prob <- 2 * log(max(ages) / sum(t$edge.length))

    paths <- list()
    for (i in 1:2) {
      x <- path_to_node(t, ye[[i]])
      paths[[i]] <- x
    }

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      births <- list()
      for (i in 1:2) {
        if (p[i] > ages[ye[[i]]]) {
          return(Inf)
        }
        n <- get_node_by_age(paths[[i]], a, p[i])
        e <- xt$branch[n] - a[n] + p[i]
        b <- list(
          node = n,
          age = e
        )
        births[[i]] <- b
      }
      Q <- from_model_to_Q(m, p[3:length(p)])
      l <- sinba_birth_like(t, Q, root, births, xt, cond, ev_prob)
      return(-l)
    })
  }
  fn <- like_func(tree, xy, x, y, mQ, youngest, ages)
  par <- c(1e-6, 1e-6, rep(1, max(mQ)))
  res <- nloptr::nloptr(
    x0 = par,
    eval_f = fn,
    opts = opts
  )
  rates <- res$solution[3:length(res$solution)]
  q <- from_model_to_Q(mQ, rates)
  q <- normalize_Q(q)
  rownames(q) <- xy_levels
  colnames(q) <- xy_levels
  pi <- rep(0.25, 4)
  names(pi) <- xy_levels

  xt <- tree_to_cpp(tree)
  paths <- list()
  for (i in 1:2) {
    x <- path_to_node(tree, youngest[[i]])
    paths[[i]] <- x
  }
  births <- list()
  for (i in 1:2) {
    n <- get_node_by_age(paths[[i]], ages, res$solution[i])
    e <- xt$branch[n] - ages[n] + res$solution[i]
    b <- list(
      node = n,
      age = e
    )
    births[[i]] <- b
  }

  obj <- list(
    logLik = -res$objective,
    model = model,
    k = k + 2, # include the birth events
    rates = rates,
    index.matrix = mQ,
    Q = q,
    births = births,
    states = xy_levels,
    pi = pi,
    root.prior = "flat",
    data = xy,
    tree = tree
  )
  class(obj) <- "fit_sinba_births"
  return(obj)
}

# sinba_birth_like calculates the likelihood
# with fixed births.
sinba_birth_like <- function(t, Q, root, births, xt, cond, t_prob) {
  # make sure that Q matrix is valid
  Q[1, 4] <- 0
  Q[2, 3] <- 0
  Q[3, 2] <- 0
  Q[4, 1] <- 0
  Q <- normalize_Q(Q)

  root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))
  root_prior <- log(0.25)

  root_id <- length(t$tip.label) + 1
  likes <- c()
  for (r in seq_len(length(root))) {
    if (root[r] != 1) {
      likes <- c(likes, -Inf)
      next
    }

    b1 <- births[[1]]
    b2 <- births[[2]]

    ev <- list()
    n1 <- b1$node
    n2 <- b2$node
    a1 <- node_age(t, n1)
    a2 <- node_age(t, n2)
    if (n1 == n2) {
      a1 <- b1$age
      a2 <- b2$age
    }

    if (a1 < a2) {
      # first trait is the oldest one
      sc <- scenario(r, 1)
      semi <- semi_active_Q(sc, Q)
      semi <- normalize_Q(semi)
      ev <- list(
        root = r,
        first = list(
          # birth of trait 1
          node = n1,
          age = b1$age,
          trait = 1,
          Q = semi
        ),
        second = list(
          # birth of trait 2
          node = n2,
          age = b2$age,
          trait = 2,
          Q = Q
        )
      )
    } else {
      # second trait is the oldest one
      sc <- scenario(r, 2)
      semi <- semi_active_Q(sc, Q)
      semi <- normalize_Q(semi)
      ev <- list(
        root = r,
        first = list(
          # birth of trait 2
          node = n2,
          age = b2$age,
          trait = 2,
          Q = semi
        ),
        second = list(
          # birth of trait 1
          node = n1,
          age = b1$age,
          trait = 1,
          Q = Q
        )
      )
    }

    st <- as.integer(active_status(xt$parent, ev$first$node, ev$second$node))

    l <- full_sinba_conditionals(
      xt$parent, xt$nodes, st, xt$branch,
      cond,
      ev$first$age, ev$second$age,
      root_Q, ev$first$Q, ev$second$Q
    )
    likes <- c(likes, l[root_id, r] + root_prior)
  }

  mx <- max(likes)
  lk <- log(sum(exp(likes - mx))) + t_prob + mx
  return(lk)
}
