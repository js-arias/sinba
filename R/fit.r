#' @export
#' @title
#' Estimate the Likelihood for a Given Sinba Transition Matrices in Two Traits
#'
#' @description
#' `fixed_sinba()` calculates the likelihood
#' for a given pair of transition matrices
#' for two traits
#' under the Sinba model.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param rate_mat Rate matrix for the traits with the full process.
#' @param births A list with two element,
#'   each element being a list that define a birth event,
#'   with a fields `node` indicating the birth of the trait
#'   and `age` indicating the time from the start of the edge
#'   in which the event happens.
#' @param semi_mat Rate matrix for the traits with the semi-active process.
#'   If it is NULL,
#'   it will use the same parameter values
#'   from the full process.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
fixed_sinba <- function(
    tree, data, rate_mat, births,
    semi_mat = NULL, root = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fixed_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  cond <- init_conditionals(t, data, 2)

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fixed_sinba: invalid size for `root` vector.")
  }
  root <- root / sum(root)

  if (is.null(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != 4) {
    stop("fixed_sinba: `rate_mat` should be 4x4")
  }

  if (is.null(semi_mat)) {
    semi_mat <- rate_mat
  }
  if (nrow(semi_mat) != ncol(semi_mat)) {
    stop("fixed_sinba: `semi_mat` must be a square matrix")
  }
  if (nrow(semi_mat) != 4) {
    stop("fixed_sinba: `semi_mat` should be 4x4")
  }

  if (is.null(births)) {
    stop("fixed_sinba: undefined `births` events")
  }
  if (length(births) < 2) {
    stop("fixed_sinba: two events required in `births`")
  }
  for (i in 1:2) {
    e <- births[[i]]
    if (is.null(e$node)) {
      stop(sprintf(
        "fixed_sinba: expecting field `node` in `births` event %d", i
      ))
    }
    if (is.null(e$age)) {
      stop(sprintf(
        "fixed_sinba: expecting field `age` in `births` event %d", i
      ))
    }
    if (e$age > t$br_len[e$node]) {
      stop(sprintf(
        "fixed_sinba: invalid value for field `age` in `births` event %d", i
      ))
    }
  }

  # if events are not in the same root path
  # the likelihood is 0.
  b1 <- births[[1]]
  b2 <- births[[2]]
  if (b1$node != b2$node) {
    if (!is_parent(t, b1$node, b2$node) && !(is_parent(t, b2$node, b1$node))) {
      obj <- list(
        logLik = -Inf,
        Q = normalize_Q(rate_mat),
        semi_Q = semi_mat,
        births = births,
        states = c("00", "01", "10", "11"),
        root_prior = root,
        data = data,
        tree = tree
      )
      class(obj) <- "fixed_sinba"
      return(obj)
    }
  }

  # check for valid events
  root_prob <- root
  youngest <- youngest_birth_event(t, cond)
  y1 <- youngest[[1]]
  y2 <- youngest[[2]]
  if (root[1] > 0) {
    # ancestral state is 00
    # check is trait 1-1 can be born
    if (y1[2] != b1$node) {
      if (!is_parent(t, b1$node, y1[2])) {
        # invalid root state
        root_prob[1] <- 0
      }
    }
    # check if trait 2-1 can be born
    if (y2[2] != b2$node) {
      if (!is_parent(t, b2$node, y2[2])) {
        # invalid root state
        root_prob[1] <- 0
      }
    }
  }
  if (root[2] > 0) {
    # ancestral state is 01
    # check is trait 1-1 can be born
    if (y1[2] != b1$node) {
      if (!is_parent(t, b1$node, y1[2])) {
        # invalid root state
        root_prob[2] <- 0
      }
    }
    # check if trait 2-0 can be born
    if (y2[1] != b2$node) {
      if (!is_parent(t, b2$node, y2[1])) {
        # invalid root state
        root_prob[2] <- 0
      }
    }
  }
  if (root[3] > 0) {
    # ancestral state is 10
    # check is trait 1-0 can be born
    if (y1[1] != b1$node) {
      if (!is_parent(t, b1$node, y1[1])) {
        # invalid root state
        root_prob[3] <- 0
      }
    }
    # check if trait 2-1 can be born
    if (y2[2] != b2$node) {
      if (!is_parent(t, b2$node, y2[2])) {
        # invalid root state
        root_prob[3] <- 0
      }
    }
  }
  if (root[4] > 0) {
    # ancestral state is 11
    # check is trait 1-0 can be born
    if (y1[1] != b1$node) {
      if (!is_parent(t, b1$node, y1[1])) {
        # invalid root state
        root_prob[4] <- 0
      }
    }
    # check if trait 2-0 can be born
    if (y2[1] != b2$node) {
      if (!is_parent(t, b2$node, y2[1])) {
        # invalid root state
        root_prob[4] <- 0
      }
    }
  }

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (all(root_prob == 0)) {
    obj <- list(
      logLik = -Inf,
      Q = normalize_Q(rate_mat),
      semi_Q = semi_mat,
      births = births,
      states = c("00", "01", "10", "11"),
      root_prior = root,
      data = data,
      tree = tree
    )
    class(obj) <- "fixed_sinba"
    return(obj)
  }

  ev_prob <- 2 * log(prob_birth(t))

  xt <- tree_to_cpp(t)
  lk <- sinba_like(
    t, rate_mat, semi_mat, births, xt, cond,
    log(root_prob), ev_prob
  )

  obj <- list(
    logLik = lk,
    Q = normalize_Q(rate_mat),
    semi_Q = semi_mat,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
    data = data,
    tree = tree
  )
  class(obj) <- "fixed_sinba"
  return(obj)
}

# sinba_like calculates the likelihood
# of the sinba model.
sinba_like <- function(t, Q, semi_Q, births, xt, cond, root_prior, ev_prob) {
  # make sure that Q matrix is valid
  Q[1, 4] <- 0
  Q[2, 3] <- 0
  Q[3, 2] <- 0
  Q[4, 1] <- 0
  Q <- normalize_Q(Q)

  root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))

  likes <- c()
  for (r in seq_len(length(root_prior))) {
    if (is.infinite(root_prior[r])) {
      next
    }

    b1 <- births[[1]]
    b2 <- births[[2]]

    ev <- list()
    n1 <- b1$node
    n2 <- b2$node
    a1 <- t$ages[n1]
    a2 <- t$ages[n2]
    if (n1 == n2) {
      a1 <- b1$age
      a2 <- b2$age
    }

    if (a1 < a2) {
      # first trait is the oldest one
      sc <- scenario(r, 1)
      semi <- semi_active_Q(sc, semi_Q)
      semi <- normalize_Q(semi)
      ev <- list(
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
      semi <- semi_active_Q(sc, semi_Q)
      semi <- normalize_Q(semi)
      ev <- list(
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

    st <- as.integer(active_status(t, ev$first$node, ev$second$node))

    l <- full_sinba_conditionals(
      xt$parent, xt$nodes, st, xt$branch,
      cond,
      ev$first$age, ev$second$age,
      root_Q, ev$first$Q, ev$second$Q
    )
    likes <- c(likes, l[t$root_id, r] + root_prior[r])
  }
  mx <- max(likes)
  lk <- log(sum(exp(likes - mx))) + ev_prob + mx
  return(lk)
}
