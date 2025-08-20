#' @import stats

#' @export
#' @title Simulate A Pair of Traits Under the Sinba Model
#'
#' @description
#' `sim_sinba` simulates a pair of traits
#' for a given phylogenetic tree
#' a given Q transition matrix,
#' and a pair of birth events.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param scenario A predefined scenario.
#'   Valid values are:
#'     "darwin": for the Darwin´s scenario.
#'     "unrep":  for the Unreplicated burst scenario.
#' @param rate_mat Rate matrix for the traits with the full process.
#' @param semi_mat Rate matrix for the traits with the semi-active process.
#' @param births A list with the two birth events.
#'   Each element of the list must contain the field
#'   `node` with the node at the end of the edge
#'   in which the event happens.
#'   Optionally it can have a field `age`
#'   with the time from the start of the edge
#'   in which the event happens.
#'   If not defined,
#'   the age will selected at random point of the edge.
#'   One of the nodes in the birth pair
#'   must be an ancestor of other node.
#'   If it is NULL,
#'   a random pair of birth events will be selected.
#'
#' @return A list with the tree,
#'   a data frame with the data,
#'   the births,
#'   and the Q and semi-active
#'   matrices used for the simulation.
sim_sinba <- function(
    tree, rate_mat = NULL, semi_mat = NULL, births = NULL,
    scenario = "") {
  if (!inherits(tree, "phylo")) {
    stop("sim_sinba: `tree` must be an object of class \"phylo\".")
  }

  if (scenario == "darwin") {
    return(darwin_scenario(tree, births))
  }
  if (scenario == "unrep") {
    return(unreplicated_burst_scenario(tree, births))
  }

  t <- phylo_to_sinba(tree)

  if (is.null(rate_mat)) {
    stop("sim_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("sim_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != 4) {
    stop("sim_sinba: `rate_mat` should be 4x4")
  }
  if (is.null(semi_mat)) {
    semi_mat <- rate_mat
  }
  if (nrow(semi_mat) != ncol(semi_mat)) {
    stop("sim_sinba: `semi_mat` must be a square matrix")
  }
  if (nrow(semi_mat) != 4) {
    stop("sim_sinba: `semi_mat` should be 4x4")
  }
  rate_mat <- normalize_Q(rate_mat)
  semi_mat <- normalize_Q(semi_active_Q("12", semi_mat))

  if (is.null(births)) {
    max_size <- length(t$tip) * 0.75
    min_size <- length(t$tip) * 0.25
    n2 <- 0
    while (TRUE) {
      n2 <- sample(seq_len(length(t$parent)), size = 1)
      sz <- node_size(t, n2)
      if (sz > min_size && sz <= max_size) {
        break
      }
    }
    path <- path_to_node(t, n2)
    n1 <- path[1]
    if (length(path) > 1) {
      n1 <- sample(path, size = 1)
    }
    births <- list()
    births[[1]] <- list(node = n1)
    births[[2]] <- list(node = n2)
  }
  if (length(births) < 2) {
    stop("sim_sinba: `births` should have at least two elements")
  }
  for (i in 1:2) {
    if (is.null(births[[i]]$node)) {
      stop("sim_sinba: expecting `node` in an element of `births`")
    }
    n <- births[[i]]$node
    if (is.null(births[[i]]$age)) {
      a <- runif(1, max = t$br_len[n])
      births[[i]]$age <- a
    }
    if (births[[i]]$age > t$br_len[n]) {
      stop(sprintf("sim_sinba: expecting `age` value for birth %d", i))
    }
  }
  if (births[[1]]$node == births[[2]]$node) {
    if (births[[1]]$age > births[[2]]$age) {
      a <- births[[2]]$age
      births[[2]]$age <- births[[1]]$age
      births[[1]]$age <- a
    }
  }

  states <- rep(1, length(t$parent))
  states <- evolve_states(t, rate_mat, semi_mat, states, t$root_id, births)

  x <- rep(0, length(t$tip))
  y <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    if (states[i] == 1) {
      next
    }
    if (states[i] == 2) {
      y[i] <- 1
      next
    }
    if (states[i] == 3) {
      x[i] <- 1
      next
    }
    # states[i] is 4
    x[i] <- 1
    y[i] <- 1
  }

  return(list(
    tree = tree,
    data = data.frame(t$tip, x, y),
    births = births,
    Q = rate_mat,
    semi_Q = semi_mat
  ))
}

evolve_states <- function(t, Q, semi_Q, states, n, births) {
  if (n != t$root_id) {
    a <- t$parent[n]
    b1 <- births[[1]]
    b2 <- births[[2]]

    if (is_parent(t, b2$node, n)) {
      states[n] <- sim_evolution(Q, t$br_len[n], states[a])
    } else if (is_parent(t, b1$node, n)) {
      states[n] <- sim_evolution(semi_Q, t$br_len[n], states[a])
    } else if (n == b2$node) {
      # birth of the full process
      l <- b2$node
      if (n == b1$node) {
        # birth of the semi active process
        l <- b2$age - b1$age
      }
      s <- sim_evolution(semi_Q, l, states[a])
      states[n] <- sim_evolution(Q, t$br_len[n] - b2$age, s)
    } else if (n == b1$node) {
      # birth of the semi active process
      states[n] <- sim_evolution(Q, t$br_len[n] - b1$age, states[a])
    }
  }

  children <- which(t$parent == n)
  for (c in children) {
    states <- evolve_states(t, Q, semi_Q, states, c, births)
  }
  return(states)
}

sim_evolution <- function(Q, len, s) {
  t <- 0
  while (Q[s, s] != 0) {
    dt <- stats::rexp(1, rate = -Q[s, s])
    if (t + dt >= len) {
      break
    }
    t <- t + dt
    nx <- s
    while (TRUE) {
      nx <- sample(seq_len(nrow(Q)), size = 1)
      if (nx != s) {
        if (runif(1) < -Q[s, nx] / Q[s, s]) {
          break
        }
      }
    }
    s <- nx
  }
  return(s)
}

darwin_scenario <- function(tree, births) {
  t <- phylo_to_sinba(tree)
  if (is.null(births)) {
    max_size <- length(t$tip) * 0.75
    min_size <- length(t$tip) * 0.25
    n1 <- 0
    while (TRUE) {
      n1 <- sample(seq_len(length(t$parent)), size = 1)
      sz <- node_size(t, n1)
      if (sz > min_size && sz <= max_size) {
        break
      }
    }
    births[[1]] <- list(node = n1)
    births[[2]] <- list(node = n1)
  }
  if (length(births) < 2) {
    stop("sim_sinba: `births` should have at least two elements")
  }
  n <- births[[1]]$node

  x <- rep(0, length(t$tip))
  y <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    if (!is_parent(t, n, i)) {
      next
    }
    x[i] <- 1
    y[i] <- 1
  }

  return(list(
    tree = tree,
    data = data.frame(t$tip, x, y),
    births = births
  ))
}

unreplicated_burst_scenario <- function(tree, births) {
  t <- phylo_to_sinba(tree)
  if (is.null(births)) {
    max_size <- length(t$tip) * 0.75
    min_size <- length(t$tip) * 0.25
    n1 <- 0
    while (TRUE) {
      n1 <- sample(seq_len(length(t$parent)), size = 1)
      sz <- node_size(t, n1)
      if (sz > min_size && sz <= max_size) {
        break
      }
    }
    births[[1]] <- list(node = n1)
    births[[2]] <- list(node = n1)
  }
  if (length(births) < 2) {
    stop("sim_sinba: `births` should have at least two elements")
  }
  n <- births[[1]]$node

  x <- rep(0, length(t$tip))
  y <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    if (!is_parent(t, n, i)) {
      next
    }
    x[i] <- 1
    if (runif(1) < 0.5) {
      y[i] <- 1
    }
  }

  return(list(
    tree = tree,
    data = data.frame(t$tip, x, y),
    births = births
  ))
}
