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
#' @param Q Transition matrix for simulation.
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
#'   births,
#'   and Q matrix used for the simulation,
#'   and factors x and y
#'   with the simulated states on the tips.
sim_sinba <- function(tree, Q, births = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("sim_sinba: `tree` must be an object of class \"phylo\".")
  }

  if (is.null(Q)) {
    stop("sim_sinba: undefined `Q` matrix")
  }
  if (nrow(Q) != ncol(Q)) {
    stop("sim_sinba: `Q` must be an square matrix.")
  }
  if (nrow(Q) != 4) {
    stop("sim_sinba: `Q` must be an square matrix.")
  }

  if (is.null(births)) {
    max_size <- tree$Nnode * 0.75
    min_size <- tree$Nnode * 0.25
    n1 <- 0
    while (TRUE) {
      n1 <- sample(1:tree$Nnode, size = 1) + length(tree$tip.label)
      sz <- node_size(tree, n1)
      if (sz > min_size && sz <= max_size) {
        break
      }
    }
    root <- length(tree$tip.label) + 1
    n2 <- n1
    path <- path_to_node(tree, n1)
    while (TRUE) {
      n2 <- sample(path, size = 1)
      if (n2 != root) {
        break
      }
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
    e <- which(tree$edge[, 2] == births[[i]]$node)
    if (is.null(births[[i]]$age)) {
      a <- runif(1, max = tree$edge.length[e])
      births[[i]]$age <- a
    }
    if (births[[i]]$age > tree$edge.length[e]) {
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

  Q <- normalize_Q(Q)
  semiQ <- normalize_Q(semi_active_Q("12", Q))

  states <- rep(1, length(tree$tip.label) + tree$Nnode)
  states <- evolve_states(
    tree, semiQ, Q, states,
    length(tree$tip.label) + 1, births
  )

  x <- rep("a", length(tree$tip.label))
  y <- rep("c", length(tree$tip.label))
  for (i in seq_len(length(tree$tip.label))) {
    if (states[i] == 1) {
      x[i] <- "a"
      y[i] <- "c"
      next
    }
    if (states[i] == 2) {
      x[i] <- "a"
      y[i] <- "d"
      next
    }
    if (states[i] == 3) {
      x[i] <- "b"
      y[i] <- "c"
      next
    }
    # states[i] is 4
    x[i] <- "b"
    y[i] <- "d"
  }
  x <- as.factor(stats::setNames(x, tree$tip.label))
  y <- as.factor(stats::setNames(y, tree$tip.label))

  return(list(
    tree = tree,
    births = births,
    Q = Q,
    x = x,
    y = y
  ))
}

evolve_states <- function(t, semiQ, Q, states, n, births) {
  e <- which(t$edge[, 2] == n)
  if (length(e) > 0) {
    a <- t$edge[e, 1]
    b1 <- births[[1]]
    b2 <- births[[2]]

    if (is_parent_in_tree(t, b2$node, n)) {
      states[n] <- sim_evolution(Q, t$edge.length[e], states[a])
    } else if (is_parent_in_tree(t, b1$node, n)) {
      states[n] <- sim_evolution(semiQ, t$edge.length[e], states[a])
    } else if (n == b2$node) {
      # birth of the full process
      l <- b2$age
      if (n == b1$node) {
        # birth of the semi active process
        l <- b2$age - b1$age
      }
      s <- sim_evolution(semiQ, l, states[a])
      states[n] <- sim_evolution(Q, t$edge.length[e] - b2$age, s)
    } else if (n == b1$node) {
      # birth of the semi active process
      states[n] <- sim_evolution(Q, t$edge.length[e] - b1$age, states[a])
    }
  }

  children <- which(t$edge[, 1] == n)
  for (e in children) {
    states <- evolve_states(t, semiQ, Q, states, t$edge[e, 2], births)
  }
  return(states)
}

sim_evolution <- function(Q, len, s) {
  t <- 0
  while (TRUE) {
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
