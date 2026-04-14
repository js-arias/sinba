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
#'     "burst":  for the Replicated burst scenario.
#'     "codis":  for the Replicated co-distribution scenario.
#' @param rate_mat Rate matrix for the traits with the full process.
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
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param pi_x The transition probability at the birth of x trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param pi_y The transition probability at the birth of y trait.
#'   If NULL it will set 1.0 for the state 1.
#'
#' @return A list with the tree,
#'   a data frame with the data,
#'   the births,
#'   and the Q matrix used for the simulation.
sim_sinba <- function(
    tree, rate_mat = NULL, births = NULL,
    model = NULL,
    pi_x = NULL, pi_y = NULL,
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
  if (scenario == "codis") {
    return(replicated_codistribution_scenario(tree))
  }
  if (scenario == "burst") {
    return(replicated_burst_scenario(tree))
  }

  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "sim_sinba: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits == 1) {
    return(sim_sinba_single(tree, rate_mat, births, model, pi_x))
  }

  t <- phylo_to_sinba(tree)

  if (is.null(rate_mat)) {
    stop("sim_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("sim_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop("sim_sinba: `rate_mat` must have the same size as the `model`")
  }

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$trait_states[["x"]]$states)
  }
  if (length(pi_y) == 0) {
    pi_y <- default_pi_vector(model$trait_states[["y"]]$states)
  }
  if (length(pi_x) != length(model$trait_states[["x"]]$states)) {
    stop("sim_sinba: invalid pi_x: size different to number of states")
  }
  if (length(pi_y) != length(model$trait_states[["y"]]$states)) {
    stop("sim_sinba: invalid pi_y: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }
  if (sum(pi_y) != 0) {
    pi_y <- pi_y / sum(pi_y)
  }

  if (is.null(births)) {
    max_size <- length(t$tip) * 0.75
    min_size <- length(t$tip) * 0.25
    n2 <- 0
    while (TRUE) {
      n2 <- sample(seq_len(length(t$parent)), size = 1)
      # skip root
      if (n2 == t$root_id) {
        next
      }
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
    if (n == t$root_id) {
      next
    }
    if (is.null(births[[i]]$age)) {
      a <- runif(1, max = t$br_len[n])
      births[[i]]$age <- a
    }
    if (births[[i]]$age > t$br_len[n]) {
      stop(sprintf("sim_sinba: expecting `age` value for birth %d", i))
    }
  }
  if ((births[[1]]$node != t$root_id) &&
    (births[[1]]$node == births[[2]]$node)) {
    if (births[[1]]$age > births[[2]]$age) {
      a <- births[[2]]$age
      births[[2]]$age <- births[[1]]$age
      births[[1]]$age <- a
    }
  }

  rate_mat <- normalize_Q(rate_mat)
  semi_mat <- build_semi_active_Q(model, "13", rate_mat)
  semi_mat <- normalize_Q(semi_mat)

  m_PI_act <- build_pi_matrix(model, "y", active_ancestor_vector("13"))
  root_vector <- rep(0, 4)
  root_vector[1] <- 1
  m_PI_semi <- build_pi_matrix(model, "x", root_vector)

  states <- rep(1, length(t$parent))
  states <- evolve_states(
    t, rate_mat, semi_mat, states, t$root_id, births,
    pi_x, pi_y, m_PI_semi, m_PI_act
  )

  x <- rep(0, length(t$tip))
  y <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    obs <- model$observed[[states[i]]]
    if (obs == "00") {
      next
    }
    if (obs == "01") {
      y[i] <- 1
      next
    }
    if (obs == "10") {
      x[i] <- 1
      next
    }
    # observed state is "11"
    x[i] <- 1
    y[i] <- 1
  }

  return(list(
    tree = tree,
    data = data.frame(t$tip, x, y),
    births = births,
    model = model,
    Q = rate_mat,
    pi_x = pi_x,
    pi_y = pi_y
  ))
}

sim_sinba_single <- function(
    tree, rate_mat, births, model, pi_x) {
  t <- phylo_to_sinba(tree)

  if (is.null(rate_mat)) {
    stop("sim_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("sim_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop("sim_sinba: `rate_mat` must have the same size as the `model`")
  }

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$states)
  }
  if (length(pi_x) != length(model$states)) {
    stop("sim_sinba: invalid pi_x: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }

  if (is.null(births)) {
    max_size <- length(t$tip) * 0.75
    min_size <- length(t$tip) * 0.25
    n <- 0
    while (TRUE) {
      n <- sample(seq_len(length(t$parent)), size = 1)
      # skip root
      if (n == t$root_id) {
        next
      }
      sz <- node_size(t, n)
      if (sz > min_size && sz <= max_size) {
        break
      }
    }
    births <- list()
    births[[1]] <- list(node = n)
  }
  if (length(births) < 1) {
    stop("sim_sinba: `births` should have at least one element")
  }
  if (is.null(births[[1]]$node)) {
    stop("sim_sinba: expecting `node` in an element of `births`")
  }
  n <- births[[1]]$node
  if (n != t$root_id) {
    if (is.null(births[[1]]$age)) {
      a <- runif(1, max = t$br_len[n])
      births[[1]]$age <- a
    }
    if (births[[1]]$age > t$br_len[n]) {
      stop(sprintf("sim_sinba: expecting `age` value for birth"))
    }
  }

  rate_mat <- normalize_Q(rate_mat)
  states <- rep(1, length(t$parent))
  states <- evolve_single_states(
    t, rate_mat, states, t$root_id, births, pi_x
  )

  x <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    obs <- model$observed[[states[i]]]
    if (obs == "0") {
      next
    }
    x[i] <- 1
  }

  return(list(
    tree = tree,
    data = data.frame(t$tip, x),
    births = births,
    model = model,
    Q = rate_mat,
    pi_x = pi_x
  ))
}

evolve_states <- function(
    t, Q, semi_Q, states, n, births,
    pi_semi, pi_act, m_PI_semi, m_PI_act) {
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
      s <- states[a]
      if (n == b1$node) {
        # birth of the semi active process
        l <- b2$age - b1$age
        s <- pick_state_at_birth(s, m_PI_semi, pi_semi)
      }
      s <- sim_evolution(semi_Q, l, s)
      s <- pick_state_at_birth(s, m_PI_act, pi_act)
      states[n] <- sim_evolution(Q, t$br_len[n] - b2$age, s)
    } else if (n == b1$node) {
      # birth of the semi active process
      s <- pick_state_at_birth(s, m_PI_semi, pi_semi)
      states[n] <- sim_evolution(Q, t$br_len[n] - b1$age, s)
    }
  }

  children <- which(t$parent == n)
  for (c in children) {
    states <- evolve_states(
      t, Q, semi_Q, states, c, births,
      pi_semi, pi_act, m_PI_semi, m_PI_act
    )
  }
  return(states)
}

evolve_single_states <- function(
    t, Q, states, n, birth, pi_act) {
  if (n != t$root_id) {
    a <- t$parent[n]
    b <- birth[[1]]

    if (is_parent(t, b$node, n)) {
      states[n] <- sim_evolution(Q, t$br_len[n], states[a])
    } else if (n == b$node) {
      # birth of the process
      s <- 1
      while (TRUE) {
        i <- sample(seq_len(length(pi_act)), size = 1)
        if (runif(1) < pi_act[i]) {
          s <- i
          break
        }
      }
      states[n] <- sim_evolution(Q, t$br_len[n] - b$age, s)
    }
  }

  children <- which(t$parent == n)
  for (c in children) {
    states <- evolve_single_states(
      t, Q, states, c, birth, pi_act
    )
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

pick_state_at_birth <- function(s, m_PI, pi_state) {
  while (TRUE) {
    ns <- sample(seq_len(ncol(m_PI)), size = 1)
    i <- m_PI[s, ns]
    if (i == 0) {
      next
    }
    if (runif(1) < pi_state[i]) {
      return(ns)
    }
  }
}

darwin_scenario <- function(tree, births) {
  t <- phylo_to_sinba(tree)
  if (is.null(births)) {
    max_size <- length(t$tip) * 0.75
    min_size <- length(t$tip) * 0.25
    n1 <- 0
    while (TRUE) {
      n1 <- sample(seq_len(length(t$parent)), size = 1)
      # skip root
      if (n1 == t$root_id) {
        next
      }
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
      # skip root
      if (n1 == t$root_id) {
        next
      }
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

replicated_codistribution_scenario <- function(tree) {
  t <- phylo_to_sinba(tree)

  max_size <- length(t$tip) * 0.30
  min_size <- length(t$tip) * 0.05
  valid_size <- c()
  for (n in seq_len(length(t$parent))) {
    sz <- node_size(t, n)
    if (sz >= min_size && sz <= max_size) {
      valid_size <- c(valid_size, n)
    }
  }

  nodes <- sample(valid_size, size = length(valid_size) * 0.4)
  x <- rep(0, length(t$tip))
  y <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    for (n in nodes) {
      if (!is_parent(t, n, i)) {
        next
      }
      x[i] <- 1
      y[i] <- 1
    }
  }

  n0 <- sample(valid_size, size = 2)
  for (i in seq_len(length(t$tip))) {
    for (n in n0) {
      if (!is_parent(t, n, i)) {
        next
      }
      x[i] <- 0
      y[i] <- 0
    }
  }

  births <- list()
  births[[1]] <- list(node = t$root_id)
  births[[2]] <- list(node = t$root_id)

  return(list(
    tree = tree,
    data = data.frame(t$tip, x, y),
    births = births
  ))
}

replicated_burst_scenario <- function(tree) {
  t <- phylo_to_sinba(tree)

  max_size <- length(t$tip) * 0.30
  min_size <- length(t$tip) * 0.05
  valid_size <- c()
  for (n in seq_len(length(t$parent))) {
    sz <- node_size(t, n)
    if (sz >= min_size && sz <= max_size) {
      valid_size <- c(valid_size, n)
    }
  }

  nodes <- sample(valid_size, size = length(valid_size) * 0.4)
  x <- rep(0, length(t$tip))
  y <- rep(0, length(t$tip))
  for (i in seq_len(length(t$tip))) {
    for (n in nodes) {
      if (!is_parent(t, n, i)) {
        next
      }
      x[i] <- 1
      if (runif(1) < 0.5) {
        y[i] <- 1
      }
    }
  }

  n0 <- sample(valid_size, size = 2)
  for (i in seq_len(length(t$tip))) {
    for (n in n0) {
      if (!is_parent(t, n, i)) {
        next
      }
      x[i] <- 0
      y[i] <- 0
    }
  }

  births <- list()
  births[[1]] <- list(node = t$root_id)
  births[[2]] <- list(node = t$root_id)

  return(list(
    tree = tree,
    data = data.frame(t$tip, x, y),
    births = births
  ))
}
