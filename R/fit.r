#' @import stats

#' @export
#' @title Maximum Likelihood Estimation of the Sinba Model Parameters
#'
#' @description
#' `fit_sinba()` searches for the maximum likelihood estimate
#' of a given model for two traits,
#' assuming a single birth for each trait.
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
#' @param reps Number of replicates for the Monte Carlo integration.
fit_sinba <- function(tree, x, y, model = "IND", fixed = NULL, reps = 100) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
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
    stop("fit_sinba: `x` and `y` must be binary traits.")
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
      stop("fit_sinba: `fixed` must be an square matrix.")
    }
    if (nrow(fixed) != 4) {
      stop("fit_sinba: `fixed` must be an square matrix.")
    }
    rownames(fixed) <- xy_levels
    colnames(fixed) <- xy_levels

    youngest <- youngest_event(tree, x, y)
    root <- sinba_root(tree, x, y)
    births <- list()
    births[[1]] <- birth_event(tree, youngest[[1]], reps)
    births[[2]] <- birth_event(tree, youngest[[2]], reps)

    xt <- tree_to_cpp(tree)
    cond <- init_tree_conditionals(tree, xy)

    l <- sinba_mc_like(tree, fixed, root, births, reps, xt, cond)

    pi <- rep(0.25, 4)
    names(pi) <- xy_levels
    obj <- list(
      logLik = l,
      states = xy_levels,
      pi = pi,
      root.prior = "flat",
      data = xy,
      tree = tree
    )
    class(obj) <- "fixed_sinba"
    return(obj)
  }

  mQ <- model_matrix(model)
  k <- max(mQ)
  rownames(mQ) <- xy_levels
  colnames(mQ) <- xy_levels
}

# sinba_mc_like calculates the likelihood
# using a Monte Carlo integration.
sinba_mc_like <- function(t, Q, root, births, reps, xt, cond) {
  # make sure that Q matrix is valid
  Q[1, 4] <- 0
  Q[2, 3] <- 0
  Q[3, 2] <- 0
  Q[4, 1] <- 0
  Q <- normalize_Q(Q)

  root_id <- length(t$tip.label) + 1
  likes <- c()
  for (rp in seq_len(reps)) {
    for (r in seq_len(length(root))) {
      if (root[r] != 1) {
        likes <- c(likes, -Inf)
        next
      }

      b1 <- births[[1]][[rp]]
      b2 <- births[[2]][[rp]]

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
      st <- as.integer(active_status(t, ev$first$node, ev$second$node))
      root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))

      l <- sinba_conditionals(
        xt$anc, xt$desc, st, xt$branch,
        cond,
        ev$first$age, ev$second$age,
        root_Q, ev$first$Q, ev$second$Q
      )
      likes <- c(likes, l[root_id, r])
    }
  }

  mx <- max(likes)
  l <- log(sum(exp(likes - mx))) - log(length(likes)) + mx
  return(l)
}

# sinba_root return the valid states
# for the root in the sinba model.
sinba_root <- function(t, x, y) {
  root <- c()
  for (i in levels(x)) {
    n <- youngest_birth_node(t, x, i)
    for (j in levels(y)) {
      s <- paste(i, j, sep = "|")
      root[s] <- 0
      if (n != length(t$tip.label) + 1) {
        next
      }
      m <- youngest_birth_node(t, y, j)
      if (m != length(t$tip.label) + 1) {
        next
      }
      root[s] <- 1
    }
  }
  return(root)
}

# youngest_event returns the youngest possible node
# for each trait.
youngest_event <- function(t, x, y) {
  youngest <- list()
  youngest[[1]] <- youngest_node(t, x)
  youngest[[2]] <- youngest_node(t, y)
  return(youngest)
}

# youngest_node returns the youngest event node
# for a given trait.
youngest_node <- function(t, x) {
  n <- length(t$tip.label) + 1
  for (i in levels(x)) {
    e <- youngest_birth_node(t, x, i)
    if (node_age(t, n) < node_age(t, e)) {
      n <- e
    }
  }
  return(n)
}

# youngest_birth_node returns the earliest node
# in which the process can born.
youngest_birth_node <- function(t, x, state) {
  ns <- rep(FALSE, nrow(t$edge) + 1)
  for (i in seq(from = nrow(t$edge), to = 1)) {
    n <- t$edge[i, 2]
    if (n <= length(t$tip.label)) {
      obs <- x[t$tip.label[n]]
      if (is.na(obs)) {
        # unobserved states are assumed as presence of the state
        ns[n] <- TRUE
      } else if (obs == state) {
        ns[n] <- TRUE
      }
    }
    a <- t$edge[i, 1]
    ns[a] <- ns[a] | ns[n]
  }

  x <- dollo_uppass(t, length(t$tip.label) + 1, ns)
  if (x < 1) {
    x <- length(t$tip.label) + 1
  }
  return(x)
}

# dollo_uppass returns the first node
# in which the state born.
dollo_uppass <- function(t, n, ns) {
  # the node does not have the state
  if (!ns[n]) {
    return(0)
  }

  # the node is a tip
  if (n <= length(t$tip.label)) {
    return(0)
  }

  desc <- which(t$edge[, 1] == n)
  ws <- 0
  for (i in desc) {
    if (ns[t$edge[i, 2]]) {
      ws <- ws + 1
    }
  }
  if (ws > 1) {
    # more than one descendant has the state
    # so the state in the node is present
    return(n)
  }

  for (i in desc) {
    x <- dollo_uppass(t, t$edge[i, 2], ns)
    if (x > 0) {
      return(x)
    }
  }
  return(0)
}

# birth_event returns the age and node locations
# of one or more birth events.
birth_event <- function(t, n, num = 1) {
  e <- which(t$edge[, 2] == n)

  births <- list()

  if (length(e) == 0) {
    # all events are on the root
    for (i in seq(num)) {
      b <- list(
        node = n,
        age = 0
      )
      births[[i]] <- b
    }
    return(births)
  }

  # first event is always in the same branch
  # as the youngest node
  b <- list(
    node = n,
    age = stats::runif(1, max = t$edge.length[e])
  )
  births[[1]] <- b
  if (num == 1) {
    return(births)
  }

  age <- node_age(t, n)
  path <- path_to_node(t, n)
  for (i in seq(from = 2, to = num)) {
    st <- stats::runif(1, max = age)
    cn <- path[0]
    for (j in path) {
      cn <- j
      ce <- which(t$edge[, 2] == j)
      if (length(ce) == 0) {
        next
      }
      if (t$edge.length[ce] > st) {
        break
      }
      st <- st - t$edge.length[ce]
    }
    b <- list(
      node = cn,
      age = st
    )
    births[[i]] <- b
  }
  return(births)
}

# scenario returns the trait activation scenario
# given the root state
# and the first active trait,
scenario <- function(r, tr) {
  if (r == 1) {
    # root 00
    if (tr == 1) {
      # first trait active: 00 <-> 10
      return("13")
    }
    # second trait active: 00 <-> 01
    return("12")
  }
  if (r == 2) {
    # root 01
    if (tr == 1) {
      # first trait active: 01 <-> 11
      return("24")
    }
    # second trait active: 01 <-> 00
    return("12")
  }
  if (r == 3) {
    # root 10
    if (tr == 1) {
      # first trait active: 10 <-> 00
      return("13")
    }
    # second trait active: 10 <-> 11
    return("34")
  }
  # root 11
  if (tr == 1) {
    # first trait active: 11 <-> 01
    return("24")
  }
  # second trait active: 11 <-> 10
  return("34")
}

# semi_active_Q returns the semiactive Q
# for a given scenario.
semi_active_Q <- function(sc, Q) {
  semi <- matrix(rep(0, length(Q)), nrow = nrow(Q))
  if (sc == "12") {
    # second trait active: 00 <-> 01
    semi[1, 2] <- Q[1, 2]
    semi[2, 1] <- Q[2, 1]
    return(semi)
  }
  if (sc == "13") {
    # first trait active: 00 <-> 10
    semi[1, 3] <- Q[1, 3]
    semi[3, 1] <- Q[3, 1]
    return(semi)
  }
  if (sc == "24") {
    # first trait active: 01 <-> 11
    semi[2, 4] <- Q[2, 4]
    semi[4, 2] <- Q[4, 2]
    return(semi)
  }
  if (sc == "34") {
    # second trait active: 10 <-> 11
    semi[3, 4] <- Q[3, 4]
    semi[4, 3] <- Q[4, 3]
    return(semi)
  }
  msg <- sprintf("semi_active_Q: unknown scenario '%s'", sc)
  stop(msg)
}

# tree_to_cpp transforms a tree into a set of vectors
# to be used in a CPP function.
tree_to_cpp <- function(t) {
  return(list(
    anc = as.integer(t$edge[, 1] - 1),
    desc = as.integer(t$edge[, 2] - 1),
    branch = t$edge.length
  ))
}

# init_tree_conditionals creates the tree conditionals
# from a tree
# and its observations.
init_tree_conditionals <- function(t, xy) {
  cond <- matrix(0, nrow = length(t$tip.label) + t$Nnode, ncol = 4)
  xy_levels <- levels(xy)
  for (i in seq_len(length(t$tip.label))) {
    cond[i, ] <- -Inf
    v <- xy[t$tip.label[i]]
    for (j in seq_len(length(xy_levels))) {
      if (v == xy_levels[j]) {
        cond[i, j] <- 0
      }
    }
  }
  colnames(cond) <- levels(xy)
  return(cond)
}

# active_status returns a vector with the active status
# of each node.
active_status <- function(t, n1, n2) {
  st <- rep(0, length(t$tip.label) + t$Nnode)
  for (i in seq_len(length(st))) {
    if (i == n2) {
      st[i] <- 2
      next
    }
    if (i == n1) {
      st[i] <- 1
      next
    }
    if (is_parent(t, n2, i)) {
      st[i] <- 2
      next
    }
    if (is_parent(t, n1, i)) {
      st[i] <- 1
    }
  }
  return(st)
}
