# model_matrix returns a design matrix for a given model.
model_matrix <- function(model = "") {
  if (model == "ER") {
    # equal rates model
    return(matrix(c(
      0, 1, 1, 0,
      1, 0, 0, 1,
      1, 0, 0, 1,
      0, 1, 1, 0
    ), nrow = 4, byrow = TRUE))
  }
  if (model == "ER2") {
    # equal rates model
    return(matrix(c(
      0, 1, 2, 0,
      1, 0, 0, 2,
      2, 0, 0, 1,
      0, 2, 1, 0
    ), nrow = 4, byrow = TRUE))
  }
  if (model == "SYM") {
    # symmetric model
    return(matrix(c(
      0, 1, 2, 0,
      1, 0, 0, 3,
      2, 0, 0, 4,
      0, 3, 4, 0
    ), nrow = 4, byrow = TRUE))
  }
  if ((model == "xy") || (model == "ARD")) {
    # correlated model
    return(matrix(c(
      0, 1, 2, 0,
      3, 0, 0, 4,
      5, 0, 0, 6,
      0, 7, 8, 0
    ), nrow = 4, byrow = TRUE))
  }
  if (model == "x") {
    # correlated model,
    # x depends on y.
    return(matrix(c(
      0, 1, 2, 0,
      3, 0, 0, 4,
      5, 0, 0, 1,
      0, 6, 3, 0
    ), nrow = 4, byrow = TRUE))
  }
  if (model == "y") {
    # correlated model,
    # y depends on x.
    return(matrix(c(
      0, 1, 2, 0,
      3, 0, 0, 2,
      4, 0, 0, 5,
      0, 4, 6, 0
    ), nrow = 4, byrow = TRUE))
  }

  # by default returns the independent model
  return(matrix(c(
    0, 1, 2, 0,
    3, 0, 0, 2,
    4, 0, 0, 1,
    0, 4, 3, 0
  ), nrow = 4, byrow = TRUE))
}

# normalize_Q set the diagonal of the Q matrix.
normalize_Q <- function(Q) {
  for (i in seq_len(nrow(Q))) {
    s <- 0
    for (j in seq_len(ncol(Q))) {
      if (i == j) {
        next
      }
      s <- s + Q[i, j]
    }
    Q[i, i] <- -s
  }
  return(Q)
}

# from_model_to_Q sets a Q matrix
# from a model matrix
# and a set of parameter values.
from_model_to_Q <- function(model, par) {
  Q <- matrix(nrow = nrow(model), ncol = ncol(model))
  for (i in seq_len(nrow(model))) {
    for (j in seq_len(ncol(model))) {
      Q[i, j] <- 0
      v <- model[i, j]
      if (v == 0) {
        next
      }
      Q[i, j] <- par[v]
    }
  }
  return(Q)
}

# node_age returns the elapsed time
# since the root time.
node_age <- function(t, n) {
  s <- 0
  while (TRUE) {
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    s <- s + t$edge.length[e]
    n <- t$edge[e, 1]
  }
  return(s)
}

# path_to_node returns the path from the root
# to a given node.
path_to_node <- function(t, n) {
  x <- c()
  while (TRUE) {
    x <- c(n, x)
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    n <- t$edge[e, 1]
  }
  return(x)
}

# length_to_root returns the total length from a node
# towards the root.
length_to_root <- function(t, n) {
  s <- 0
  while (TRUE) {
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    s <- s + t$edge.length[e]
    n <- t$edge[e, 1]
  }
  return(s)
}

# is_parent returns true if node p is parent of n.
is_parent <- function(anc, p, n) {
  if (p == n) {
    # a node cannot be parent of itself
    return(FALSE)
  }
  while (TRUE) {
    a <- anc[n] + 1
    if (a == 0) {
      break
    }
    if (a == p) {
      return(TRUE)
    }
    n <- a
  }
  return(FALSE)
}

# node_length returns a list with the total distance
# from the node to the root.
node_lengths <- function(t) {
  nl <- c()
  for (i in seq_len(length(t$tip.label) + t$Nnode)) {
    l <- length_to_root(t, i)
    nl <- c(nl, l)
  }
  return(nl)
}

# get_node_by_age returns a node in a path
# using the age a of the node
get_node_by_age <- function(p, ages, a) {
  for (i in p) {
    if (ages[i] >= a) {
      return(i)
    }
  }
  return(-1)
}
