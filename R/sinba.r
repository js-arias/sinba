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

# is_parent returns true if node p is parent of n.
is_parent <- function(t, p, n) {
  if (p == n) {
    # a node cannot be parent of itself
    return(FALSE)
  }
  while (TRUE) {
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    a <- t$edge[e, 1]
    if (a == p) {
      return(TRUE)
    }
    n <- a
  }
  return(FALSE)
}
