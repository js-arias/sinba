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
