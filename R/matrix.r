#' @export
#' @title Build a Model Matrix
#'
#' @description
#' `new_model()` creates a new model matrix.
#'
#' @param model Model of evolution for the traits.
#'   By default it returns the Pagel's independent model ("IND")
#'   (i.e., this is equivalent to using two Q matrices).
#'   The "CORR" model is the Pagel's correlated model,
#'   in which each transition has a different parameter.
#'   In the "ER" model both traits have equal rates
#'   in any direction;
#'   the "ER2" model also has equal rates,
#'   but rates are different for each trait;
#'   in the "ERs" model the rates of state transitions are equal
#'   but can be different depending on the state.
#'   If the "SYM" model changes between states are equal.
#'   There a two full dependant models,
#'   "x" for a model in which trait x
#'   (the first trait)
#'   depends on y;
#'   and "y" in which trait y
#'   (the second trait)
#'   depends on x.
#'   In the "sCORR" model,
#'   rates are correlated by the state of the other trait.
#'   In the "CMK" model,
#'   is like the symmetrical model,
#'   but changes of more than one state are allowed.
new_model <- function(model = "") {
  names <- c("IND", "CORR", "ER", "ER2", "ERs", "SYM", "sCORR", "x", "y", "CMK")
  if (!(model %in% names)) {
    model <- "IND"
  }
  states <- c("0,0", "0,1", "1,0", "1,1")
  m <- model_matrix(model)
  obj <- list(
    name = model,
    model = m,
    states = states
  )
  class(obj) <- "sinba_model"
  return(obj)
}

# model_matrix returns a design matrix for a given model.
model_matrix <- function(model = "") {
  if (model == "CORR") {
    # correlated model
    return(matrix(c(
      0, 1, 2, 0,
      3, 0, 0, 4,
      5, 0, 0, 6,
      0, 7, 8, 0
    ), nrow = 4, byrow = TRUE))
  }
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
    # equal rates model for each trait
    return(matrix(c(
      0, 1, 2, 0,
      1, 0, 0, 2,
      2, 0, 0, 1,
      0, 2, 1, 0
    ), nrow = 4, byrow = TRUE))
  }
  if (model == "ERs") {
    # equal rates model for each state
    return(matrix(c(
      0, 1, 1, 0,
      2, 0, 0, 1,
      2, 0, 0, 1,
      0, 2, 2, 0
    ), nrow = 4, byrow = TRUE))
  }
  if (model == "IND") {
    # independent model
    return(matrix(c(
      0, 1, 2, 0,
      3, 0, 0, 2,
      4, 0, 0, 1,
      0, 4, 3, 0
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
  if (model == "sCORR") {
    # simplified correlated model
    return(matrix(c(
      0, 1, 1, 0,
      1, 0, 0, 2,
      1, 0, 0, 2,
      0, 2, 2, 0
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
  if (model == "CMK") {
    # correlated Mk model
    # the correlation is given by simultaneous change
    # in both states.
    return(matrix(c(
      0, 1, 2, 3,
      1, 0, 3, 2,
      2, 3, 0, 1,
      3, 2, 1, 0
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

#' @export
#' @title Add a New Parameter to a Model
#'
#' @description
#' `model_add()` reads a model definition
#' and add a new parameter at the indicated cell
#' of the model matrix.
#'
#' @param model A model definition as build with `new_model()`.
#' @param cell A vector with the row and column of the new parameter.
model_add <- function(model, cell) {
  if (!inherits(model, "sinba_model")) {
    stop("model_add: `model` must be an object of class \"sinba_model\".")
  }
  m <- model$model
  if (length(cell) < 2) {
    stop("model_equate: `cell` should at least have two elements.")
  }
  if (cell[1] < 1 || cell[1] > nrow(m)) {
    return(model)
  }
  if (cell[2] < 1 || cell[2] > ncol(m)) {
    return(model)
  }
  if (cell[1] == cell[2]) {
    return(model)
  }
  m[cell[1], cell[2]] <- max(m) + 1
  model$name <- "user defined"
  model$model <- m
  return(model)
}

#' @export
#' @title Make Equivalent Two or More Parameters
#'
#' @description
#' `model_equate()` read a model definition
#' and make the indicated parameters
#' to be identical.
#'
#' @param model A model definition as build with `new_model()`.
#' @param params A vector with the parameter IDs to be equated.
model_equate <- function(model, params) {
  if (!inherits(model, "sinba_model")) {
    stop("model_equate: `model` must be an object of class \"sinba_model\".")
  }
  m <- model$model
  if (length(params) < 2) {
    stop("model_equate: `params` should at least have two elements.")
  }

  p <- params[1]
  for (cc in seq_len(ncol(m))) {
    for (r in seq_len(nrow(m))) {
      if (m[r, cc] %in% params) {
        m[r, cc] <- p
      }
    }
  }

  m <- format_model_matrix(m)
  model$name <- "user defined"
  model$model <- m
  return(model)
}

#' @export
#' @title Remove a Parameter From a Model
#'
#' @description
#' `model_drop()` removes one o more parameters
#' from a model.
#'
#' @param model A model definition as build with `new_model()`.
#' @param params A vector with the parameter IDs to be removed.
model_drop <- function(model, params) {
  if (!inherits(model, "sinba_model")) {
    stop("model_equate: `model` must be an object of class \"sinba_model\".")
  }
  m <- model$model
  if (length(params) < 1) {
    return(model)
  }

  for (cc in seq_len(ncol(m))) {
    for (r in seq_len(nrow(m))) {
      if (m[r, cc] %in% params) {
        m[r, cc] <- 0
      }
    }
  }

  m <- format_model_matrix(m)
  model$name <- "user defined"
  model$model <- m
  return(model)
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

# format_model_matrix checks the parameter labels
# in a model matrix
format_model_matrix <- function(model) {
  param <- list()
  states <- 0
  for (i in seq_len(nrow(model))) {
    for (j in seq_len(ncol(model))) {
      if (i == j) {
        model[i, j] <- 0
        next
      }
      if (model[i, j] == 0) {
        next
      }
      x <- param[model[i, j]][[1]]
      if (is.null(x)) {
        states <- states + 1
        x <- states
        param[[model[i, j]]] <- states
      }
      model[i, j] <- x
    }
  }
  return(model)
}

# collapse_model removes parameters from unobserved states.
collapse_model <- function(obs) {
  model <- model_matrix("ARD")
  for (i in seq_len(length(obs))) {
    if (obs[i]) {
      next
    }
    model[i, ] <- 0
    model[, i] <- 0
  }
  model <- format_model_matrix(model)
  if (max(model) == 0) {
    stop("collapse_model: model cannot be collapsed")
  }
  return(model)
}
