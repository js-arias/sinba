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
#' @param traits The number of traits in the model.
#'   By default it is two.
#'   It can also be a single trait.
new_model <- function(model = "", traits = 2) {
  if (traits != 2) {
    states <- c("x0", "x1")
    observed <- list()
    observed[[states[1]]] <- "0"
    observed[[states[2]]] <- "1"
    m <- matrix(0, nrow = 2, ncol = 2)
    m[1, 2] <- 1
    m[2, 1] <- 2
    obj <- list(
      name = "single trait",
      model = m,
      traits = 1,
      states = states,
      observed = observed
    )
    class(obj) <- "sinba_model"
    return(obj)
  }
  names <- c("IND", "CORR", "ER", "ER2", "ERs", "SYM", "sCORR", "x", "y", "CMK")
  if (!(model %in% names)) {
    model <- "IND"
  }
  states <- c("x0,y0", "x0,y1", "x1,y0", "x1,y1")
  observed <- list()
  observed[[states[1]]] <- "00"
  observed[[states[2]]] <- "01"
  observed[[states[3]]] <- "10"
  observed[[states[4]]] <- "11"
  m <- model_matrix(model)
  obj <- list(
    name = model,
    model = m,
    traits = 2,
    states = states,
    observed = observed
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
#' @title Create a New Model With Hidden States
#'
#' @description
#' `new_hidden_model()` creates a new model matrix
#' with hidden states.
#' The hidden states are defined in a list
#' with three elements:
#' `trait` to indicates the trait
#' (either "x", the first trait,
#' or "y", the second trait);
#' `state` to indicate the observed state
#' (either 0 or 1);
#' and `hidden` with a vector of the labels
#' for the hidden states
#' (e.g. "a", "b", "c").
#' By default,
#' the model is the maximal independent model
#' for the given states.
#'
#' @param states A list with one or more elements,
#'   each element containing a trait state
#'   and its hidden states.
#' @param traits The number of traits in the model.
#'   By default it is two.
#'   It can also be a single trait.
#'
#' @return A model matrix with hidden states.
new_hidden_model <- function(states, traits = 2) {
  if (is.null(states)) {
    return(new_model(), traits)
  }

  if (traits != 2) {
    return(single_hidden_model(states))
  }

  # states with hidden states
  hs <- c(FALSE, FALSE, FALSE, FALSE)
  x_states <- c()
  y_states <- c()
  for (e in states) {
    if (is.null(e$trait)) {
      next
    }
    if (e$trait != "x" && e$trait != "y") {
      next
    }
    if (is.null(e$state)) {
      next
    }
    if (e$state != 0 && e$state != 1) {
      next
    }
    if (length(e$hidden) <= 1) {
      next
    }
    obs <- 0
    if (e$trait == "x") {
      obs <- 1
      if (e$state == 1) {
        obs <- 2
      }
    } else {
      obs <- 3
      if (e$state == 1) {
        obs <- 4
      }
    }
    if (hs[obs]) {
      next
    }
    for (i in seq_len(length(e$hidden))) {
      s <- sprintf("%s%d%s", e$trait, e$state, e$hidden[i])
      if (e$trait == "x") {
        x_states <- c(x_states, s)
      } else {
        y_states <- c(y_states, s)
      }
    }
    hs[obs] <- TRUE
  }
  if (!any(hs)) {
    return(new_model())
  }
  for (i in seq_len(length(hs))) {
    if (hs[i]) {
      next
    }
    if (i == 1) {
      x_states <- c(x_states, "x0")
    }
    if (i == 2) {
      x_states <- c(x_states, "x1")
    }
    if (i == 3) {
      y_states <- c(y_states, "y0")
    }
    if (i == 4) {
      y_states <- c(y_states, "y1")
    }
  }
  x_states <- sort(x_states)
  y_states <- sort(y_states)

  # combine individual states
  observed <- list()
  states <- c()
  for (i in seq_len(length(x_states))) {
    for (j in seq_len(length(y_states))) {
      s <- sprintf("%s,%s", x_states[i], y_states[j])
      o <- ""
      if (startsWith(x_states[i], "x0")) {
        o <- "00"
        if (startsWith(y_states[j], "y1")) {
          o <- "01"
        }
      } else {
        o <- "10"
        if (startsWith(y_states[j], "y1")) {
          o <- "11"
        }
      }
      observed[[s]] <- o
      states <- c(states, s)
    }
  }
  states <- sort(states)

  # build the matrix
  m <- matrix(0, nrow = length(states), ncol = length(states))
  rownames(m) <- states
  colnames(m) <- states
  p <- 1

  # x trait
  for (i in seq_len(length(x_states))) {
    for (j in seq_len(length(x_states))) {
      if (i == j) {
        next
      }
      for (k in seq_len(length(y_states))) {
        from <- sprintf("%s,%s", x_states[i], y_states[k])
        to <- sprintf("%s,%s", x_states[j], y_states[k])
        m[from, to] <- p
      }
      p <- p + 1
    }
  }

  # y trait
  for (i in seq_len(length(y_states))) {
    for (j in seq_len(length(y_states))) {
      if (i == j) {
        next
      }
      for (k in seq_len(length(x_states))) {
        from <- sprintf("%s,%s", x_states[k], y_states[i])
        to <- sprintf("%s,%s", x_states[k], y_states[j])
        m[from, to] <- p
      }
      p <- p + 1
    }
  }

  rownames(m) <- NULL
  colnames(m) <- NULL
  m <- format_model_matrix(m)
  obj <- list(
    name = "hidden",
    model = m,
    traits = 2,
    states = states,
    observed = observed
  )
  class(obj) <- "sinba_model"
  return(obj)
}

# single_hidden_model creates a model
# for a single trait
# with hidden states.
single_hidden_model <- function(states) {
  # states with hidden states
  hs <- c(FALSE, FALSE)
  x_states <- c()
  for (e in states) {
    if (is.null(e$trait)) {
      next
    }
    if (e$trait != "x") {
      next
    }
    if (is.null(e$state)) {
      next
    }
    if (e$state != 0 && e$state != 1) {
      next
    }
    if (length(e$hidden) <= 1) {
      next
    }
    obs <- e$state + 1
    if (hs[obs]) {
      next
    }
    for (i in seq_len(length(e$hidden))) {
      s <- sprintf("%s%d%s", e$trait, e$state, e$hidden[i])
      x_states <- c(x_states, s)
    }
    hs[obs] <- TRUE
  }

  if (!any(hs)) {
    return(new_model(traits = 1))
  }
  for (i in seq_len(length(hs))) {
    if (hs[i]) {
      next
    }
    if (i == 1) {
      x_states <- c(x_states, "x0")
    }
    if (i == 2) {
      x_states <- c(x_states, "x1")
    }
  }
  x_states <- sort(x_states)

  # combine individual states
  observed <- list()
  for (i in seq_len(length(x_states))) {
    o <- ""
    if (startsWith(x_states[i], "x0")) {
      o <- "0"
    } else {
      o <- "1"
    }
    observed[[x_states[i]]] <- o
  }

  # build the matrix
  m <- matrix(0, nrow = length(x_states), ncol = length(x_states))
  p <- 1

  for (i in seq_len(length(x_states))) {
    for (j in seq_len(length(x_states))) {
      if (i == j) {
        next
      }
      m[i, j] <- p
      p <- p + 1
    }
  }

  m <- format_model_matrix(m)
  obj <- list(
    name = "single hidden",
    model = m,
    traits = 1,
    states = x_states,
    observed = observed
  )
  class(obj) <- "sinba_model"
  return(obj)
}

#' @export
#' @title Create a New Model With Hidden Rates
#'
#' @description
#' `new_rates_model()` creates a new model matrix
#' with hidden rates.
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
#' @param rates The labels to identify the rates.
#' @param traits The number of traits in the model.
#'   By default it is two.
#'   It can also be a single trait.
new_rates_model <- function(model = "", rates = NULL, traits = 2) {
  if (length(rates) < 2) {
    return(new_model(model, traits))
  }

  names <- c("IND", "CORR", "ER", "ER2", "ERs", "SYM", "sCORR", "x", "y", "CMK")
  if (!(model %in% names)) {
    model <- "IND"
  }
  states <- rep(0, 4 * length(rates))
  base <- model_matrix(model)
  if (traits != 2) {
    states <- rep(0, 2 * length(rates))
    base <- matrix(0, nrow = 2, ncol = 2)
    base[1, 2] <- 1
    base[2, 1] <- 2
    if (model == "ER") {
      base[2, 1] <- 1
    }
    traits <- 1
  }
  m <- matrix(0, nrow = length(states), ncol = length(states))

  # state transitions
  for (i in seq_len(length(rates))) {
    mx <- max(m)
    offset <- (i - 1) * nrow(base)
    for (j in seq_len(nrow(base))) {
      for (k in seq_len(nrow(base))) {
        v <- base[j, k]
        if (v > 0) {
          v <- v + mx
        }
        m[j + offset, k + offset] <- v
      }
    }
  }

  # rate transitions
  rt <- max(m) + 1
  for (i in seq_len(length(rates))) {
    off_i <- (i - 1) * nrow(base)
    for (j in seq_len(length(rates))) {
      if (i == j) {
        next
      }
      off_j <- (j - 1) * nrow(base)
      for (k in seq_len(nrow(base))) {
        m[k + off_i, k + off_j] <- rt
      }
    }
  }

  obs <- c("00", "01", "10", "11")
  sts <- c("x0,y0", "x0,y1", "x1,y0", "x1,y1")
  if (traits != 2) {
    obs <- c("0", "1")
    sts <- c("x0", "x1")
  }
  states <- c()
  observed <- list()
  for (i in seq_len(length(rates))) {
    for (j in seq_len(length(obs))) {
      s <- sprintf("%s[%s]", sts[j], rates[i])
      states <- c(states, s)
      observed[[s]] <- obs[j]
    }
  }

  m <- format_model_matrix(m)
  obj <- list(
    name = "hidden rates",
    model = m,
    traits = traits,
    states = states,
    observed = observed
  )
  class(obj) <- "sinba_model"
  return(obj)
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
  m <- format_model_matrix(m)
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

# build_semi_active_Q returns the semi-active Q
# using a model,
# the birth scenario,
# and the Q matrix,
build_semi_active_Q <- function(model, sc, Q) {
  sm <- matrix(0, nrow = 4, ncol = 4)
  if (sc == "12") {
    # second trait active: 00 <-> 01
    sm[1, 1] <- 1
    sm[1, 2] <- 1
    sm[2, 1] <- 1
    sm[2, 2] <- 1
  }
  if (sc == "13") {
    # first trait active: 00 <-> 10
    sm[1, 1] <- 1
    sm[1, 3] <- 1
    sm[3, 1] <- 1
    sm[3, 3] <- 1
  }
  if (sc == "24") {
    # first trait active: 01 <-> 11
    sm[2, 2] <- 1
    sm[2, 4] <- 1
    sm[4, 2] <- 1
    sm[4, 4] <- 1
  }
  if (sc == "34") {
    # second trait active: 10 <-> 11
    sm[3, 3] <- 1
    sm[3, 4] <- 1
    sm[4, 3] <- 1
    sm[4, 4] <- 1
  }
  rownames(sm) <- c("00", "01", "10", "11")
  colnames(sm) <- c("00", "01", "10", "11")

  m <- model$model
  semi <- matrix(0, nrow = nrow(m), ncol = ncol(m))
  for (i in seq_len(nrow(m))) {
    from <- model$observed[[model$states[i]]]
    for (j in seq_len(ncol(m))) {
      if (i == j) {
        next
      }
      if (m[i, j] == 0) {
        next
      }
      to <- model$observed[[model$states[j]]]
      if (sm[from, to] == 0) {
        next
      }
      semi[i, j] <- Q[i, j]
    }
  }
  return(semi)
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

# reduce_matrix creates a new matrix without zero cols-rows.
reduce_matrix <- function(m) {
  n <- colnames(m)
  states <- c()
  for (i in seq_len(length(n))) {
    if (any(m[i, ] != 0)) {
      states <- c(states, n[i])
    }
  }
  r <- matrix(0, ncol = length(states), nrow = length(states))
  colnames(r) <- states
  rownames(r) <- states
  for (i in seq_len(length(states))) {
    from <- states[i]
    for (j in seq_len(length(states))) {
      if (i == j) {
        next
      }
      to <- states[j]
      r[from, to] <- m[from, to]
    }
  }
  return(r)
}
