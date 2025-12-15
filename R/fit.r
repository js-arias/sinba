#' @import nloptr

#' @export
#' @title
#' Maximum Likelihood Estimation of the Sinba Model
#'
#' @description
#' `fit_sinba()` searches for the maximum likelihood estimate
#' with the Sinba model for two traits.
#' It estimate the values of the parameters for the Q matrices
#' as well as the birth events.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param ev_prob Set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param single_birth If true,
#'   it will set the birth of both traits in the same location
#'   (reducing one parameter).
#' @param root Set the root state.
#'   By default,
#'   the root state will be optimized as a parameter.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba <- function(
    tree, data, model = NULL,
    single_birth = FALSE,
    ev_prob = 1,
    root = "",
    opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop("fit_sinba: `model` must be an object of class \"sinba_model\".")
  }
  mQ <- model$model
  k <- max(mQ) + 2
  if (single_birth) {
    k <- max(mQ) + 1
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  root_states <- c("00", "01", "10", "11")
  if (root != "") {
    if (!(root %in% root_states)) {
      stop("fit_sinba: invalid root state.")
    }
  }

  ev_prob <- 2 * log(ev_prob)

  # separate birth parameters
  # from transition (traditional) parameters
  transition_start <- 3
  if (single_birth) {
    transition_start <- 2
  }

  # closure for the likelihood function
  like_func <- function(yn, r) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[transition_start:length(p)] > 1000)) {
        return(Inf)
      }

      births <- list()
      for (i in 1:2) {
        j <- i
        if (single_birth) {
          j <- 1
        }

        n <- get_node_by_len(t, p[j], yn[j])
        if (n <= 0) {
          return(Inf)
        }
        b <- list(
          node = n,
          age = t$br_len[n] + p[j] - t$age[n]
        )
        births[[i]] <- b
      }

      Q <- from_model_to_Q(mQ, p[transition_start:length(p)])
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        r, ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      xtol_abs = rep(v, k),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (root != "") {
      if (root != root_states[r]) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    # check all possible birth events
    ev <- birth_events(t, youngest)
    for (i in sample(seq_len(length(ev)))) {
      e <- ev[[i]]
      fn <- like_func(e, r)
      par <- c(
        runif(1, max = t$age[e[1]]),
        runif(1, max = t$age[e[2]]), runif(k - 2)
      )
      if (single_birth) {
        max_age <- t$age[e[1]]
        if (max_age > t$age[e[2]]) {
          max_age <- t$age[e[2]]
        }
        par <- c(runif(1, max = max_age), runif(k - 1))
      }
      rr <- nloptr::nloptr(
        x0 = par,
        eval_f = fn,
        opts = opts
      )
      if (rr$objective < res$objective) {
        rr$root <- r
        rr$ev_nodes <- e
        res <- rr
      }
    }
  }

  q <- from_model_to_Q(mQ, res$solution[transition_start:length(res$solution)])
  q <- normalize_Q(q)
  births <- list()
  for (i in 1:2) {
    j <- i
    if (single_birth) {
      j <- 1
    }
    n <- get_node_by_len(t, res$solution[j], res$ev_node[j])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[j] - t$age[n]
    )
    births[[i]] <- b
  }
  root_state <- root_states[res$root]

  # retrieve the scenario
  sc <- NULL
  if (!single_birth) {
    v <- scenario(res$root, 1)
    if (res$solution[2] < res$solution[1]) {
      # second trait is the oldest one
      v <- scenario(res$root, 2)
    }
    sc <- v
  }

  # if we have to calculate the root,
  # we add a parameter.
  if (root == "") {
    k <- k + 1
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenario = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

#' @export
#' @title Extract Log-Likelihood From a "fit_sinba" Object
#'
#' @description
#' This method implements the `logLik` method
#' on a "fit_sinba" object.
#'
#' @param object An object of type "fit_sinba".
#' @param ... Additional arguments are unused.
logLik.fit_sinba <- function(object, ...) {
  l <- object$logLik
  attr(l, "df") <- object$k
  attr(l, "nobs") <- 2 * length(object$tree$tip.label)
  class(l) <- "logLik"
  return(l)
}

#' @export
#' @title Basic Print For a "fit_sinba" Object
#'
#' @description
#' This method implements the `print` method
#' on a `fit_sinba` object.
#'
#' @param x An of type "fit_sinba".
#' @param digits The number of digits for decimal output.
#' @param ... Additional arguments are unused.
print.fit_sinba <- function(x, digits = 6, ...) {
  cat("Sinba: Fit\n")

  states <- x$model$states
  mm <- x$model$model
  rownames(mm) <- states
  colnames(mm) <- states
  cat("Model:\n")
  print(mm)
  cat(paste("Free parameters = ", x$k, ".\n", sep = ""))

  aic <- 2 * x$k - 2 * x$logLik
  aicc <- aic + (2 * x$k * x$k + 2 * x$k) /
    (2 * length(x$tree$tip.label) - x$k - 1)
  fit <- c(x$logLik, aic, aicc)
  names(fit) <- c("logLik", "AIC", "AICc")
  print(fit)

  cat("Root state: ", x$root, "\n", sep = "")

  cat("Birth events:\n")
  n <- colnames(x$data)
  b1 <- x$births[[1]]
  cat(paste("- Trait ", n[2], " Node ",
    b1$node, " time ", round(b1$age, digits), "\n",
    sep = ""
  ))
  b2 <- x$births[[2]]
  cat(paste("- Trait ", n[3], " Node ",
    b2$node, " time ", round(b2$age, digits), "\n",
    sep = ""
  ))
  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- states
  colnames(Q) <- states
  print(Q)
  age1 <- phylo_node_age(x$tree, b1$node)
  age2 <- phylo_node_age(x$tree, b2$node)
  if (b1$node == b2$node) {
    age1 <- b1$age
    age2 <- b2$age
  }
  if ((age1 != age2) && (length(x$scenarios) > 0)) {
    cat("Semi-active process:\n")
    for (i in seq_len(length(x$scenarios))) {
      sc <- x$scenarios[i]
      sQ <- build_semi_active_Q(x$model, sc, x$Q)
      if (sc == "12") {
        cat("second trait active: 00 <-> 01\n")
      } else if (sc == "13") {
        cat("first trait active: 00 <-> 10\n")
      } else if (sc == "24") {
        cat("first trait active: 01 <-> 11\n")
      } else if (sc == "34") {
        cat("second trait active: 10 <-> 11\n")
      }
      rownames(sQ) <- states
      colnames(sQ) <- states
      sQ <- reduce_matrix(sQ)
      print(normalize_Q(sQ))
    }
  }
}

#' @export
#' @title
#' Maximum Likelihood Estimation of a Transition Matrix With Fixed Births
#'
#' @description
#' `fit_fixed_births()` searches for the maximum likelihood estimate
#' of the transition matrix
#' with the Sinba model with two traits
#' with a fixed birth events.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param births A list with two element,
#'   each element being a list that define a birth event,
#'   with a fields `node` indicating the birth of the trait
#'   and `age` indicating the time from the start of the edge
#'   in which the event happens.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param ev_prob Set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param root Set the root state.
#'   By default,
#'   the root state will be optimized as a parameter.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_births <- function(
    tree, data, births, model = NULL,
    ev_prob = 1,
    root = "",
    opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_births: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_fixed_births: `model` must be an object of class \"sinba_model\"."
    )
  }

  mQ <- model$model
  k <- max(mQ)

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  root_states <- c("00", "01", "10", "11")
  if (root != "") {
    if (!(root %in% root_states)) {
      stop("fit_fixed_births: invalid root state.")
    }
  }

  ev_prob <- 2 * log(ev_prob)

  if (is.null(births)) {
    stop("fit_fixed_births: undefined `births` events")
  }
  if (length(births) < 2) {
    stop("fit_fixed_births: two events required in `births`")
  }
  for (i in 1:2) {
    e <- births[[i]]
    if (is.null(e$node)) {
      stop(sprintf(
        "fit_fixed_births: expecting field `node` in `births` event %d", i
      ))
    }
    if (is.null(e$age)) {
      stop(sprintf(
        "fit_fixed_births: expecting field `age` in `births` event %d", i
      ))
    }
    if (e$age > t$br_len[e$node]) {
      stop(sprintf(
        "fit_fixed_births: invalid value for field `age` in `births` event %d",
        i
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
        k = k,
        model = model,
        Q = matrix(nrow = 4, ncol = 4),
        births = births,
        root = "NA",
        data = data,
        tree = tree
      )
      class(obj) <- "fit_sinba"
      return(obj)
    }
  }

  # closure for the likelihood function
  like_func <- function(r) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p > 1000)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p)
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        r, ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      xtol_abs = rep(v, k),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (root != "") {
      if (root != root_states[r]) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    if (!is_valid_birth(t, births[[1]]$node, youngest[[1]])) {
      next
    }
    if (!is_valid_birth(t, births[[2]]$node, youngest[[2]])) {
      next
    }
    fn <- like_func(r)
    par <- c(runif(k))
    rr <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    if (rr$objective < res$objective) {
      rr$root <- r
      res <- rr
    }
  }

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (is.infinite(res$objective)) {
    obj <- list(
      logLik = -Inf,
      k = k,
      model = model,
      Q = matrix(nrow = 4, ncol = 4),
      births = births,
      root = "NA",
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba"
    return(obj)
  }

  q <- from_model_to_Q(mQ, res$solution)
  q <- normalize_Q(q)
  root_state <- root_states[res$root]

  # retrieve the scenario
  age1 <- phylo_node_age(tree, births[[1]]$node)
  age2 <- phylo_node_age(tree, births[[1]]$node)
  if (births[[1]]$node == births[[2]]$node) {
    age1 <- b1$age
    age2 <- b2$age
  }
  sc <- scenario(res$root, 1)
  if (age2 < age1) {
    # second trait is the oldest one
    sc <- scenario(res$root, 2)
  }

  # if we have to calculate the root,
  # we add a parameter.
  if (root == "") {
    k <- k + 1
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenario = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

#' @export
#' @title
#' Maximum Likelihood Estimation of Births With Fixed Transition Matrix
#'
#' @description
#' `fit_fixed_matrix()` searches for the maximum likelihood estimate
#' of birth events of two traits
#' given a transition matrix
#' under the Sinba model.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param rate_mat Rate matrix for the traits with the full process.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param ev_prob Set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param root Set the root state.
#'   By default,
#'   the root state will be optimized as a parameter.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_matrix <- function(
    tree, data, rate_mat,
    model = NULL,
    ev_prob = 1,
    root = "",
    opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_matrix: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_fixed_matrix: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (is.null(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop("fit_fixed_matrix: `rate_mat` must have the same size as the `model`")
  }

  k <- 2

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  root_states <- c("00", "01", "10", "11")
  if (root != "") {
    if (!(root %in% root_states)) {
      stop("fit_fixed_matrix: invalid root state.")
    }
  }

  ev_prob <- 2 * log(ev_prob)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }

      births <- list()
      for (i in 1:2) {
        n <- get_node_by_len(t, p[i], yn[i])
        if (n <= 0) {
          return(Inf)
        }
        b <- list(
          node = n,
          age = t$br_len[n] + p[i] - t$age[n]
        )
        births[[i]] <- b
      }

      lk <- sinba_like(
        t, rate_mat, model, births, xt, cond,
        r, ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      xtol_abs = rep(v, k),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (root != "") {
      if (root != root_states[r]) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    # check all possible birth events
    ev <- birth_events(t, youngest)
    for (i in sample(seq_len(length(ev)))) {
      e <- ev[[i]]
      fn <- like_func(e, r)
      par <- c(
        runif(1, max = t$age[e[1]]),
        runif(1, max = t$age[e[2]])
      )
      rr <- nloptr::nloptr(
        x0 = par,
        eval_f = fn,
        opts = opts
      )
      if (rr$objective < res$objective) {
        rr$root <- r
        rr$ev_nodes <- e
        res <- rr
      }
    }
  }

  q <- normalize_Q(rate_mat)
  births <- list()
  for (i in 1:2) {
    j <- i
    n <- get_node_by_len(t, res$solution[j], res$ev_node[j])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[j] - t$age[n]
    )
    births[[i]] <- b
  }
  root_state <- root_states[res$root]

  # retrieve the scenarios
  sc <- scenario(res$root, 1)
  if (res$solution[2] < res$solution[1]) {
    # second trait is the oldest one
    sc <- scenario(res$root, 2)
  }

  # if we have to calculate the root,
  # we add a parameter.
  if (root == "") {
    k <- k + 1
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenario = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

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
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param ev_prob Set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param root Set the root state.
#'   By default,
#'   the root state will be optimized as a parameter.
fixed_sinba <- function(
    tree, data, rate_mat, births,
    model = NULL,
    ev_prob = 1,
    root = "") {
  if (!inherits(tree, "phylo")) {
    stop("fixed_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop("fixed_sinba: `model` must be an object of class \"sinba_model\".")
  }
  if (is.null(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop("fixed_sinba: `rate_mat` must have the same size as the `model`")
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  root_states <- c("00", "01", "10", "11")
  if (root != "") {
    if (!(root %in% root_states)) {
      stop("fixed_sinba: invalid root state.")
    }
  }

  ev_prob <- 2 * log(ev_prob)

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
        k = 0,
        model = model,
        Q = normalize_Q(rate_mat),
        births = births,
        root = "NA",
        data = data,
        tree = tree
      )
      class(obj) <- "fit_sinba"
      return(obj)
    }
  }

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (root != "") {
      if (root != root_states[r]) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    if (!is_valid_birth(t, births[[1]]$node, youngest[[1]])) {
      next
    }
    if (!is_valid_birth(t, births[[2]]$node, youngest[[2]])) {
      next
    }
    xt <- tree_to_cpp(t)
    lk <- -sinba_like(
      t, rate_mat, model, births, xt, cond,
      r, ev_prob
    )
    if (lk < res$objective) {
      res <- list(
        objective = lk,
        root = r
      )
    }
  }

  k <- 0
  # if we have to calculate the root,
  # we add a parameter.
  if (root == "") {
    k <- 1
  }

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (is.infinite(res$objective)) {
    obj <- list(
      logLik = -Inf,
      k = k,
      model = model,
      Q = normalize_Q(rate_mat),
      births = births,
      root = "NA",
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba"
    return(obj)
  }

  q <- normalize_Q(rate_mat)
  root_state <- root_states[res$root]

  # retrieve the scenario
  age1 <- phylo_node_age(tree, births[[1]]$node)
  age2 <- phylo_node_age(tree, births[[1]]$node)
  if (births[[1]]$node == births[[2]]$node) {
    age1 <- b1$age
    age2 <- b2$age
  }
  sc <- scenario(res$root, 1)
  if (age2 < age1) {
    # second trait is the oldest one
    sc <- scenario(res$root, 2)
  }

  obj <- list(
    logLik = -res$objective,
    k = 0,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenario = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_like calculates the likelihood
# of the sinba model.
sinba_like <- function(
    t, Q, model, births, xt, cond, root, ev_prob) {
  Q <- normalize_Q(Q)

  root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))

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
    sc <- scenario(root, 1)
    semi <- build_semi_active_Q(model, sc, Q)
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
    sc <- scenario(root, 2)
    semi <- build_semi_active_Q(model, sc, Q)
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

  # takes into account the hidden states
  root_states <- c("00", "01", "10", "11")
  likes <- c()
  for (i in seq_len(ncol(l))) {
    obs <- model$observed[[model$states[i]]]
    if (obs == root_states[root]) {
      likes <- c(likes, l[t$root_id, i])
    }
  }
  scaled <- exp(likes - max(likes))
  fitz <- scaled / sum(scaled)
  like <- log(sum(fitz * scaled)) + max(likes)
  return(like)
}
