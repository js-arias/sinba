#' @export
#' @title
#' Maximum Likelihood Estimation of the Sinba Model in a Single Trait
#'
#' @description
#' `fit_sinba_single()` searches for the maximum likelihood estimate
#' with the Sinba model
#' for a single trait.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   A second should contain the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   The model must be defined for a single trait.
#' @param ev_prob set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba_single <- function(
    tree, data, model = NULL,
    ev_prob = 1,
    opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba_single: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(model)) {
    model <- new_model(traits = 1)
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_sinba_single: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits != 1) {
    stop("fit_sinba_single: `model` must be for a single trait.")
  }
  mQ <- model$model
  k <- max(mQ) + 1

  et <- encode_traits(t, data, 1)
  cond <- set_conditionals(t, et, model)

  root <- c("0", "1")

  ev_prob <- log(ev_prob)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    max_len <- t$age[yn]
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[2:length(p)] > 1000)) {
        return(Inf)
      }
      if (p[1] > max_len) {
        return(Inf)
      }
      n <- get_node_by_len(t, p[1], yn)
      if (n <= 0) {
        return(Inf)
      }
      birth <- list(
        node = n,
        age = t$br_len[n] + p[1] - t$age[n]
      )
      Q <- from_model_to_Q(mQ, p[2:length(p)])
      lk <- sinba_single_like(
        t, Q, model, birth, xt, cond,
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

  res <- list()
  res[[1]] <- list(objective = Inf)
  for (r in sample(seq_len(length(root)))) {
    youngest <- youngest_birth_event(t, et, root[r])
    yn <- youngest[[1]]
    fn <- like_func(yn, r)
    par <- c(runif(1, max = t$age[yn]), runif(max(mQ)))
    rr <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    if (rr$objective < res[[1]]$objective) {
      rr$root <- r
      rr$yn <- yn
      res <- list()
      res[[1]] <- rr
    } else if (rr$objective == res[[1]]$objective) {
      rr$root <- r
      rr$yn <- yn
      res[[length(res) + 1]] <- rr
    }
  }

  to_ret <- list()
  for (i in seq_len(length(res))) {
    rr <- res[[i]]
    q <- from_model_to_Q(mQ, rr$solution[2:length(rr$solution)])
    q <- normalize_Q(q)
    n <- get_node_by_len(t, rr$solution[1], rr$yn)
    birth <- list(
      node = n,
      age = t$br_len[n] + rr$solution[1] - t$age[n]
    )
    root_state <- root[rr$root]
    obj <- list(
      logLik = -rr$objective,
      k = k,
      model = model,
      Q = q,
      birth = birth,
      root = root_state,
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba_single"
    to_ret[[length(to_ret) + 1]] <- obj
  }
  return(to_ret)
}

#' @export
#' @title Extract Log-Likelihood From a "fit_sinba_single" Object
#'
#' @description
#' This method implements the `logLik` method
#' on a "fit_sinba_single" object.
#'
#' @param object An object of type "fit_sinba_single".
#' @param ... Additional arguments are unused.
logLik.fit_sinba_single <- function(object, ...) {
  l <- object$logLik
  attr(l, "df") <- object$k
  attr(l, "nobs") <- length(object$tree$tip.label)
  class(l) <- "logLik"
  return(l)
}


#' @export
#' @title Basic Print For a "fit_sinba_single" Object
#'
#' @description
#' This method implements the `print` method
#' on a `fit_sinba_single` object.
#'
#' @param x An object of type "fit_sinba_single".
#' @param digits The number of digits for decimal output.
#' @param ... Additional arguments are unused.
print.fit_sinba_single <- function(x, digits = 6, ...) {
  cat("Single Sinba: Fit\n")

  states <- x$model$states
  mm <- x$model$model
  rownames(mm) <- states
  colnames(mm) <- states
  cat("Model:\n")
  print(mm)
  cat(paste("Free parameters = ", x$k, ".\n", sep = ""))

  aic <- 2 * x$k - 2 * x$logLik
  aicc <- aic + (2 * x$k * x$k + 2 * x$k) / (length(x$tree$tip.label) - x$k - 1)
  fit <- c(x$logLik, aic, aicc)
  names(fit) <- c("logLik", "AIC", "AICc")
  print(fit)

  cat("Root state: ", x$root, "\n", sep = "")

  cat("Birth event:\n")
  b <- x$birth
  cat(paste("- Node ", b$node, " time ", round(b$age, digits), "\n", sep = ""))
  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- states
  colnames(Q) <- states
  print(Q)
}

#' @export
#' @title
#' Estimate the Likelihood for a Given Sinba Transition Matrix in a Single Trait
#'
#' @description
#' `fixed_sinba_single()` calculates the likelihood
#' for a given transition matrix
#' for a single trait
#' under the Sinba model.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   A second should contain the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param rate_mat Rate matrix for the trait.
#' @param birth A list with the field `node`
#'   indicated the birth of the trait,
#'   and `age` the time from the start of the edge
#'   in which the event happens.
#' @param model A model build with `new_model()`
#'   or `new_hidden_model()`.
#' @param ev_prob set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
fixed_sinba_single <- function(
    tree, data, rate_mat, birth,
    model = NULL,
    ev_prob = 1) {
  if (!inherits(tree, "phylo")) {
    stop("fixed_sinba_single: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(model)) {
    model <- new_model(traits = 1)
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fixed_sinba_single: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits != 1) {
    stop("fixed_sinba_single: `model` must be for a single trait.")
  }
  if (is.null(rate_mat)) {
    stop("fixed_sinba_single: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fixed_sinba_single: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop(
      "fixed_sinba_single: `rate_mat` must have the same size as the `model`"
    )
  }

  et <- encode_traits(t, data, 1)
  cond <- set_conditionals(t, et, model)

  root <- c("0", "1")

  ev_prob <- log(ev_prob)

  if (is.null(birth)) {
    stop("fixed_sinba_single: undefined `birth` event")
  }
  if (is.null(birth$node)) {
    stop("fixed_sinba_singe: expecting field `node` of `birth`")
  }
  if (is.null(birth$age)) {
    stop("fixed_sinba_singe: expecting field `age` of `birth`")
  }
  if (birth$age > t$br_len[birth$node]) {
    stop("fixed_sinba_singe: invalid value for field `age` of `birth`")
  }

  res <- list()
  res[[1]] <- list(objective = Inf)
  for (r in sample(seq_len(length(root)))) {
    youngest <- youngest_birth_event(t, et, root[r])
    yn <- youngest[[1]]
    if (!is_valid_birth(t, birth$node, yn)) {
      next
    }
    xt <- tree_to_cpp(t)
    lk <- sinba_single_like(
      t, rate_mat, model, birth, xt, cond,
      r, ev_prob
    )
    if (lk < res[[1]]$objective) {
      rr <- list(
        objective = lk,
        root = r
      )
      res <- list()
      res[[1]] <- rr
    } else if (lk == res[[1]]$objective) {
      rr <- list(
        objective = lk,
        root = r
      )
      res[[length(res) + 1]] <- rr
    }
  }

  # if no birth sequence is compatible with the birth event
  # the likelihood is 0
  if (is.infinite(res[[1]]$objective)) {
    obj <- list(
      logLik = -Inf,
      k = 0,
      model = model,
      Q = normalize_Q(rate_mat),
      birth = birth,
      root = "NA",
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba_single"
    return(obj)
  }

  to_ret <- list()
  for (i in seq_len(length(res))) {
    rr <- res[[i]]
    root_state <- root[rr$root]
    obj <- list(
      logLik = rr$objective,
      k = 0,
      model = model,
      Q = normalize_Q(rate_mat),
      birth = birth,
      root = root_state,
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba_single"
    to_ret[[length(to_ret) + 1]] <- obj
  }
  return(to_ret)
}

# sinba_single calculates the likelihood of a single trait
# under the sinba model.
sinba_single_like <- function(
    t, Q, model, birth,
    xt, cond, root, ev_prob) {
  # make sure the Q matrix is valid
  Q <- normalize_Q(Q)
  root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))

  st <- as.integer(active_status(t, birth$node, length(t$parent) + 1))
  l <- full_sinba_conditionals(
    xt$parent, xt$nodes, st, xt$branch,
    cond,
    birth$age, 0,
    root_Q, Q, root_Q
  )

  # takes into account the hidden states
  root_states <- c("0", "1")
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
