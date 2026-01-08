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
#' @param ev_prob Det the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param root Set the root state.
#'   By default,
#'   the root state will be optimized as a parameter.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba_single <- function(
    tree, data, model = NULL,
    ev_prob = 1,
    root = "",
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

  root_states <- c("0", "1")
  if (root != "") {
    if (!(root %in% root_states)) {
      stop("fit_sinba_single: invalid root state.")
    }
  }

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
    opts <- def_nloptr_opts(k)
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
    yn <- youngest[[1]]
    fn <- like_func(yn, r)
    par <- c(runif(1, max = t$age[yn]), runif(max(mQ)))
    rr <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    if (rr$objective < res$objective) {
      rr$root <- r
      rr$yn <- yn
      res <- rr
    }
  }

  q <- from_model_to_Q(mQ, res$solution[2:length(res$solution)])
  q <- normalize_Q(q)
  n <- get_node_by_len(t, res$solution[1], res$yn)
  birth <- list(
    node = n,
    age = t$br_len[n] + res$solution[1] - t$age[n]
  )
  root_state <- root_states[res$root]

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
    birth = birth,
    root = root_state,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba_single"
  return(obj)
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
#' @param ev_prob Set the probability of a birth event.
#'   By default is 1
#'   (i.e., we have observed different states in the traits).
#' @param root Set the root state.
#'   By default,
#'   the root state will be optimized as a parameter.
fixed_sinba_single <- function(
    tree, data, rate_mat, birth,
    model = NULL,
    ev_prob = 1,
    root = "") {
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

  root_states <- c("0", "1")
  if (root != "") {
    if (!(root %in% root_states)) {
      stop("fixed_sinba_single: invalid root state.")
    }
  }

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

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (root != "") {
      if (root != root_states[r]) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    yn <- youngest[[1]]
    if (!is_valid_birth(t, birth$node, yn)) {
      next
    }
    xt <- tree_to_cpp(t)
    lk <- sinba_single_like(
      t, rate_mat, model, birth, xt, cond,
      r, ev_prob
    )
    if (lk < res$objective) {
      res <- list(
        objective = lk,
        root = r
      )
    }
  }

  # if no birth sequence is compatible with the birth event
  # the likelihood is 0
  if (is.infinite(res$objective)) {
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

  root_state <- root_states[res$root]
  k <- 0
  # if we have to calculate the root,
  # we add a parameter.
  if (root == "") {
    k <- 1
  }

  obj <- list(
    logLik = res$objective,
    k = k,
    model = model,
    Q = normalize_Q(rate_mat),
    birth = birth,
    root = root_state,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba_single"
  return(obj)
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
