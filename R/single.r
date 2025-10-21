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
#' @param model A model build with `new_model()`
#'   or `new_hidden_model()`.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba_single <- function(
    tree, data, model = NULL,
    root = NULL, opts = NULL) {
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

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_sinba_single: invalid size for `root` vector.")
  }
  root <- root / sum(root)

  youngest <- youngest_birth_node(t, et, 1)
  ev_prob <- log(prob_birth(t))

  # closure for the likelihood function
  like_func <- function(yn) {
    obs_prior <- c(1, 1)
    y <- youngest[[1]]
    if (y[2] != yn) {
      if (!is_parent(t, yn, y[2])) {
        obs_prior[1] <- 0
      }
    }
    y <- youngest[[1]]
    if (y[1] != yn) {
      if (!is_parent(t, yn, y[1])) {
        obs_prior[2] <- 0
      }
    }
    root_prob <- set_root_prior_single(model, root, obs_prior)

    max_len <- t$age[yn]
    xt <- tree_to_cpp(t)
    log_root <- log(root_prob)

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
        log_root, ev_prob
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

  y <- youngest[[1]]
  res <- list(objective = Inf)
  if (any(y == t$root_id)) {
    # if at least one state born at root
    yn <- y[1]
    if (t$age[y[2]] > t$age[yn]) {
      yn <- y[2]
    }

    fn <- like_func(yn)
    par <- c(runif(1, max = t$age[yn]), runif(max(mQ)))
    res <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    res$yn <- yn
  } else {
    # for rare cases in which both states born in a node
    # different from the root
    # for example ((a,(a, a)),(b,(b, b)))
    for (yn in y) {
      fn <- like_func(yn)
      par <- c(runif(1, max = t$age[yn]), runif(max(mQ)))
      r <- nloptr::nloptr(
        x0 = par,
        eval_f = fn,
        opts = opts
      )
      if (r$objective < res$objective) {
        res <- r
        res$yn <- yn
      }
    }
  }

  q <- from_model_to_Q(mQ, res$solution[2:length(res$solution)])
  q <- normalize_Q(q)
  n <- get_node_by_len(t, res$solution[1], res$yn)
  birth <- list(
    node = n,
    age = t$br_len[n] + res$solution[1] - t$age[n]
  )
  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    birth = birth,
    root_prior = root,
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

  cat("Birth event:\n")
  b <- x$birth
  cat(paste("- Node ", b$node, " time ", round(b$age, digits), "\n", sep = ""))
  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- states
  colnames(Q) <- states
  print(Q)
  cat("Root prior:\n")
  root <- x$root_prior
  names(root) <- states
  print(root)
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
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
fixed_sinba_single <- function(
    tree, data, rate_mat, birth,
    model = NULL,
    root = NULL) {
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

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fixed_sinba_single: invalid size for `root` vector.")
  }
  root <- root / sum(root)

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

  youngest <- youngest_birth_node(t, et, 1)
  y <- youngest[[1]]
  obs_prior <- c(1, 1)
  for (yn in seq_len(length(y))) {
    if (y[yn] == birth$node) {
      next
    }
    if (is_parent(t, y[yn], birth$node)) {
      next
    }
    obs_prior[yn] <- 0
  }
  ev_prob <- log(prob_birth(t))

  root_prob <- set_root_prior_single(model, root, obs_prior)

  if (all(root_prob == 0)) {
    obj <- list(
      logLik = -Inf,
      model = model,
      Q = normalize_Q(rate_mat),
      birth = birth,
      root_prior = root,
      data = data,
      tree = tree
    )
    class(obj) <- "fixed_sinba_single"
    return(obj)
  }

  xt <- tree_to_cpp(t)
  lk <- sinba_single_like(
    t, rate_mat, model, birth, xt, cond,
    log(root_prob), ev_prob
  )

  obj <- list(
    logLik = lk,
    Q = normalize_Q(rate_mat),
    birth = birth,
    root_prior = root,
    data = data,
    tree = tree
  )
  class(obj) <- "fixed_sinba_single"
  return(obj)
}

#' @export
#' @title Basic Print For a "fixed_sinba_single" Object
#'
#' @description
#' This method implements the `print` method
#' on a `fixed_sinba_single` object.
#'
#' @param x An object of type "fixed_sinba_single".
#' @param digits The number of digits for decimal output.
#' @param ... Additional arguments are unused.
print.fixed_sinba_single <- function(x, digits = 6, ...) {
  cat("Single Sinba: Fixed Rate Matrix\n")
  cat(paste("Log-Likelihood = ", round(x$logLik, digits), "\n", sep = ""))
  cat("Birth event:\n")
  b <- x$birth
  cat(paste("- Node ", b$node, " time ", round(b$age, digits), "\n", sep = ""))
  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- c(0, 1)
  colnames(Q) <- c(0, 1)
  print(Q)
  cat("Root prior:\n")
  root <- x$root_prior
  names(root) <- c(0, 1)
  print(root)
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
  likes <- l[t$root_id, ] + root
  mx <- max(likes)
  lk <- log(sum(exp(likes - mx))) + ev_prob + mx
  return(lk)
}
