#' @export
#' @title
#' Maximum Likelihood Estimation of the Pagel's Model Under the Sinba Model
#'
#' @description
#' `fit_pagel()` searches for the maximum likelihood estimate
#' of a given model for two traits
#' using the Pagel (1994) model.
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
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_pagel <- function(
    tree, data, model = NULL,
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
  k <- max(mQ)

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  # closure for the likelihood function
  like_func <- function() {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p > 1000)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p)
      lk <- mk_like(t, Q, xt, cond)
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

  fn <- like_func()
  par <- c(runif(k))
  rr <- nloptr::nloptr(
    x0 = par,
    eval_f = fn,
    opts = opts
  )

  q <- from_model_to_Q(mQ, rr$solution)
  q <- normalize_Q(q)
  obj <- list(
    logLik = -rr$objective,
    k = k,
    model = model,
    Q = q,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_mk"
  return(obj)
}

#' @export
#' @title Extract Log-Likelihood from a "fit_mk" Object
#'
#' @description
#' This method implements the `logLik` method
#' on a "fit_mk" object.
#'
#' @param object An object of type "fit_mk".
#' @param ... Additional arguments are unused.
logLik.fit_mk <- function(object, ...) {
  l <- object$logLik
  attr(l, "df") <- object$k
  attr(l, "nobs") <- 2 * length(object$tree$tip.label)
  class(l) <- "logLik"
  return(l)
}

#' @export
#' @title Basic Print For a "fit_mk" Object
#'
#' @description
#' This method implements the `print` method
#' on a `fit_mk` object.
#'
#' @param x An of type "fit_mk".
#' @param digits The number of digits for decimal output.
#' @param ... Additional arguments are unused.
print.fit_mk <- function(x, digits = 6, ...) {
  cat("Mk: Fit\n")

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

  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- states
  colnames(Q) <- states
  print(Q)
}

mk_like <- function(t, Q, xt, cond) {
  Q <- normalize_Q(Q)
  l <- full_conditionals(
    xt$parent, xt$nodes, xt$branch, cond, Q
  )
  scaled <- exp(l[t$root_id, ] - max(l[t$root_id, ]))
  fitz <- scaled / sum(scaled)
  like <- log(sum(fitz * scaled)) + max(l[t$root_id, ])
  return(like)
}
