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
#' @param model Model of evolution for the traits.
#'   By default it uses the independent model ("IND").
#'   The standard model for correlated traits
#'   sn the "DEP" model.
#'   In the "xy" model it is assumed that traits are correlated
#'   and the process for both traits
#'   start simultaneously,
#'   so there is no semi-active process,
#'   and only a single birth.
#'   In the "ER" model both traits have equal rates
#'   in any direction;
#'   the "ER2" model also has equal rates,
#'   but rates are different for each trait;
#'   If the "SYM" model changes between states are equal.
#'   There a two full dependant models,
#'   "x" for a model in which trait x depends on y;
#'   and "y" in which trait y depends on x.
#'   The "coll" model collapse (i.e., removes)
#'   entries for unobserved traits.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba <- function(tree, data, model = "IND", root = NULL, opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  cond <- init_conditionals(t, data, 2)

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_sinba: invalid size for `root` vector.")
  }
  root <- root / sum(root)

  if (model == "x" || model == "y") {
    obj <- sinba_dep(t, cond, model, root, opts)
    obj$data <- data
    obj$tree <- tree
    return(obj)
  }
  if (model == "ARD" || model == "xy") {
    obj <- sinba_corr(t, cond, model, root, opts)
    obj$data <- data
    obj$tree <- tree
    return(obj)
  }

  mQ <- model_matrix(model)
  if (model == "coll") {
    mQ <- collapse_model(observed(t, data))
  }
  k <- max(mQ) + 2

  youngest <- youngest_birth_event(t, cond)
  ev_prob <- 2 * log(prob_birth(t))

  # closure for the likelihood function
  like_func <- function(yn) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[3:length(p)] > 1000)) {
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

      # update root priors
      root_prob <- update_root(t, root, births, youngest)
      if (all(root_prob == 0)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p[3:length(p)])
      semi_Q <- Q
      lk <- sinba_like(
        t, Q, semi_Q, births, xt, cond,
        log(root_prob), ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      # set the upper bound using the number of replicates
      xtol_abs = rep(v, 3),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  # check all possible birth events
  ev <- birth_events(t, youngest)

  res <- list(objective = Inf)
  for (i in sample(seq_len(length(ev)))) {
    e <- ev[[i]]
    fn <- like_func(e)
    par <- c(
      runif(1, max = t$age[e[1]]),
      runif(1, max = t$age[e[2]]), runif(k - 2)
    )
    r <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    if (r$objective < res$objective) {
      res <- r
      res$ev_nodes <- e
    }
  }

  q <- from_model_to_Q(mQ, res$solution[3:length(res$solution)])
  q <- normalize_Q(q)
  births <- list()
  for (i in 1:2) {
    n <- get_node_by_len(t, res$solution[i], res$ev_node[i])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[i] - t$age[n]
    )
    births[[i]] <- b
  }
  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
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
  cat(paste("Model: ", x$model, ".\n", sep = ""))
  cat(paste("Free parameters = ", x$k, ".\n", sep = ""))

  aic <- 2 * x$k - 2 * x$logLik
  aicc <- aic + (2 * x$k * x$k + 2 * x$k) /
    (2 * length(x$tree$tip.label) - x$k - 1)
  fit <- c(x$logLik, aic, aicc)
  names(fit) <- c("logLik", "AIC", "AICc")
  print(fit)

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
  rownames(Q) <- x$states
  colnames(Q) <- x$states
  print(Q)
  age1 <- phylo_node_age(x$tree, b1$node)
  age2 <- phylo_node_age(x$tree, b2$node)
  if (b1$node == b2$node) {
    age1 <- b1$age
    age2 <- b2$age
  }
  if (age1 != age2) {
    cat("Semi-active process:\n")
    sQ <- x$Q
    if (!is.null(x$semi_Q)) {
      sQ <- x$semi_Q
    }
    semi <- matrix(0, nrow = 2, ncol = 2)
    if (age1 < age2) {
      # first trait born first
      semi[1, 2] <- sQ[1, 3]
      if (sQ[1, 3] == 0) {
        semi[1, 2] <- sQ[2, 4]
      }
      semi[2, 1] <- sQ[3, 1]
      if (sQ[3, 1] == 0) {
        semi[2, 1] <- sQ[4, 2]
      }
      rownames(semi) <- c("0*", "1*")
      colnames(semi) <- c("0*", "1*")
    } else {
      # second trait born first
      semi[1, 2] <- sQ[1, 2]
      if (sQ[1, 2] == 0) {
        semi[1, 2] <- sQ[3, 4]
      }
      semi[2, 1] <- sQ[2, 1]
      if (sQ[2, 1] == 0) {
        semi[2, 1] <- sQ[4, 3]
      }
      rownames(semi) <- c("*0", "*1")
      colnames(semi) <- c("*0", "*1")
    }
    semi <- normalize_Q(semi)
    print(semi)
  }
  cat("Root prior:\n")
  root <- x$root_prior
  names(root) <- x$states
  print(root)
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
#' @param model Model of evolution for the traits.
#'   By default it uses the independent model ("IND").
#'   Other valid values are "ARD"
#'   or "xy" for a fully correlated model;
#'   "ER" for a model in which both traits have equal rates
#'   in any direction;
#'   "ER2" for an equal rates model,
#'   but rates are different for each trait;
#'   "SYM" for the symmetric model
#'   in which changes between states are equal;
#'   "x" for a model in which trait x depends on y;
#'   and "y" in which trait y depends on x.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_births <- function(
    tree, data, births, model = "IND", root = NULL,
    opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_births: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  cond <- init_conditionals(t, data, 2)

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_fixed_births: invalid size for `root` vector.")
  }
  root <- root / sum(root)

  mQ <- model_matrix(model)
  if (model == "coll") {
    mQ <- collapse_model(observed(t, data))
  }
  k <- max(mQ)

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
        states = c("00", "01", "10", "11"),
        root_prior = root,
        data = data,
        tree = tree
      )
      class(obj) <- "fit_sinba"
      return(obj)
    }
  }

  youngest <- youngest_birth_event(t, cond)

  # check for valid events
  root_prob <- update_root(t, root, births, youngest)

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (all(root_prob == 0)) {
    obj <- list(
      logLik = -Inf,
      k = k,
      model = model,
      Q = matrix(nrow = 4, ncol = 4),
      births = births,
      states = c("00", "01", "10", "11"),
      root_prior = root,
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba"
    return(obj)
  }

  # if birth sequence is not compatible with the model
  # the likelihood is 0
  if (model == "x") {
    # the first trait is the dependant trait
    # so it must born after second trait
    if (b1$node == b2$node) {
      if (b1$age < b2$age) {
        obj <- list(
          logLik = -Inf,
          k = k,
          model = model,
          Q = matrix(nrow = 4, ncol = 4),
          births = births,
          states = c("00", "01", "10", "11"),
          root_prior = root,
          data = data,
          tree = tree
        )
        class(obj) <- "fit_sinba"
        return(obj)
      }
    } else {
      if (t$age[b1$node] < t$age[b2$node]) {
        obj <- list(
          logLik = -Inf,
          k = k,
          model = model,
          Q = matrix(nrow = 4, ncol = 4),
          births = births,
          states = c("00", "01", "10", "11"),
          root_prior = root,
          data = data,
          tree = tree
        )
        class(obj) <- "fit_sinba"
        return(obj)
      }
    }
  }
  if (model == "y") {
    # the second trait is the dependant trait
    # so it must born after first trait
    if (b1$node == b2$node) {
      if (b2$age < b1$age) {
        obj <- list(
          logLik = -Inf,
          k = k,
          model = model,
          Q = matrix(nrow = 4, ncol = 4),
          births = births,
          states = c("00", "01", "10", "11"),
          root_prior = root,
          data = data,
          tree = tree
        )
        class(obj) <- "fit_sinba"
        return(obj)
      }
    } else {
      if (t$age[b2$node] < t$age[b1$node]) {
        obj <- list(
          logLik = -Inf,
          k = k,
          model = model,
          Q = matrix(nrow = 4, ncol = 4),
          births = births,
          states = c("00", "01", "10", "11"),
          root_prior = root,
          data = data,
          tree = tree
        )
        class(obj) <- "fit_sinba"
        return(obj)
      }
    }
  }

  ev_prob <- 2 * log(prob_birth(t))

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
      semi_Q <- Q
      lk <- sinba_like(
        t, Q, semi_Q, births, xt, cond,
        log(root_prob), ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      # set the upper bound using the number of replicates
      xtol_abs = rep(v, 3),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  fn <- like_func()
  par <- c(runif(k))
  res <- nloptr::nloptr(
    x0 = par,
    eval_f = fn,
    opts = opts
  )

  q <- from_model_to_Q(mQ, res$solution)
  q <- normalize_Q(q)
  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
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
#' @param semi_mat Rate matrix for the traits with the semi-active process.
#'   If it is NULL,
#'   it will use the same parameter values
#'   from the full process.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_matrix <- function(
    tree, data, rate_mat, semi_mat = NULL,
    root = NULL, opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_matrix: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  cond <- init_conditionals(t, data, 2)

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_fixed_matrix: invalid size for `root` vector.")
  }
  root <- root / sum(root)

  if (is.null(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != 4) {
    stop("fit_fixed_matrix: `rate_mat` should be 4x4")
  }

  if (is.null(semi_mat)) {
    semi_mat <- rate_mat
  }
  if (nrow(semi_mat) != ncol(semi_mat)) {
    stop("fit_fixed_matrix: `semi_mat` must be a square matrix")
  }
  if (nrow(semi_mat) != 4) {
    stop("fit_fixed_matrix: `semi_mat` should be 4x4")
  }

  k <- 2

  youngest <- youngest_birth_event(t, cond)
  ev_prob <- 2 * log(prob_birth(t))

  # closure for the likelihood function
  like_func <- function(yn) {
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

      # update root priors
      root_prob <- update_root(t, root, births, youngest)
      if (all(root_prob == 0)) {
        return(Inf)
      }

      lk <- sinba_like(
        t, rate_mat, semi_mat, births, xt, cond,
        log(root_prob), ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      # set the upper bound using the number of replicates
      xtol_abs = rep(v, 3),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  # check all possible birth events
  ev <- birth_events(t, youngest)

  res <- list(objective = Inf)
  for (i in sample(seq_len(length(ev)))) {
    e <- ev[[i]]
    fn <- like_func(e)
    par <- c(
      runif(1, max = t$age[e[1]]),
      runif(1, max = t$age[e[2]]), runif(k - 2)
    )
    r <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    if (r$objective < res$objective) {
      res <- r
      res$ev_nodes <- e
    }
  }

  q <- normalize_Q(rate_mat)
  sq <- normalize_Q(semi_mat)
  births <- list()
  for (i in 1:2) {
    n <- get_node_by_len(t, res$solution[i], res$ev_node[i])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[i] - t$age[n]
    )
    births[[i]] <- b
  }
  obj <- list(
    logLik = -res$objective,
    k = k,
    model = "fixed",
    Q = q,
    semi_Q = sq,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
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
#' @param semi_mat Rate matrix for the traits with the semi-active process.
#'   If it is NULL,
#'   it will use the same parameter values
#'   from the full process.
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
fixed_sinba <- function(
    tree, data, rate_mat, births,
    semi_mat = NULL, root = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fixed_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  cond <- init_conditionals(t, data, 2)

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fixed_sinba: invalid size for `root` vector.")
  }
  root <- root / sum(root)

  if (is.null(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != 4) {
    stop("fixed_sinba: `rate_mat` should be 4x4")
  }

  if (is.null(semi_mat)) {
    semi_mat <- rate_mat
  }
  if (nrow(semi_mat) != ncol(semi_mat)) {
    stop("fixed_sinba: `semi_mat` must be a square matrix")
  }
  if (nrow(semi_mat) != 4) {
    stop("fixed_sinba: `semi_mat` should be 4x4")
  }

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
        Q = normalize_Q(rate_mat),
        semi_Q = semi_mat,
        births = births,
        states = c("00", "01", "10", "11"),
        root_prior = root,
        data = data,
        tree = tree
      )
      class(obj) <- "fixed_sinba"
      return(obj)
    }
  }

  youngest <- youngest_birth_event(t, cond)

  # check for valid events
  root_prob <- update_root(t, root, births, youngest)

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (all(root_prob == 0)) {
    obj <- list(
      logLik = -Inf,
      Q = normalize_Q(rate_mat),
      semi_Q = semi_mat,
      births = births,
      states = c("00", "01", "10", "11"),
      root_prior = root,
      data = data,
      tree = tree
    )
    class(obj) <- "fixed_sinba"
    return(obj)
  }

  ev_prob <- 2 * log(prob_birth(t))

  xt <- tree_to_cpp(t)
  lk <- sinba_like(
    t, rate_mat, semi_mat, births, xt, cond,
    log(root_prob), ev_prob
  )

  obj <- list(
    logLik = lk,
    Q = normalize_Q(rate_mat),
    semi_Q = semi_mat,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
    data = data,
    tree = tree
  )
  class(obj) <- "fixed_sinba"
  return(obj)
}

#' @export
#' @title Basic Print For a "fixed_sinba" Object
#'
#' @description
#' This method implements the `print` method
#' on a `fixed_sinba` object.
#'
#' @param x An object of type "fixed_sinba".
#' @param digits The number of digits for decimal output.
#' @param ... Additional arguments are unused.
print.fixed_sinba <- function(x, digits = 6, ...) {
  cat("Sinba: Fixed Rate Matrix\n")
  cat(paste("Log-Likelihood = ", round(x$logLik, digits), "\n", sep = ""))
  cat("Birth events:\n")
  n <- colnames(x$data)
  b1 <- x$births[[1]]
  cat(paste("- Trait ", n[2], " Node ",
    b1$node, " time ", round(b1$age, digits), "\n",
    sep = ""
  ))
  b3 <- x$births[[2]]
  cat(paste("- Trait ", n[3], " Node ",
    b3$node, " time ", round(b3$age, digits), "\n",
    sep = ""
  ))
  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- x$states
  colnames(Q) <- x$states
  print(Q)
  cat("Semi active rate matrix:\n")
  semi <- x$semi_Q
  rownames(semi) <- x$states
  colnames(semi) <- x$states
  print(semi)
  cat("Root prior:\n")
  root <- x$root_prior
  names(root) <- x$states
  print(root)
}

# sinba_corr makes the maximum likelihood estimation
# of the correlated model.
sinba_corr <- function(t, cond, model, root, opts) {
  mQ <- model_matrix(model)

  youngest <- youngest_birth_event(t, cond)
  ev_prob <- 2 * log(prob_birth(t))

  if (model == "xy") {
    k <- max(mQ) + 1
    # closure for the likelihood function with the fully
    # correlated model.
    # Births happen in the same place.
    xy_like_func <- function(yn) {
      xt <- tree_to_cpp(t)

      return(function(p) {
        if (any(p < 0)) {
          return(Inf)
        }
        if (any(p[2:length(p)] > 1000)) {
          return(Inf)
        }

        n <- get_node_by_len(t, p[1], yn)
        if (n <= 0) {
          return(Inf)
        }
        b <- list(
          node = n,
          age = t$br_len[n] + p[1] - t$age[n]
        )
        births <- list()
        births[[1]] <- b
        births[[2]] <- b

        # update root priors
        root_prob <- update_root(t, root, births, youngest)
        if (all(root_prob == 0)) {
          return(Inf)
        }

        Q <- from_model_to_Q(mQ, p[2:length(p)])
        semi_Q <- Q
        lk <- sinba_like(
          t, Q, semi_Q, births, xt, cond,
          log(root_prob), ev_prob
        )
        return(-lk)
      })
    }

    if (is.null(opts)) {
      v <- 1e-06
      opts <- list(
        "algorithm" = "NLOPT_LN_SBPLX",
        # set the upper bound using the number of replicates
        xtol_abs = rep(v, 3),
        maxeval = 10000
      )
    }
    if (is.null(opts$algorithm)) {
      opts$algorithm <- "NLOPT_LN_SBPLX"
    }

    # check all possible birth events
    ev <- birth_events(t, youngest)

    res <- list(objective = Inf)
    for (i in sample(seq_len(length(ev)))) {
      e <- ev[[i]]
      yn <- e[1]
      if (t$age[yn] > t$age[e[2]]) {
        yn <- e[2]
      }

      fn <- xy_like_func(yn)
      par <- c(
        runif(1, max = t$age[e[1]]), runif(k - 1)
      )
      r <- nloptr::nloptr(
        x0 = par,
        eval_f = fn,
        opts = opts
      )
      if (r$objective < res$objective) {
        res <- r
        res$ev_node <- yn
      }
    }

    q <- from_model_to_Q(mQ, res$solution[2:length(res$solution)])
    q <- normalize_Q(q)
    n <- get_node_by_len(t, res$solution[1], res$ev_node)
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[i] - t$age[n]
    )
    births <- list()
    births[[1]] <- b
    births[[2]] <- b
    obj <- list(
      logLik = -res$objective,
      k = k,
      model = model,
      Q = q,
      births = births,
      states = c("00", "01", "10", "11"),
      root_prior = root
    )
    class(obj) <- "fit_sinba"
    return(obj)
  }

  k <- max(mQ) + 2 + 2 # two trait states and two births

  # closure for the likelihood model
  # with an independent semi-active part.
  ard_like_func <- function(yn, m) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[3:length(p)] > 1000)) {
        return(Inf)
      }

      semi <- matrix()

      # check that births are in agreement with the model
      if (m == "x") {
        # the first trait is the dependent trait
        # so it must born after the second trait
        if (p[1] < p[2]) {
          return(Inf)
        }
        semi <- matrix(c(
          0, 1, 0, 0,
          2, 0, 0, 0,
          0, 0, 0, 1,
          0, 0, 2, 0
        ), nrow = 4, byrow = TRUE)
      }
      if (m == "y") {
        # the second trait is the dependent trait
        # so it must born after the first trait
        if (p[2] < p[1]) {
          return(Inf)
        }
        semi <- matrix(c(
          0, 0, 1, 0,
          0, 0, 0, 1,
          2, 0, 0, 0,
          0, 2, 0, 0
        ), nrow = 4, byrow = TRUE)
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

      # update root priors
      root_prob <- update_root(t, root, births, youngest)
      if (all(root_prob == 0)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p[3:(2 + max(mQ))])
      semi_Q <- from_model_to_Q(semi, p[(3 + max(mQ)):length(p)])
      lk <- sinba_like(
        t, Q, semi_Q, births, xt, cond,
        log(root_prob), ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      # set the upper bound using the number of replicates
      xtol_abs = rep(v, 3),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  # check all possible birth events
  ev <- birth_events(t, youngest)

  res <- list(objective = Inf)
  # check each birth sequence
  for (i in sample(c("x", "y"))) {
    for (j in sample(seq_len(length(ev)))) {
      e <- ev[[j]]
      # check if events are consistent with the model constraint
      if (i == "x") {
        # the first trait is the dependent trait
        # so it must born after the second trait
        if (t$age[e[1]] < t$age[e[2]]) {
          e[2] <- e[1]
        }
      }
      if (i == "y") {
        # the second trait is the dependent trait
        # so it must born after the first trait
        if (t$age[e[2]] < t$age[e[1]]) {
          e[1] <- e[2]
        }
      }

      fn <- ard_like_func(e, i)
      par <- c(
        runif(1, max = t$age[e[1]]),
        runif(1, max = t$age[e[2]]), runif(k - 2)
      )
      r <- nloptr::nloptr(
        x0 = par,
        eval_f = fn,
        opts = opts
      )
      if (r$objective < res$objective) {
        res <- r
        res$ev_nodes <- e
        res$model <- i
      }
    }
  }
  q <- from_model_to_Q(mQ, res$solution[3:(2 + max(mQ))])
  q <- normalize_Q(q)
  semi_q <- matrix()
  if (res$model == "x") {
    semi <- matrix(c(
      0, 1, 0, 0,
      2, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 2, 0
    ), nrow = 4, byrow = TRUE)
    semi_q <- from_model_to_Q(
      semi,
      res$solution[(3 + max(mQ)):length(res$solution)]
    )
  }
  if (res$model == "y") {
    semi <- matrix(c(
      0, 0, 1, 0,
      0, 0, 0, 1,
      2, 0, 0, 0,
      0, 2, 0, 0
    ), nrow = 4, byrow = TRUE)
    semi_q <- from_model_to_Q(
      semi,
      res$solution[(3 + max(mQ)):length(res$solution)]
    )
  }
  semi_q <- normalize_Q(semi_q)

  births <- list()
  for (i in 1:2) {
    n <- get_node_by_len(t, res$solution[i], res$ev_nodes[i])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[i] - t$age[n]
    )
    births[[i]] <- b
  }
  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    semi_Q = semi_q,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_dep makes the maximum likelihood estimation
# of dependant models.
sinba_dep <- function(t, cond, model, root, opts) {
  mQ <- model_matrix(model)
  k <- max(mQ) + 2

  youngest <- youngest_birth_event(t, cond)
  ev_prob <- 2 * log(prob_birth(t))

  # closure for the likelihood function
  like_func <- function(yn, m) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[3:length(p)] > 1000)) {
        return(Inf)
      }

      # check that births are in agreement with the model
      if (m == "x") {
        # the first trait is the dependent trait
        # so it must born after the second trait
        if (p[1] < p[2]) {
          return(Inf)
        }
      }
      if (m == "y") {
        # the second trait is the dependent trait
        # so it must born after the first trait
        if (p[2] < p[1]) {
          return(Inf)
        }
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

      # update root priors
      root_prob <- update_root(t, root, births, youngest)
      if (all(root_prob == 0)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p[3:length(p)])
      semi_Q <- Q
      lk <- sinba_like(
        t, Q, semi_Q, births, xt, cond,
        log(root_prob), ev_prob
      )
      return(-lk)
    })
  }

  if (is.null(opts)) {
    v <- 1e-06
    opts <- list(
      "algorithm" = "NLOPT_LN_SBPLX",
      # set the upper bound using the number of replicates
      xtol_abs = rep(v, 3),
      maxeval = 10000
    )
  }
  if (is.null(opts$algorithm)) {
    opts$algorithm <- "NLOPT_LN_SBPLX"
  }

  # check all possible birth events
  ev <- birth_events(t, youngest)

  res <- list(objective = Inf)
  for (i in sample(seq_len(length(ev)))) {
    e <- ev[[i]]
    # check if events are consistent with the model constraint
    if (model == "x") {
      # the first trait is the dependent trait
      # so it must born after the second trait
      if (t$age[e[1]] < t$age[e[2]]) {
        e[2] <- e[1]
      }
    }
    if (model == "y") {
      # the second trait is the dependent trait
      # so it must born after the first trait
      if (t$age[e[2]] < t$age[e[1]]) {
        e[1] <- e[2]
      }
    }

    fn <- like_func(e, model)
    par <- c(
      runif(1, max = t$age[e[1]]),
      runif(1, max = t$age[e[2]]), runif(k - 2)
    )
    r <- nloptr::nloptr(
      x0 = par,
      eval_f = fn,
      opts = opts
    )
    if (r$objective < res$objective) {
      res <- r
      res$ev_nodes <- e
    }
  }

  q <- from_model_to_Q(mQ, res$solution[3:length(res$solution)])
  q <- normalize_Q(q)
  births <- list()
  for (i in 1:2) {
    n <- get_node_by_len(t, res$solution[i], res$ev_nodes[i])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[i] - t$age[n]
    )
    births[[i]] <- b
  }
  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_like calculates the likelihood
# of the sinba model.
sinba_like <- function(t, Q, semi_Q, births, xt, cond, root_prior, ev_prob) {
  # make sure that Q matrix is valid
  Q[1, 4] <- 0
  Q[2, 3] <- 0
  Q[3, 2] <- 0
  Q[4, 1] <- 0
  Q <- normalize_Q(Q)

  root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))

  likes <- c()
  for (r in seq_len(length(root_prior))) {
    if (is.infinite(root_prior[r])) {
      next
    }

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
      sc <- scenario(r, 1)
      semi <- semi_active_Q(sc, semi_Q)
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
      sc <- scenario(r, 2)
      semi <- semi_active_Q(sc, semi_Q)
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
    likes <- c(likes, l[t$root_id, r] + root_prior[r])
  }
  mx <- max(likes)
  lk <- log(sum(exp(likes - mx))) + ev_prob + mx
  return(lk)
}
