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
#' @param model A model build with `new_model()`.
#'   By default it uses the independent model
#' @param root Root prior probabilities.
#'   By default,
#'   all states will have the same probability.
#' @param root_method Method for root calculation at the root.
#'   By default it use the root prior.
#'   If set as "FitzJohn" it will use the FitzJohn et al. (2009)
#'   method,
#'   in which ancestral states are weighted by its own likelihood.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba <- function(
    tree, data, model = NULL,
    root = NULL, root_method = "prior", opts = NULL) {
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

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_sinba: invalid size for `root` vector.")
  }
  root <- root / sum(root)
  if (root_method == "FitzJohn") {
    root <- rep(1, ncol(cond))
  }

  youngest <- youngest_birth_node(t, et, 2)
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
      root_prob <- set_root_prior(t, model, root, births, youngest)
      if (all(root_prob == 0)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p[3:length(p)])
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        log(root_prob), root_method, ev_prob
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

  # retrieve the scenarios
  root_states <- update_root(t, rep(1, 4), births, youngest)
  root_names <- c("00", "01", "10", "11")
  sc <- c()
  for (i in seq_len(length(root_states))) {
    if (root_states[i] == 0) {
      next
    }
    has_prior <- FALSE
    for (j in seq_len(length(root))) {
      obs <- model$observed[[model$states[j]]]
      if (obs != root_names[i]) {
        next
      }
      if (root[j] != 0) {
        has_prior <- TRUE
      }
    }
    if (!has_prior) {
      next
    }
    v <- scenario(root_states[i], 1)
    if (res$solution[2] < res$solution[1]) {
      # second trait is the oldest one
      v <- scenario(root_states[i], 2)
    }
    sc <- c(sc, v)
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root_prior = root,
    root_method = root_method,
    scenarios = sc,
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
  if (age1 != age2) {
    cat("Semi-active process:\n")
    if (length(x$scenarios) == 0) {
      sQ <- x$semi_Q
      rownames(sQ) <- states
      colnames(sQ) <- states
      sQ <- reduce_matrix(sQ)
      print(normalize_Q(sQ))
    } else {
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

  if (x$root_method == "FitzJohn") {
    cat("Root method: FitzJohn\n")
  } else {
    cat("Root prior:\n")
    root <- x$root_prior
    names(root) <- states
    print(root)
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
#' @param root_method Method for root calculation at the root.
#'   By default it use the root prior.
#'   If set as "FitzJohn" it will use the FitzJohn et al. (2009)
#'   method,
#'   in which ancestral states are weighted by its own likelihood.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_births <- function(
    tree, data, births, model = NULL,
    root = NULL, root_method = "FitzJohn",
    opts = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_births: `tree` must be an object of class \"phylo\".")
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

  if (is.null(root)) {
    root <- rep(1, ncol(cond))
  }
  if (length(root) != ncol(cond)) {
    stop("fit_fixed_births: invalid size for `root` vector.")
  }
  root <- root / sum(root)
  if (root_method == "FitzJohn") {
    root <- rep(1, ncol(cond))
  }

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
        root_prior = root,
        root_method = root_method,
        data = data,
        tree = tree
      )
      class(obj) <- "fit_sinba"
      return(obj)
    }
  }

  youngest <- youngest_birth_node(t, et, 2)

  # check for valid events
  root_prob <- set_root_prior(t, model, root, births, youngest)

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (all(root_prob == 0)) {
    obj <- list(
      logLik = -Inf,
      k = k,
      model = model,
      Q = matrix(nrow = 4, ncol = 4),
      births = births,
      root_prior = root,
      root_method = root_method,
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba"
    return(obj)
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
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        log(root_prob), root_method, ev_prob
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

  fn <- like_func()
  par <- c(runif(k))
  res <- nloptr::nloptr(
    x0 = par,
    eval_f = fn,
    opts = opts
  )

  q <- from_model_to_Q(mQ, res$solution)
  q <- normalize_Q(q)
  # retrieve the scenarios
  root_states <- update_root(t, rep(1, 4), births, youngest)
  root_names <- c("00", "01", "10", "11")
  sc <- c()
  for (i in seq_len(length(root_states))) {
    if (root_states[i] == 0) {
      next
    }
    has_prior <- FALSE
    for (j in seq_len(length(root))) {
      obs <- model$observed[[model$states[j]]]
      if (obs != root_names[i]) {
        next
      }
      if (root[j] != 0) {
        has_prior <- TRUE
      }
    }
    if (!has_prior) {
      next
    }
    v <- scenario(root_states[i], 1)
    if (res$solution[2] < res$solution[1]) {
      # second trait is the oldest one
      v <- scenario(root_states[i], 2)
    }
    sc <- c(sc, v)
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
    root_method = root_method,
    scenarios = sc,
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
#' @param root_method Method for root calculation at the root.
#'   By default it use the root prior.
#'   If set as "FitzJohn" it will use the FitzJohn et al. (2009)
#'   method,
#'   in which ancestral states are weighted by its own likelihood.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_matrix <- function(
    tree, data, rate_mat, semi_mat = NULL,
    root = NULL, root_method = "prior", opts = NULL) {
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
  if (root_method == "FitzJohn") {
    root <- rep(1, ncol(cond))
  }

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
        log(root_prob), root_method, ev_prob
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
    root_method = root_method,
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
#' @param root_method Method for root calculation at the root.
#'   By default it use the root prior.
#'   If set as "FitzJohn" it will use the FitzJohn et al. (2009)
#'   method,
#'   in which ancestral states are weighted by its own likelihood.
fixed_sinba <- function(
    tree, data, rate_mat, births,
    semi_mat = NULL, root = NULL, root_method = "prior") {
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
  if (root_method == "FitzJohn") {
    root <- rep(1, ncol(cond))
  }

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
        root_method = root_method,
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
      root_method = root_method,
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
    log(root_prob), root_method, ev_prob
  )

  obj <- list(
    logLik = lk,
    Q = normalize_Q(rate_mat),
    semi_Q = semi_mat,
    births = births,
    states = c("00", "01", "10", "11"),
    root_prior = root,
    root_method = root_method,
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
  if (x$root_method == "FitzJohn") {
    cat("Root method: FitzJohn\n")
  } else {
    cat("Root prior:\n")
    root <- x$root_prior
    names(root) <- x$states
    print(root)
  }
}

# sinba_xy makes the maximum likelihood estimation
# of the correlated model.
sinba_xy <- function(t, cond, model, root, root_method, opts) {
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
          log(root_prob), root_method, ev_prob
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
      root_prior = root,
      root_method = root_method
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
        log(root_prob), root_method, ev_prob
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
    root_prior = root,
    root_method = root_method
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_dep makes the maximum likelihood estimation
# of dependant models.
sinba_dep <- function(t, cond, model, root, root_method, opts) {
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
        log(root_prob), root_method, ev_prob
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
    root_prior = root,
    root_method = root_method
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_like calculates the likelihood
# of the sinba model.
sinba_like <- function(
    t, Q, model, births, xt, cond,
    root_prior, root_method, ev_prob) {
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
      sc <- scenario(r, 2)
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
    likes <- c(likes, l[t$root_id, r] + root_prior[r])
  }
  mx <- max(likes)
  if (root_method == "FitzJohn") {
    l <- exp(likes - mx)
    d <- sum(l)
    lk <- 0
    for (i in seq_len(length(l))) {
      lk <- lk + l[i] * l[i] / d
    }
    return(log(lk) + ev_prob + mx)
  }
  lk <- log(sum(exp(likes - mx))) + ev_prob + mx
  return(lk)
}
