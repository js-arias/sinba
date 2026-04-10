#' @import nloptr

# Maximum transition rate.
# This maximum is used as very large transition rates
# produce matrices that cannot be exponentiated.
maximum_transition_rate <- 1000

#' @export
#' @title
#' Maximum Likelihood Estimation of the Sinba Model
#'
#' @description
#' `fit_sinba()` searches for the maximum likelihood estimate
#' with the Sinba model for one or two traits.
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
#' @param pi_x The transition probability at the birth of x trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param pi_y The transition probability at the birth of y trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param root Root prior probabilities.
#'   By default,
#'   it uses a FitzJohn prior.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_sinba <- function(
    tree, data, model = NULL,
    pi_x = NULL, pi_y = NULL,
    root = NULL,
    opts = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop("fit_sinba: `model` must be an object of class \"sinba_model\".")
  }
  if (model$traits == 1) {
    return(fit_sinba_single(tree, data, model, pi_x, root, opts))
  }
  mQ <- model$model
  k <- max(mQ) + 2

  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$trait_states[["x"]]$states)
  }
  if (length(pi_y) == 0) {
    pi_y <- default_pi_vector(model$trait_states[["y"]]$states)
  }
  if (length(pi_x) != length(model$trait_states[["x"]]$states)) {
    stop("fit_sinba: invalid pi_x: size different to number of states")
  }
  if (length(pi_y) != length(model$trait_states[["y"]]$states)) {
    stop("fit_sinba: invalid pi_y: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }
  if (sum(pi_y) != 0) {
    pi_y <- pi_y / sum(pi_y)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  # separate birth parameters
  # from transition (traditional) parameters
  transition_start <- 3

  max_rate <- maximum_transition_rate / max(t$ages)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[transition_start:length(p)] > max_rate)) {
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

      Q <- from_model_to_Q(mQ, p[transition_start:length(p)])
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        r, pi_x, pi_y, root
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

  root_states <- c("00", "01", "10", "11")

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (sum(root) != 0) {
      is_valid_root <- FALSE
      for (i in seq_len(length(model$states))) {
        obs <- model$observed[[model$states[i]]]
        if ((obs == root_states[r]) && (root[r] > 0)) {
          is_valid_root <- TRUE
        }
      }
      if (!is_valid_root) {
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
  root_state <- root_states[res$root]

  # retrieve the scenario
  sc <- scenario(res$root, 1)
  if (res$solution[2] < res$solution[1]) {
    # second trait is the oldest one
    sc <- scenario(res$root, 2)
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
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
  cat("Traits: ", x$model$traits, ".\n", sep = "")

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

  if (x$model$traits == 1) {
    b <- x$birth
    ed <- 0
    if (b$node > length(x$tree$tip.label) + 1) {
      ed <- which(x$tree$edge[, 2] == b$node)
    }

    cat(paste("- Edge ", ed,
      " (leads to node ", b$node, ") time ", round(b$age, digits), "\n",
      sep = ""
    ))
  } else {
    cat("Birth events:\n")
    n <- colnames(x$data)
    b1 <- x$births[[1]]
    e1 <- 0
    if (b1$node > length(x$tree$tip.label) + 1) {
      e1 <- which(x$tree$edge[, 2] == b1$node)
    }
    cat(paste("- Trait ", n[2], " Edge ", e1,
      " (leads to node ", b1$node, ") time ", round(b1$age, digits), "\n",
      sep = ""
    ))
    b2 <- x$births[[2]]
    e2 <- 0
    if (b2$node > length(x$tree$tip.label) + 1) {
      e2 <- which(x$tree$edge[, 2] == b2$node)
    }
    cat(paste("- Trait ", n[3], " Edge ", e2,
      " (leads to node ", b2$node, ") time ", round(b2$age, digits), "\n",
      sep = ""
    ))
  }

  if (is.infinite(x$logLik)) {
    return()
  }

  cat("Rates:\n")
  Q <- x$Q
  rownames(Q) <- states
  colnames(Q) <- states
  print(Q)

  if (x$model$traits == 1) {
    return()
  }
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
#' Maximum Likelihood Estimation With a Simultaneous Birth
#'
#' @description
#' `fit_simultaneous()` searches for the maximum likelihood estimate
#' with the Sinba model assuming that the birth of two traits
#' is simultaneous.
#' It estimate the values of the parameters for the Q matrix
#' as well as the birth event.
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
#' @param pi_xy The transition probability at the birth
#'   of x and y traits.
#'   If NULL it will set 1.0 for the state 11.
#' @param root Root prior probabilities.
#'   By default is uses FitzJohn prior.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_simultaneous <- function(
    tree, data, model = NULL,
    pi_xy = NULL,
    root = NULL,
    opts = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_simultaneous: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits != 2) {
    stop("fit_simultaneous: `model` should be a two traits model")
  }
  mQ <- model$model
  k <- max(mQ) + 1

  if (!inherits(tree, "phylo")) {
    stop("fit_simultaneous: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (length(pi_xy) == 0) {
    pi_xy <- rep(0, length(model$states))
    for (i in seq_len(length(model$states))) {
      if (model$observed[[model$states[i]]] == "11") {
        pi_xy[i] <- 1
      }
    }
  }
  if (length(pi_xy) != length(model$states)) {
    stop("fit_simultaneous: invalid pi_xy: size different to number of states")
  }
  if (sum(pi_xy) != 0) {
    pi_xy <- pi_xy / sum(pi_xy)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  # simultaneous birth parameter
  transition_start <- 2

  max_rate <- maximum_transition_rate / max(t$ages)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[transition_start:length(p)] > max_rate)) {
        return(Inf)
      }

      births <- list()
      for (i in 1:2) {
        n <- get_node_by_len(t, p[1], yn[1])
        if (n <= 0) {
          return(Inf)
        }
        b <- list(
          node = n,
          age = t$br_len[n] + p[1] - t$age[n]
        )
        births[[i]] <- b
      }

      Q <- from_model_to_Q(mQ, p[transition_start:length(p)])
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        r, pi_xy, pi_xy, root, TRUE
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

  root_states <- c("00", "01", "10", "11")

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (sum(root) != 0) {
      is_valid_root <- FALSE
      for (i in seq_len(length(model$states))) {
        obs <- model$observed[[model$states[i]]]
        if ((obs == root_states[r]) && (root[r] > 0)) {
          is_valid_root <- TRUE
        }
      }
      if (!is_valid_root) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    # check all possible birth events
    ev <- birth_events(t, youngest)
    for (i in sample(seq_len(length(ev)))) {
      e <- min(ev[[i]])
      fn <- like_func(e, r)
      par <- c(runif(1, max = t$age[e]), runif(k - 1))
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
    n <- get_node_by_len(t, res$solution[1], res$ev_node[1])
    if (n <= 0) {
      return(Inf)
    }
    b <- list(
      node = n,
      age = t$br_len[n] + res$solution[1] - t$age[n]
    )
    births[[i]] <- b
  }
  root_state <- root_states[res$root]

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

#' @export
#' @title
#' Maximum Likelihood Estimation With Mixed Births
#'
#' @description
#' `fit_mixed()` searches for the maximum likelihood estimate
#' with the Sinba model assuming that one of the births
#' is fixed at the root,
#' and the other is free.
#' It estimates the values of the parameters for the Q matrix
#' as well as the location of the birth event
#' of the free trait.
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
#' @param trait The free trait
#'   (i.e., the birth to be estimated).
#'   By default is "x",
#'   for the first trait.
#'   Use "y" for the second trait.
#' @param pi_trait The transition probability for the birth
#'   of the free trait.
#' @param root Root prior probabilities.
#'   By default is uses FitzJohn prior.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_mixed <- function(
    tree, data, model = NULL,
    trait = "x",
    pi_trait = NULL,
    root = NULL,
    opts = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_mixed: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits != 2) {
    stop("fit_mixed: `model` should be a two traits model")
  }
  mQ <- model$model
  k <- max(mQ) + 1

  if (!inherits(tree, "phylo")) {
    stop("fit_mixed: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  pi_x <- default_pi_vector(model$trait_states[["x"]]$states)
  pi_y <- default_pi_vector(model$trait_states[["y"]]$states)

  if (trait == "x") {
    if (length(pi_trait) > 0) {
      pi_x <- pi_trait
    }
  } else if (trait == "y") {
    if (length(pi_trait) > 0) {
      pi_y <- pi_trait
    }
  } else {
    stop("fit_mixed: unknown trait")
  }
  if (length(pi_x) != length(model$trait_states[["x"]]$states)) {
    stop("fit_mixed: invalid pi_trait: size different to number of states")
  }
  if (length(pi_y) != length(model$trait_states[["y"]]$states)) {
    stop("fit_mixed: invalid pi_trait: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }
  if (sum(pi_y) != 0) {
    pi_y <- pi_y / sum(pi_y)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  # separate birth parameters
  # from transition (traditional) parameters
  transition_start <- 2

  max_rate <- maximum_transition_rate / max(t$ages)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[transition_start:length(p)] > max_rate)) {
        return(Inf)
      }

      births <- list()
      n <- get_node_by_len(t, p[1], yn[1])
      if (n <= 0) {
        return(Inf)
      }
      b_trait <- list(
        node = n,
        age = t$br_len[n] + p[1] - t$age[n]
      )
      b_root <- list(
        node = t$root_id,
        age = 0
      )
      if (trait == "x") {
        births[[1]] <- b_trait
        births[[2]] <- b_root
      } else {
        births[[1]] <- b_root
        births[[2]] <- b_trait
      }

      Q <- from_model_to_Q(mQ, p[transition_start:length(p)])
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        r, pi_x, pi_y, root
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

  root_states <- c("00", "01", "10", "11")

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (sum(root) != 0) {
      is_valid_root <- FALSE
      for (i in seq_len(length(model$states))) {
        obs <- model$observed[[model$states[i]]]
        if ((obs == root_states[r]) && (root[r] > 0)) {
          is_valid_root <- TRUE
        }
      }
      if (!is_valid_root) {
        next
      }
    }

    youngest <- youngest_birth_event(t, et, root_states[r])
    # check all possible birth events
    ev <- birth_events(t, youngest)
    for (i in sample(seq_len(length(ev)))) {
      e <- ev[[i]][1]
      if (trait == "y") {
        e <- ev[[i]][2]
      }
      fn <- like_func(e, r)
      par <- c(runif(1, max = t$age[e]), runif(k - 1))
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
  n <- get_node_by_len(t, res$solution[1], res$ev_node[1])
  if (n <= 0) {
    return(Inf)
  }
  b_trait <- list(
    node = n,
    age = t$br_len[n] + res$solution[1] - t$age[n]
  )
  b_root <- list(
    node = t$root_id,
    age = 0
  )
  if (trait == "x") {
    births[[1]] <- b_trait
    births[[2]] <- b_root
  } else {
    births[[1]] <- b_root
    births[[2]] <- b_trait
  }
  root_state <- root_states[res$root]

  # retrieve the scenario
  sc <- scenario(res$root, 1)
  if (trait == "x") {
    # second trait is the oldest one
    sc <- scenario(res$root, 2)
  }

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenarios = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
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
#' @param births A list with one or two element,
#'   each element being a list that define a birth event,
#'   with a fields `node` indicating the birth of the trait
#'   and `age` indicating the time from the start of the edge
#'   in which the event happens.
#'   If the model is for two traits,
#'   and there is single defined birth,
#'   the same birth will be used for both traits.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param pi_x The transition probability at the birth of x trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param pi_y The transition probability at the birth of y trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param root Root prior probabilities.
#'   By default,
#'   it uses a FitzJohn prior.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_births <- function(
    tree, data, births, model = NULL,
    pi_x = NULL, pi_y = NULL,
    root = NULL,
    opts = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_fixed_births: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits == 1) {
    return(fit_sinba_single_fixed_birth(
      tree, data, births, model, pi_x, root, opts
    ))
  }
  mQ <- model$model
  k <- max(mQ)

  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_births: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$trait_states[["x"]]$states)
  }
  if (length(pi_y) == 0) {
    pi_y <- default_pi_vector(model$trait_states[["y"]]$states)
  }
  if (length(pi_x) != length(model$trait_states[["x"]]$states)) {
    stop("fit_fixed_births: invalid pi_x: size different to number of states")
  }
  if (length(pi_y) != length(model$trait_states[["y"]]$states)) {
    stop("fit_fixed_births: invalid pi_y: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }
  if (sum(pi_y) != 0) {
    pi_y <- pi_y / sum(pi_y)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  max_rate <- maximum_transition_rate / max(t$ages)

  if (length(births) < 1) {
    stop("fit_fixed_births: at least ine birth event is required")
  }
  if (length(births) < 2) {
    births[[2]] <- births[[1]]
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
      if (any(p > max_rate)) {
        return(Inf)
      }

      Q <- from_model_to_Q(mQ, p)
      lk <- sinba_like(
        t, Q, model, births, xt, cond,
        r, pi_x, pi_y, root
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

  root_states <- c("00", "01", "10", "11")

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (sum(root) != 0) {
      is_valid_root <- FALSE
      for (i in seq_len(length(model$states))) {
        obs <- model$observed[[model$states[i]]]
        if ((obs == root_states[r]) && (root[r] > 0)) {
          is_valid_root <- TRUE
        }
      }
      if (!is_valid_root) {
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

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
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
#' of birth events
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
#' @param pi_x The transition probability at the birth of x trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param pi_y The transition probability at the birth of y trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param root Root prior probabilities.
#'   By default,
#'   it uses a FitzJohn prior.
#' @param opts User defined parameters for the optimization
#'   with the `nloptr` package.
#'   By default it attempts a reasonable set of options.
fit_fixed_matrix <- function(
    tree, data, rate_mat,
    model = NULL,
    pi_x = NULL, pi_y = NULL,
    root = NULL,
    opts = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fit_fixed_matrix: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits == 1) {
    return(fit_sinba_single_fixed_matrix(
      tree, data, rate_mat, model, pi_x, root, opts
    ))
  }

  k <- 2

  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_matrix: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop("fit_fixed_matrix: `rate_mat` must have the same size as the `model`")
  }

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$trait_states[["x"]]$states)
  }
  if (length(pi_y) == 0) {
    pi_y <- default_pi_vector(model$trait_states[["y"]]$states)
  }
  if (length(pi_x) != length(model$trait_states[["x"]]$states)) {
    stop("fit_fixed_matrix: invalid pi_x: size different to number of states")
  }
  if (length(pi_y) != length(model$trait_states[["y"]]$states)) {
    stop("fit_fixed_matrix: invalid pi_y: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }
  if (sum(pi_y) != 0) {
    pi_y <- pi_y / sum(pi_y)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

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
        r, pi_x, pi_y, root
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

  root_states <- c("00", "01", "10", "11")

  res <- list(objective = Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (sum(root) != 0) {
      is_valid_root <- FALSE
      for (i in seq_len(length(model$states))) {
        obs <- model$observed[[model$states[i]]]
        if ((obs == root_states[r]) && (root[r] > 0)) {
          is_valid_root <- TRUE
        }
      }
      if (!is_valid_root) {
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

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenarios = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

#' @export
#' @title
#' Estimate the Likelihood for a Fixed Transition Matrix
#' and Fixed Births
#'
#' @description
#' `fixed_sinba()` calculates the likelihood
#' for a fixed transition matrix
#' and fixed birth events
#' under the Sinba model.
#'
#' @param tree A phylogenetic tree of class "phylo".
#' @param data A data frame with the data.
#'   The first column should contain the taxon names,
#'   The second and third column contains the data,
#'   coded as 0 and 1.
#'   Any other column will be ignored.
#' @param rate_mat Rate matrix for the traits with the full process.
#' @param births A list with one or two element,
#'   each element being a list that define a birth event,
#'   with a fields `node` indicating the birth of the trait
#'   and `age` indicating the time from the start of the edge
#'   in which the event happens.
#'   If the model is for two traits,
#'   and there is single defined birth,
#'   the same birth will be used for both traits.
#' @param model A model build with `new_model()`,
#'   `new_hidden_model()`,
#'   or `new_rates_model()`.
#'   By default it uses the independent model.
#' @param pi_x The transition probability at the birth of x trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param pi_y The transition probability at the birth of y trait.
#'   If NULL it will set 1.0 for the state 1.
#' @param root Root prior probabilities.
#'   By default,
#'   it uses a FitzJohn prior.
fixed_sinba <- function(
    tree, data, rate_mat, births,
    model = NULL,
    pi_x = NULL, pi_y = NULL,
    root = NULL) {
  if (is.null(model)) {
    model <- new_model("IND")
  }
  if (!inherits(model, "sinba_model")) {
    stop(
      "fixed_sinba: `model` must be an object of class \"sinba_model\"."
    )
  }
  if (model$traits == 1) {
    return(fixed_sinba_single(
      tree, data, rate_mat, births, model, pi_x, root
    ))
  }

  if (!inherits(tree, "phylo")) {
    stop("fixed_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  if (is.null(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fixed_sinba: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop("fixed_sinba: `rate_mat` must have the same size as the `model`")
  }

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$trait_states[["x"]]$states)
  }
  if (length(pi_y) == 0) {
    pi_y <- default_pi_vector(model$trait_states[["y"]]$states)
  }
  if (length(pi_x) != length(model$trait_states[["x"]]$states)) {
    stop("fixed_sinba: invalid pi_x: size different to number of states")
  }
  if (length(pi_y) != length(model$trait_states[["y"]]$states)) {
    stop("fixed_sinba: invalid pi_y: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }
  if (sum(pi_y) != 0) {
    pi_y <- pi_y / sum(pi_y)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 2)
  cond <- set_conditionals(t, et, model)

  if (length(births) < 1) {
    stop("fixed_sinba: at least ine birth event is required")
  }
  if (length(births) < 2) {
    births[[2]] <- births[[1]]
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
        "fixed_sinba: invalid value for field `age` in `births` event %d",
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
        k = 0,
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

  root_states <- c("00", "01", "10", "11")

  res <- list(objective = -Inf)
  for (r in sample(seq_len(length(root_states)))) {
    if (sum(root) != 0) {
      is_valid_root <- FALSE
      for (i in seq_len(length(model$states))) {
        obs <- model$observed[[model$states[i]]]
        if ((obs == root_states[r]) && (root[r] > 0)) {
          is_valid_root <- TRUE
        }
      }
      if (!is_valid_root) {
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
    lk <- sinba_like(
      t, rate_mat, model, births, xt, cond,
      r, pi_x, pi_y, root
    )
    if (lk > res$objective) {
      res <- list(
        objective = lk,
        root = r
      )
    }
  }

  # if no birth sequence is compatible with birth events
  # the likelihood is 0
  if (is.infinite(res$objective)) {
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
    logLik = res$objective,
    k = 0,
    model = model,
    Q = q,
    births = births,
    root = root_state,
    scenarios = sc,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_like calculates the likelihood
# of the sinba model.
sinba_like <- function(
    t, Q, model, births, xt, cond, root, pi_x, pi_y, pi_root,
    simultaneous = FALSE) {
  l <- sinba_cond(t, Q, model, births, xt, cond, root, pi_x, pi_y, simultaneous)

  return(add_root_prior(l[t$root_id, ], pi_root))
}

add_root_prior <- function(likes, pi) {
  scaled <- exp(likes - max(likes))

  # FitzJohn root
  if (sum(pi) == 0) {
    fitz <- scaled / sum(scaled)
    like <- log(sum(fitz * scaled)) + max(likes)
    return(like)
  }

  like <- log(sum(pi * scaled)) + max(likes)
  return(like)
}

# sinba_cond return the conditional likelihoods
# under the sinba model.
sinba_cond <- function(
    t, Q, model, births, xt, cond, root, pi_x, pi_y, simultaneous) {
  Q <- normalize_Q(Q)

  root_vector <- rep(0, 4)
  root_vector[root] <- 1
  m_PI_semi <- matrix()
  m_PI_root <- matrix()
  pi_semi <- pi_x
  pi_root <- pi_y

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
    anc_vector <- active_ancestor_vector(sc)
    m_PI_semi <- build_pi_matrix(model, "y", anc_vector)
    m_PI_root <- build_pi_matrix(model, "x", root_vector)
    pi_semi <- pi_y
    pi_root <- pi_x
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
    anc_vector <- active_ancestor_vector(sc)
    m_PI_semi <- build_pi_matrix(model, "x", anc_vector)
    m_PI_root <- build_pi_matrix(model, "y", root_vector)
    pi_semi <- pi_x
    pi_root <- pi_y
  }

  st <- as.integer(active_status(t, ev$first$node, ev$second$node))

  if (simultaneous) {
    m_PI_root <- matrix(0, nrow = nrow(model$model), ncol = ncol(model$model))
    m_PI_root[root, ] <- seq_len(length(pi_x))
    l <- sinba_simultaneous(
      xt$parent, xt$nodes, st, xt$branch,
      cond,
      ev$first$age,
      ev$second$Q,
      pi_x,
      m_PI_root
    )
    return(l)
  }

  l <- sinba_conditionals(
    xt$parent, xt$nodes, st, xt$branch,
    cond,
    ev$first$age, ev$second$age,
    ev$first$Q, ev$second$Q,
    pi_semi, pi_root,
    m_PI_semi, m_PI_root
  )
  return(l)
}

# provide a default nloptr options
def_nloptr_opts <- function(k) {
  v <- 1e-12
  return(list(
    "algorithm" = "NLOPT_LN_SBPLX",
    xtol_abs = rep(v, k),
    maxeval = 100000
  ))
}
