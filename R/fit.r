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
fit_sinba <- function(tree, data, model = "IND", root = NULL, opts = NULL) {
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

  mQ <- model_matrix(model)
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
      b1 <- births[[1]]
      b2 <- births[[2]]

      root_prob <- root
      y1 <- youngest[[1]]
      y2 <- youngest[[2]]
      if (root[1] > 0) {
        # ancestral state is 00
        # check is trait 1-1 can be born
        if (y1[2] != b1$node) {
          if (!is_parent(t, b1$node, y1[2])) {
            # invalid root state
            root_prob[1] <- 0
          }
        }
        # check if trait 2-1 can be born
        if (y2[2] != b2$node) {
          if (!is_parent(t, b2$node, y2[2])) {
            # invalid root state
            root_prob[1] <- 0
          }
        }
      }
      if (root[2] > 0) {
        # ancestral state is 01
        # check is trait 1-1 can be born
        if (y1[2] != b1$node) {
          if (!is_parent(t, b1$node, y1[2])) {
            # invalid root state
            root_prob[2] <- 0
          }
        }
        # check if trait 2-0 can be born
        if (y2[1] != b2$node) {
          if (!is_parent(t, b2$node, y2[1])) {
            # invalid root state
            root_prob[2] <- 0
          }
        }
      }
      if (root[3] > 0) {
        # ancestral state is 10
        # check is trait 1-0 can be born
        if (y1[1] != b1$node) {
          if (!is_parent(t, b1$node, y1[1])) {
            # invalid root state
            root_prob[3] <- 0
          }
        }
        # check if trait 2-1 can be born
        if (y2[2] != b2$node) {
          if (!is_parent(t, b2$node, y2[2])) {
            # invalid root state
            root_prob[3] <- 0
          }
        }
      }
      if (root[4] > 0) {
        # ancestral state is 11
        # check is trait 1-0 can be born
        if (y1[1] != b1$node) {
          if (!is_parent(t, b1$node, y1[1])) {
            # invalid root state
            root_prob[4] <- 0
          }
        }
        # check if trait 2-0 can be born
        if (y2[1] != b2$node) {
          if (!is_parent(t, b2$node, y2[1])) {
            # invalid root state
            root_prob[4] <- 0
          }
        }
      }
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
  y1 <- youngest[[1]]
  y2 <- youngest[[2]]
  ev <- list()
  for (yn1 in y1) {
    if (yn1 == t$root_id) {
      next
    }
    for (yn2 in y2) {
      if (yn2 == t$root_id) {
        next
      }
      if (yn1 != yn2) {
        if (!is_parent(t, yn1, yn2) && !is_parent(t, yn2, yn1)) {
          next
        }
      }
      ev[[length(ev) + 1]] <- c(yn1, yn2)
    }
  }
  if (length(ev) == 0) {
    # at least one of the traits start at the root
    if (any(y1 != t$root_id)) {
      for (yn1 in y1) {
        if (yn1 == t$root_id) {
          next
        }
        ev[[length(ev) + 1]] <- c(yn1, t$root_id)
      }
    } else if (any(y2 != t$root_id)) {
      for (yn2 in y2) {
        if (yn2 == t$root_id) {
          next
        }
        ev[[length(ev) + 1]] <- c(t$root_id, yn2)
      }
    }
    if (length(ev) == 0) {
      # both traits start at the root
      ev[[length(ev) + 1]] <- c(t$root_id, t$root_id)
    }
  }

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
  cat("Root prior:\n")
  root <- x$root_prior
  names(root) <- x$states
  print(root)
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

  # check for valid events
  root_prob <- root
  youngest <- youngest_birth_event(t, cond)
  y1 <- youngest[[1]]
  y2 <- youngest[[2]]
  if (root[1] > 0) {
    # ancestral state is 00
    # check is trait 1-1 can be born
    if (y1[2] != b1$node) {
      if (!is_parent(t, b1$node, y1[2])) {
        # invalid root state
        root_prob[1] <- 0
      }
    }
    # check if trait 2-1 can be born
    if (y2[2] != b2$node) {
      if (!is_parent(t, b2$node, y2[2])) {
        # invalid root state
        root_prob[1] <- 0
      }
    }
  }
  if (root[2] > 0) {
    # ancestral state is 01
    # check is trait 1-1 can be born
    if (y1[2] != b1$node) {
      if (!is_parent(t, b1$node, y1[2])) {
        # invalid root state
        root_prob[2] <- 0
      }
    }
    # check if trait 2-0 can be born
    if (y2[1] != b2$node) {
      if (!is_parent(t, b2$node, y2[1])) {
        # invalid root state
        root_prob[2] <- 0
      }
    }
  }
  if (root[3] > 0) {
    # ancestral state is 10
    # check is trait 1-0 can be born
    if (y1[1] != b1$node) {
      if (!is_parent(t, b1$node, y1[1])) {
        # invalid root state
        root_prob[3] <- 0
      }
    }
    # check if trait 2-1 can be born
    if (y2[2] != b2$node) {
      if (!is_parent(t, b2$node, y2[2])) {
        # invalid root state
        root_prob[3] <- 0
      }
    }
  }
  if (root[4] > 0) {
    # ancestral state is 11
    # check is trait 1-0 can be born
    if (y1[1] != b1$node) {
      if (!is_parent(t, b1$node, y1[1])) {
        # invalid root state
        root_prob[4] <- 0
      }
    }
    # check if trait 2-0 can be born
    if (y2[1] != b2$node) {
      if (!is_parent(t, b2$node, y2[1])) {
        # invalid root state
        root_prob[4] <- 0
      }
    }
  }

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
