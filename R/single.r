fit_sinba_single <- function(
    tree, data, model,
    pi_x, root,
    opts) {
  if (!inherits(tree, "phylo")) {
    stop("fit_sinba: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  mQ <- model$model
  k <- max(mQ) + 1

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$states)
  }
  if (length(pi_x) != length(model$states)) {
    stop("fit_sinba: invalid pi_x: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 1)
  cond <- set_conditionals(t, et, model)

  max_rate <- maximum_transition_rate / max(t$ages)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    max_len <- t$age[yn]
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p[2:length(p)] > max_rate)) {
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
        r, pi_x, root
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

  root_states <- c("0", "1")

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

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    birth = birth,
    root = root,
    pi_x = pi_x,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

fit_sinba_single_fixed_birth <- function(
    tree, data, birth, model,
    pi_x, root,
    opts) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_births: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  mQ <- model$model
  k <- max(mQ)

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$states)
  }
  if (length(pi_x) != length(model$states)) {
    stop("fit_fixed_births: invalid pi_x: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 1)
  cond <- set_conditionals(t, et, model)

  max_rate <- maximum_transition_rate / max(t$ages)

  if (length(birth) < 1) {
    stop("fit_fixed_births: at least ine birth event is required")
  }
  birth <- birth[[1]]
  if (is.null(birth$node)) {
    stop(sprintf(
      "fit_fixed_births: expecting field `node` in `birth` event"
    ))
  }
  if (is.null(birth$age)) {
    stop(sprintf(
      "fit_fixed_births: expecting field `age` in `birth` event"
    ))
  }

  # closure for the likelihood function
  like_func <- function(yn, r) {
    max_len <- t$age[yn]
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (any(p < 0)) {
        return(Inf)
      }
      if (any(p > max_rate)) {
        return(Inf)
      }
      Q <- from_model_to_Q(mQ, p)
      lk <- sinba_single_like(
        t, Q, model, birth, xt, cond,
        r, pi_x, root
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

  root_states <- c("0", "1")

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
    yn <- youngest[[1]]
    if (!is_valid_birth(t, birth$node, yn)) {
      next
    }

    fn <- like_func(yn, r)
    par <- runif(max(mQ))
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

  # if no birth is compatible with the birth event
  # the likelihood is 0
  if (is.infinite(res$objective)) {
    obj <- list(
      logLik = -Inf,
      k = k,
      model = model,
      Q = matrix(nrow = 4, ncol = 4),
      birth = birth,
      root = root,
      pi_x = pi_x,
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba"
    return(obj)
  }

  q <- from_model_to_Q(mQ, res$solution)
  q <- normalize_Q(q)

  obj <- list(
    logLik = -res$objective,
    k = k,
    model = model,
    Q = q,
    birth = birth,
    root = root,
    pi_x = pi_x,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

fit_sinba_single_fixed_matrix <- function(
    tree, data, rate_mat, model,
    pi_x, root, opts) {
  if (!inherits(tree, "phylo")) {
    stop("fit_fixed_matrix: `tree` must be an object of class \"phylo\".")
  }
  t <- phylo_to_sinba(tree)

  k <- 1

  if (is.null(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a matrix")
  }
  if (nrow(rate_mat) != ncol(rate_mat)) {
    stop("fit_fixed_matrix: `rate_mat` must be a square matrix")
  }
  if (nrow(rate_mat) != nrow(model$model)) {
    stop(
      "fit_fixed_matrix: `rate_mat` must have the same size as the `model`"
    )
  }

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$states)
  }
  if (length(pi_x) != length(model$states)) {
    stop("fit_fixed_matrix: invalid pi_x: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 1)
  cond <- set_conditionals(t, et, model)

  # closure for the likelihood function
  like_func <- function(yn, r) {
    max_len <- t$age[yn]
    xt <- tree_to_cpp(t)

    return(function(p) {
      if (p[1] < 0) {
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

      lk <- sinba_single_like(
        t, rate_mat, model, birth, xt, cond,
        r, pi_x, root
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

  root_states <- c("0", "1")

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
    yn <- youngest[[1]]
    fn <- like_func(yn, r)
    par <- runif(1, max = t$age[yn])
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

  q <- normalize_Q(rate_mat)
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
    root = root,
    pi_x = pi_x,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

fixed_sinba_single <- function(
    tree, data, rate_mat, birth,
    model, pi_x, root) {
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
    stop(
      "fixed_sinba: `rate_mat` must have the same size as the `model`"
    )
  }

  if (length(pi_x) == 0) {
    pi_x <- default_pi_vector(model$states)
  }
  if (length(pi_x) != length(model$states)) {
    stop("fixed_sinba: invalid pi_x: size different to number of states")
  }
  if (sum(pi_x) != 0) {
    pi_x <- pi_x / sum(pi_x)
  }

  if ((is.null(root)) || (sum(root) == 0)) {
    root <- rep(0, length(model$states))
  }
  if (sum(root) != 0) {
    root <- root / sum(root)
  }

  et <- encode_traits(t, data, 1)
  cond <- set_conditionals(t, et, model)

  if (length(birth) < 1) {
    stop("fixed_sinba: at least ine birth event is required")
  }
  birth <- birth[[1]]
  if (is.null(birth$node)) {
    stop(sprintf(
      "fixed_sinba: expecting field `node` in `birth` event"
    ))
  }
  if (is.null(birth$age)) {
    stop(sprintf(
      "fixed_sinba: expecting field `age` in `birth` event"
    ))
  }

  root_states <- c("0", "1")

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
    yn <- youngest[[1]]
    if (!is_valid_birth(t, birth$node, yn)) {
      next
    }

    xt <- tree_to_cpp(t)
    lk <- sinba_single_like(
      t, rate_mat, model, birth, xt, cond,
      r, pi_x, root
    )

    if (lk > res$objective) {
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
      root = root,
      pi_x = pi_x,
      data = data,
      tree = tree
    )
    class(obj) <- "fit_sinba"
    return(obj)
  }

  obj <- list(
    logLik = res$objective,
    k = 0,
    model = model,
    Q = normalize_Q(rate_mat),
    birth = birth,
    root = root,
    pi_x = pi_x,
    data = data,
    tree = tree
  )
  class(obj) <- "fit_sinba"
  return(obj)
}

# sinba_single calculates the likelihood of a single trait
# under the sinba model.
sinba_single_like <- function(
    t, Q, model, birth,
    xt, cond, root, pi_x, pi_root) {
  l <- sinba_single_cond(t, Q, model, birth, xt, cond, root, pi_x)

  return(add_root_prior(l[t$root_id, ], pi_root))
}

sinba_single_cond <- function(
    t, Q, model, birth,
    xt, cond, root, pi_x) {
  # make sure the Q matrix is valid
  Q <- normalize_Q(Q)
  m_PI_root <- build_pi_matrix(model, "x", state_vector(model, "x"))

  st <- as.integer(active_status(t, birth$node, length(t$parent) + 1))

  l <- sinba_conditionals(
    xt$parent, xt$nodes, st, xt$branch,
    cond,
    birth$age, 0,
    Q, Q,
    pi_x, pi_x,
    m_PI_root, m_PI_root
  )
  return(l)
}
