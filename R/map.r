#' @import stats

#' @export
#' @title
#' Stochastic Map For the Sinba Model
#'
#' @description
#' `map_sinba()` build one or more stochastic mappings
#' using a reconstruction from the Sinba model
#' or by reading a tree,
#' a data set,
#' and a model,
#' and fitting the model on the tree
#' for each stochastic map.
#'
#' @param fitted An object of the type 'fit_sinba'
#'   (i.e., a reconstruction using the Sinba model).
#' @param n Number of stochastic maps.
#'   Default is 100.
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
map_sinba <- function(
    fitted = NULL, n = 100,
    tree = NULL, data = NULL, model = NULL, opts = NULL) {
  if (!is.null(fitted)) {
    return(map_sinba_fitted(fitted, n))
  }
  if (is.null(tree) || is.null(data)) {
    stop("map_sinba: undefined data")
  }

  maps <- list()
  for (i in seq_len(n)) {
    f <- fit_sinba(tree, data, model = model, opts = opts)
    sm <- map_sinba_fitted(f, n = 1)
    maps[[i]] <- sm[[1]]
  }
  class(maps) <- c("multiSimmap", "multiPhylo")
  return(maps)
}

map_sinba_fitted <- function(fitted, n) {
  if (!inherits(fitted, "fit_sinba")) {
    stop("map_sinba: `fitted` must be an object of class \"fit_sinba\".")
  }
  if (is.infinite(fitted$logLik)) {
    stop("map_sinba: impossible reconstruction (logLik == -Inf).")
  }

  t <- phylo_to_sinba(fitted$tree)
  et <- encode_traits(t, fitted$data, 2)
  cond <- set_conditionals(t, et, fitted$model)
  b1 <- fitted$births[[1]]
  b2 <- fitted$births[[2]]

  root <- fitted$root
  root_id <- 0
  root_states <- c("00", "01", "10", "11")
  for (i in seq_len(length(root_states))) {
    if (root_states[i] == root) {
      root_id <- i
    }
  }
  Q <- fitted$Q
  model <- fitted$model

  xt <- tree_to_cpp(t)
  root_Q <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))

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
    sc <- scenario(root_id, 1)
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
    sc <- scenario(root_id, 2)
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

  s_names <- fitted$model$states
  edges <- paste(fitted$tree$edge[, 1], ",", fitted$tree$edge[, 2], sep = "")

  maps <- list()
  for (i in seq_len(n)) {
    sm <- evolve_map_states(
      t, l, st,
      ev$first$age, ev$second$age,
      root_Q, ev$first$Q, ev$second$Q,
      model, fitted$root
    )
    mtree <- fitted$tree
    mtree$maps <- sm$maps
    mtree$mapped.edge <- sm$edges
    rownames(mtree$mapped.edge) <- edges
    colnames(mtree$mapped.edge) <- s_names
    class(mtree) <- c("simmap", "phylo")
    attr(mtree, "map.order") <- "right-to-left"
    maps[[i]] <- mtree
  }
  class(maps) <- c("multiSimmap", "multiPhylo")
  return(maps)
}

evolve_map_states <- function(
    t, cond, st,
    age_1, age_2,
    root_Q, first_Q, second_Q,
    model, root) {
  root_state <- pick_root_state(model, root, cond[t$root, ])
  states <- rep(0, length(t$parent))
  states[t$root] <- root_state
  sm <- list()
  edges <- matrix(0, nrow = length(t$edge), ncol = ncol(root_Q))
  for (e in seq_len(length(t$edge))) {
    n <- t$edge[e]
    p <- t$parent[n]
    s <- states[p]
    if (st[n] != st[p]) {
      Qs <- list()
      lens <- c()
      stages <- 1
      if (st[p] == 0) {
        Qs[[stages]] <- root_Q
        lens <- c(age_1)
        stages <- 2
      }
      Qs[[stages]] <- first_Q
      if (st[n] == 2) {
        lens <- c(lens, age_2)
        stages <- stages + 1
        Qs[[stages]] <- second_Q
      }
      lens <- c(lens, t$br_len[n])
      h <- pick_history_with_births(Qs, lens, s, cond[n, ], model$states)
      states[n] <- h$end
      sm[[e]] <- h$evs
      edges[e, ] <- h$cm
      next
    }
    Q <- root_Q
    if (st[n] == 2) {
      Q <- second_Q
    }
    if (st[n] == 1) {
      Q <- first_Q
    }
    h <- pick_history(Q, t$br_len[n], s, cond[n, ], model$states)
    states[n] <- h$end
    sm[[e]] <- h$evs
    edges[e, ] <- h$cm
  }
  return(list(
    maps = sm,
    edges = edges
  ))
}

pick_root_state <- function(model, root, end) {
  # remove non root states
  for (i in seq_len(length(end))) {
    s <- model$states[i]
    if (model$observed[[s]] != root) {
      end[i] <- -Inf
    }
  }

  end <- exp(end - max(end))
  end <- end / sum(end)
  # apply FitzJohn prior
  end <- end * end
  end <- end / max(end)

  while (TRUE) {
    s <- sample(seq_len(length(end)), 1)
    if (runif(1) < end[s]) {
      return(s)
    }
  }
}

pick_history <- function(Q, len, s, end, s_names) {
  end <- exp(end - max(end))
  end <- end / sum(end)
  while (TRUE) {
    h <- store_evolution(Q, len, s, s_names)
    if (runif(1) < end[h$end]) {
      return(h)
    }
  }
}

pick_history_with_births <- function(Qs, lens, s, end, s_names) {
  end <- exp(end - max(end))
  end <- end / sum(end)
  while (TRUE) {
    h <- store_evolution_with_births(Qs, lens, s, s_names)
    if (runif(1) < end[h$end]) {
      return(h)
    }
  }
}

# make a simulation of the evolution
# and store the event times
# and states at each event.
store_evolution <- function(Q, len, s, s_names) {
  t <- 0
  evs <- c()
  times <- c()
  cm <- rep(0, nrow(Q))
  while (Q[s, s] != 0) {
    dt <- stats::rexp(1, rate = -Q[s, s])
    if (t + dt >= len) {
      break
    }
    t <- t + dt
    times <- c(times, dt)
    evs <- c(evs, s)
    cm[s] <- cm[s] + dt
    s <- pick_next_state(Q, s)
  }
  dt <- len - t
  times <- c(times, dt)
  evs <- c(evs, s)
  cm[s] <- cm[s] + dt
  names(times) <- s_names[evs]
  return(list(
    end = s,
    evs = times,
    cm = cm
  ))
}

store_evolution_with_births <- function(Qs, lens, s, s_names) {
  t <- 0
  evs <- c()
  times <- c()

  stage <- 1
  Q <- Qs[[1]]
  len <- lens[1]
  if (Q[s, s] == 0) {
    stage <- 2
    Q <- Qs[[stage]]
    t <- len
    len <- lens[stage]
  }
  while (TRUE) {
    dt <- len - t
    if (Q[s, s] != 0) {
      dt <- stats::rexp(1, rate = -Q[s, s])
    }
    if (t + dt >= len) {
      if (stage == length(lens)) {
        break
      }
      stage <- stage + 1
      Q <- Qs[[stage]]
      t <- len
      len <- lens[stage]
      next
    }
    t <- t + dt
    # in this function
    # event times are stored
    times <- c(times, t)
    evs <- c(evs, s)
    s <- pick_next_state(Q, s)
  }
  times <- c(times, len)

  # from full times to delta times
  prev <- 0
  dts <- c()
  cm <- rep(0, nrow(Q))
  for (i in seq_len(length(times))) {
    dt <- times[i] - prev
    dts <- c(dts, dt)
    cm[evs[i]] <- cm[evs[i]] + dt
    prev <- times[i]
  }

  evs <- c(evs, s)
  names(dts) <- s_names[evs]
  return(list(
    end = s,
    evs = dts,
    cm = cm
  ))
}

pick_next_state <- function(Q, s) {
  while (TRUE) {
    nx <- sample(seq_len(nrow(Q)), size = 1)
    if (nx != s) {
      if (runif(1) < -Q[s, nx] / Q[s, s]) {
        return(nx)
      }
    }
  }
}

#' @export
#' @title
#' Stochastic Map For the Pagel´s Model
#'
#' @description
#' `map_pagel()` build one or more stochastic mappings
#' using a reconstruction from the Pagel's model
#' or by reading a tree,
#' a data set,
#' and a model,
#' and fitting the model on the tree
#' for each stochastic map.
#'
#' @param fitted An object of the type 'fit_sinba'
#'   (i.e., a reconstruction using the Sinba model).
#' @param n Number of stochastic maps.
#'   Default is 100.
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
map_pagel <- function(
    fitted = NULL, n = 100,
    tree = NULL, data = NULL, model = NULL, opts = NULL) {
  if (!is.null(fitted)) {
    return(map_pagel_fitted(fitted, n))
  }
  if (is.null(tree) || is.null(data)) {
    stop("map_pagel: undefined data")
  }
  f <- fit_pagel(tree, data, model = model, opts = opts)
  return(map_pagel_fitted(f, n))
}

map_pagel_fitted <- function(fitted, n) {
  if (!inherits(fitted, "fit_mk")) {
    stop("map_pagel: `fitted` must be an object of class \"fit_sinba\".")
  }

  t <- phylo_to_sinba(fitted$tree)
  et <- encode_traits(t, fitted$data, 2)
  cond <- set_conditionals(t, et, fitted$model)
  root <- fitted$root
  Q <- fitted$Q
  model <- fitted$model

  xt <- tree_to_cpp(t)
  l <- full_conditionals(
    xt$parent, xt$nodes, xt$branch, cond, Q
  )

  s_names <- fitted$model$states
  edges <- paste(fitted$tree$edge[, 1], ",", fitted$tree$edge[, 2], sep = "")

  maps <- list()
  for (i in seq_len(n)) {
    sm <- evolve_map_mk(t, l, Q, model, root)
    mtree <- fitted$tree
    mtree$maps <- sm$maps
    mtree$mapped.edge <- sm$edges
    rownames(mtree$mapped.edge) <- edges
    colnames(mtree$mapped.edge) <- s_names
    class(mtree) <- c("simmap", "phylo")
    attr(mtree, "map.order") <- "right-to-left"
    maps[[i]] <- mtree
  }
  class(maps) <- c("multiSimmap", "multiPhylo")
  return(maps)
}

evolve_map_mk <- function(t, cond, Q, model, root) {
  root_state <- pick_mk_root(root, cond[t$root, ])
  states <- rep(0, length(t$parent))
  states[t$root] <- root_state
  sm <- list()
  edges <- matrix(0, nrow = length(t$edge), ncol = ncol(Q))
  for (e in seq_len(length(t$edge))) {
    n <- t$edge[e]
    p <- t$parent[n]
    s <- states[p]
    h <- pick_history(Q, t$br_len[n], s, cond[n, ], model$states)
    states[n] <- h$end
    sm[[e]] <- h$evs
    edges[e, ] <- h$cm
  }
  return(list(
    maps = sm,
    edges = edges
  ))
}

pick_mk_root <- function(root, end) {
  end <- exp(end - max(end))
  end <- end / sum(end)
  if (is.null(root)) {
    # apply FitzJohn prior
    end <- end * end
  } else {
    end <- end * root
  }

  # scale to the maximum
  end <- end / max(end)

  while (TRUE) {
    s <- sample(seq_len(length(end)), 1)
    if (runif(1) < end[s]) {
      return(s)
    }
  }
}
