#' @import stats

#' @export
#' @title
#' Stochastic Map For the Sinba Model
#'
#' @description
#' `map_sinba()` build one or more stochastic mappings
#' using a reconstruction from the Sinba model.
#'
#' @param fitted An object of the type 'fit_sinba'
#'   (i.e., a reconstruction using the Sinba model).
#' @param n Number of stochastic maps.
#'   Default is 100.
map_sinba <- function(fitted, n = 100) {
  if (!inherits(fitted, "fit_sinba")) {
    stop("map_sinba: `fitted` must be an object of class \"fit_sinba\".")
  }
  if (is.infinite(fitted$logLik)) {
    stop("map_sinba: impossible reconstruction (logLik == -Inf).")
  }

  model <- fitted$model
  t <- phylo_to_sinba(fitted$tree)
  Q <- fitted$Q
  sc <- fitted$sc
  pi_x <- fitted$pi_x
  pi_y <- fitted$pi_y
  et <- encode_traits(t, fitted$data, 2)
  cond <- set_conditionals(t, et, fitted$model)
  births <- fitted$births

  m_PI_act <- matrix()
  m_PI_semi <- matrix()

  root_states <- c("00", "01", "10", "11")
  root <- fitted$root
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

    root_vector <- rep(0, 4)
    root_vector[r] <- 1
    if ((sc == "12") || (sc == "34")) {
      anc_vector <- active_ancestor_vector(sc)
      m_PI_act <- set_pi_matrix(build_pi_matrix(model, "x", anc_vector), pi_x)
      m_PI_semi <- set_pi_matrix(build_pi_matrix(model, "y", root_vector), pi_y)
    } else {
      anc_vector <- active_ancestor_vector(sc)
      m_PI_act <- set_pi_matrix(build_pi_matrix(model, "y", anc_vector), pi_y)
      m_PI_semi <- set_pi_matrix(build_pi_matrix(model, "x", root_vector), pi_x)
    }

    xt <- tree_to_cpp(t)
    sc <- sinba_cond(t, Q, model, births, xt, cond, r, pi_x, pi_y, FALSE)

    s_names <- fitted$model$states
    edges <- paste(fitted$tree$edge[, 1], ",", fitted$tree$edge[, 2], sep = "")
    l <- sc$l
    st <- sc$st
    ev <- sc$ev
    if (ev$second$n == t$root_id) {
      mx <- max(l[t$root_id, ])
      ev$second$cond <- exp(l[t$root_id, ] - mx)
    } else {
      br_len <- xt$branch[n] - ev$second$age
      ev$second$cond <- birth_conditional(
        br_len, ev$second$Q, l[ev$second$n, ],
        m_PI_act
      )
    }
    if (ev$first$n == t$root_id) {
      mx <- max(l[t$t_root_id, ])
      ev$first$cond <- exp(l[t$root_id, ] - mx)
    } else {
      br_len <- xt$branch[n] - ev$first$age
      nc <- l[ev$first$n, ]
      if (ev$first$n == ev$second$n) {
        br_len <- ev$second$age - ev$first$age
        mc <- ev$second$cond
        nc <- c()
        for (i in seq_len(nrow(mc))) {
          nc <- c(nc, sum(mc[i, ]))
        }
        nc <- log(nc)
      }
      ev$first$cond <- birth_conditional(
        br_len, ev$first$Q, nc,
        m_PI_semi
      )
    }

    maps <- list()
    for (i in seq_len(n)) {
      sm <- evolve_map_sinba(
        t, l, ev$first$cond, ev$second$cond, st,
        ev$first$age, ev$second$age,
        ev$first$Q, ev$second$Q,
        ev$first$cond, ev$second$cond,
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
}

evolve_map_sinba <- function(
    t, cond, cond_b1, cond_b2,
    st,
    age_1, age_2,
    first_Q, second_Q,
    semi_P, act_P,
    model, root) {
  root_Q <- matrix(0, nrow = nrow(first_Q), ncol = ncol(first_Q))
  root_state <- pick_root_state(root, cond[t$root, ])
  states <- rep(root_state, length(t$parent))
  states[t$root] <- root_state
  sm <- list()
  edges <- matrix(0, nrow = length(t$edge), ncol = ncol(root_Q))
  for (e in seq_len(length(t$edge))) {
    n <- t$edge[e]
    p <- t$parent[n]
    s <- states[p]
    if (st[n] != st[p]) {
      Qs <- list()
      Ps <- list()
      lens <- c()
      stages <- 1
      if (st[p] == 0) {
        Qs[[stages]] <- root_Q
        lens <- c(age_1)
        stages <- 2
        Ps[[stages]] <- semi_P
      }
      Qs[[stages]] <- first_Q
      if (st[n] == 2) {
        lens <- c(lens, age_2)
        stages <- stages + 1
        Qs[[stages]] <- second_Q
        Ps[[stages]] <- act_P
      }
      lens <- c(lens, t$br_len[n])
      h <- pick_history_with_births(Qs, Ps, lens, s, cond[n, ], model$states)
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

pick_history <- function(Q, len, s, end, s_names) {
  end <- exp(end - max(end))
  end <- end / max(end)
  while (TRUE) {
    h <- store_evolution(Q, len, s, s_names)
    if (runif(1) < end[h$end]) {
      return(h)
    }
  }
}

pick_history_with_births <- function(Qs, Ps, lens, s, end, s_names) {
  end <- exp(end - max(end))
  end <- end / max(end)
  while (TRUE) {
    h <- store_evolution_with_births(Qs, Ps, lens, s, s_names)
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

store_evolution_with_births <- function(Qs, Ps, lens, s, s_names) {
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

    # birth event
    nx <- pick_birth_state(Ps[[stage]], s)
    if (nx != s) {
      # trait change at the birth event
      times <- c(times, t)
      evs <- c(evs, s)
      s <- nx
    }
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
      # birth event
      stage <- stage + 1
      Q <- Qs[[stage]]
      t <- len
      nx <- pick_birth_state(Ps[[stage]], s)
      if (nx != s) {
        # trait change at the birth event
        times <- c(times, t)
        evs <- c(evs, s)
        s <- nx
      }
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

pick_birth_state <- function(P, s) {
  mx <- max(P[s, ])
  while (TRUE) {
    nx <- sample(seq_len(nrow(P)), size = 1)
    if (runif(1) < P[s, nx] / mx) {
      return(nx)
    }
  }
}

#' @export
#' @title
#' Stochastic Map For the Pagel´s Model
#'
#' @description
#' `map_pagel()` build one or more stochastic mappings
#' using a reconstruction from the Pagel's model.
#'
#' @param fitted An object of the type 'fit_mk'
#'   (i.e., a reconstruction produced with `fit_pagel()` function).
#' @param n Number of stochastic maps.
#'   Default is 100.
map_pagel <- function(fitted, n = 100) {
  if (!inherits(fitted, "fit_mk")) {
    stop("map_pagel: `fitted` must be an object of class \"fit_mk\".")
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
  root_state <- pick_root_state(root, cond[t$root, ])
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

pick_root_state <- function(root, end) {
  end <- exp(end - max(end))
  end <- end / sum(end)
  if (sum(root) == 0) {
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
