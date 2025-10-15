# encode_traits set up the observed traits in the data
# as the observed states.
encode_traits <- function(t, data, traits) {
  if (traits == 1) {
    et <- list()
    for (i in seq_len(length(t$tip))) {
      st <- retrieve_state(data, t$tip[i], 2)
      et[[t$tip[i]]] <- list(
        name = t$tip[i],
        state = st
      )
    }
    return(et)
  }
  et <- list()
  for (i in seq_len(length(t$tip))) {
    x <- retrieve_state(data, t$tip[i], 2)
    y <- retrieve_state(data, t$tip[i], 3)
    st <- sprintf("%s%s", x, y)
    et[[t$tip[i]]] <- list(
      name = t$tip[i],
      state = st
    )
  }
  return(et)
}

retrieve_state <- function(data, tip, trait) {
  r <- which(data[, 1] == tip)
  if (length(r) == 0) {
    # unobserved tip
    return("?")
  }
  states <- c(FALSE, FALSE)
  for (i in seq_len(length(r))) {
    if (data[r[i], trait] == 0) {
      states[1] <- TRUE
    }
    if (data[r[i], trait] == 1) {
      states[2] <- TRUE
    }
  }
  if (all(states)) {
    # polymorphic tip
    return("p")
  }
  if (!states[1] && !states[2]) {
    # unknown state is treated as unknown
    if (length(r) == 0) {
      # unobserved tip
      return("?")
    }
  }
  if (states[1]) {
    return("0")
  }
  return("1")
}


# set_conditionals set up the initial log conditionals
# (i.e., the tips)
# of a tree.
# Internal nodes are set as 0
# (i.e., with probability 1).
set_conditionals <- function(t, enc, model) {
  cond <- matrix(0, nrow = length(t$parent), ncol = length(model$states))
  for (i in seq_len(length(t$tip))) {
    cond[i, ] <- -Inf
    obs <- enc[[t$tip[i]]]$state
    if (obs == "?" || obs == "p") {
      obs <- c("0", "1")
    } else if (obs == "?0" || obs == "p0") {
      obs <- c("00", "10")
    } else if (obs == "?1" || obs == "p1") {
      obs <- c("01", "11")
    } else if (obs == "0?" || obs == "0p") {
      obs <- c("00", "01")
    } else if (obs == "1?" || obs == "1p") {
      obs <- c("10", "11")
    } else if (obs == "??" || obs == "p?" || obs == "?p" || obs == "pp") {
      obs <- c("00", "01", "10", "11")
    }
    for (j in seq_len(length(model$states))) {
      x <- model$observed[model$states[j]]
      if (any(obs == x)) {
        cond[i, j] <- 0
      }
    }
  }
  return(cond)
}

# init_conditionals set up the initial log conditionals
# (i.e., the tips)
# of a tree.
# Internal nodes are set as 0
# (i.e., with probability 1).
init_conditionals <- function(t, data, traits) {
  if (traits == 1) {
    # conditionals of a single trait
    cond <- matrix(0, nrow = length(t$parent), ncol = 2)
    for (i in seq_len(length(t$tip))) {
      cond[i, ] <- -Inf
      r <- which(data[, 1] == t$tip[i])
      if (length(r) == 0) {
        stop(sprintf(
          "init_conditional: tip '%s' not found in `data`",
          t$tip[i]
        ))
      }
      if (length(r) != 1) {
        stop(sprintf(
          "init_conditional: multiple entries for '%s' in `data`",
          t$tip[i]
        ))
      }

      if (data[r, 2] == 0) {
        cond[i, 1] <- 0
      } else if (data[r, 2] == 1) {
        cond[i, 2] <- 0
      } else {
        stop(sprintf(
          "init_conditional: unknown state %d for '%s' in `data`",
          data[r, 2], t$tip[i]
        ))
      }
    }
    return(cond)
  }

  # conditional of two traits
  cond <- matrix(0, nrow = length(t$parent), ncol = 4)
  for (i in seq_len(length(t$tip))) {
    cond[i, ] <- -Inf
    r <- which(data[, 1] == t$tip[i])
    if (length(r) == 0) {
      stop(sprintf(
        "init_conditional: tip '%s' not found in `data`",
        t$tip[i]
      ))
    }
    if (length(r) != 1) {
      stop(sprintf(
        "init_conditional: multiple entries for '%s' in `data`",
        t$tip[i]
      ))
    }

    # do all combinations
    if (data[r, 2] == 0 && data[r, 3] == 0) {
      cond[i, 1] <- 0
    } else if (data[r, 2] == 0 && data[r, 3] == 1) {
      cond[i, 2] <- 0
    } else if (data[r, 2] == 1 && data[r, 3] == 0) {
      cond[i, 3] <- 0
    } else if (data[r, 2] == 1 && data[r, 3] == 1) {
      cond[i, 4] <- 0
    } else {
      stop(sprintf(
        "init_conditional: unknown combination %d-%d for '%s' in `data`",
        data[r, 2], data[r, 3], t$tip[i]
      ))
    }
  }
  return(cond)
}

# observed returns the state combinations of two traits.
observed <- function(t, data) {
  obs <- rep(FALSE, 4)
  for (i in seq_len(length(t$tip))) {
    r <- which(data[, 1] == t$tip[i])
    # do all combinations
    if (data[r, 2] == 0 && data[r, 3] == 0) {
      obs[1] <- TRUE
    } else if (data[r, 2] == 0 && data[r, 3] == 1) {
      obs[2] <- TRUE
    } else if (data[r, 2] == 1 && data[r, 3] == 0) {
      obs[3] <- TRUE
    } else if (data[r, 2] == 1 && data[r, 3] == 1) {
      obs[4] <- TRUE
    }
  }
  return(obs)
}
