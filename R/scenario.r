# scenario returns the trait activation scenario
# given the root state
# and the first active trait,
scenario <- function(r, tr) {
  if (r == 1) {
    # root 00
    if (tr == 1) {
      # first trait active: 00 <-> 10
      return("13")
    }
    # second trait active: 00 <-> 01
    return("12")
  }
  if (r == 2) {
    # root 01
    if (tr == 1) {
      # first trait active: 01 <-> 11
      return("24")
    }
    # second trait active: 01 <-> 00
    return("12")
  }
  if (r == 3) {
    # root 10
    if (tr == 1) {
      # first trait active: 10 <-> 00
      return("13")
    }
    # second trait active: 10 <-> 11
    return("34")
  }
  # root 11
  if (tr == 1) {
    # first trait active: 11 <-> 01
    return("24")
  }
  # second trait active: 11 <-> 10
  return("34")
}

# birth_events returns the possible birth events
# for the given birth nodes.
# This is for the case in which both nodes
# are not in the same path.
birth_events <- function(t, youngest) {
  ev <- list()
  y1 <- youngest[[1]]
  y2 <- youngest[[2]]

  # the nodes are equal
  if (y1 == y2) {
    ev[[1]] <- c(y1, y2)
    return(ev)
  }

  # the nodes are on the same path towards the root
  if (is_parent(t, y1, y2) || is_parent(t, y2, y1)) {
    ev[[1]] <- c(y1, y2)
    return(ev)
  }

  # the MRCA of both nodes
  p1 <- path_to_node(t, y1)
  p2 <- path_to_node(t, y2)
  x <- t$root_id
  for (i in seq_len(length(p1))) {
    if (i > length(p2)) {
      break
    }
    if (p1[i] != p2[i]) {
      break
    }
    x <- p1[i]
  }
  ev[[1]] <- c(x, y2)
  ev[[2]] <- c(y1, x)

  return(ev)
}

# is_valid_birth returns true is a birth event
# is compatible.
is_valid_birth <- function(t, ev, yn) {
  if (ev == yn) {
    return(TRUE)
  }
  return(is_parent(t, ev, yn))
}

# set_root_prior sets the priors of the root
# taking into account the birth events.
set_root_prior <- function(t, model, root, births, youngest) {
  obs_prior <- update_root(t, rep(1, 4), births, youngest)
  names(obs_prior) <- c("00", "01", "10", "11")
  prior <- rep(0, length(root))
  for (i in seq_len(length(root))) {
    obs <- model$observed[[model$states[i]]]
    prior[i] <- root[i] * obs_prior[obs]
  }
  return(prior)
}

# set_root_prior_single sets the prior of the root
# taking into account the birth events
# for a single trait.
set_root_prior_single <- function(model, root, obs_prior) {
  names(obs_prior) <- c("0", "1")
  prior <- rep(0, length(root))
  for (i in seq_len(length(root))) {
    obs <- model$observed[[model$states[i]]]
    prior[i] <- root[i] * obs_prior[obs]
  }
  return(prior)
}

# update_root returns the root priors
# taking into account the birth events.
update_root <- function(t, root, births, youngest) {
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
  return(root_prob)
}

# creates the pi matrix for the birth of an event.
# the pi matrix is really the P matrix
# (the probability matrix)
# used in a birth event.
build_pi_matrix <- function(m, trait, root) {
  pi <- matrix(0, nrow = length(m$states), ncol = length(m$states))
  v <- state_vector(m, trait)
  for (i in seq_len(nrow(pi))) {
    pi[i, ] <- v
  }
  if (m$traits == 1) {
    states <- m$states

    # remove states not in the root
    names(root) <- c("0", "1")
    for (i in seq_len(length(states))) {
      from <- m$observed[[states[i]]]
      if (root[from] != 0) {
        next
      }
      pi[i, ] <- 0
    }

    if (m$rates == 1) {
      return(pi)
    }
    # multiple rates model
    for (i in seq_len(length(states))) {
      row_state <- strsplit(states[i], split = "[", fixed = TRUE)[[1]]
      for (j in seq_len(length(states))) {
        col_state <- strsplit(states[j], split = "[", fixed = TRUE)[[1]]
        if (row_state[2] == col_state[2]) {
          # skip if rates are the same
          next
        }
        # we only change inside the same rate
        pi[i, j] <- 0
      }
    }
    return(pi)
  }
  cm <- matrix(0, nrow = 4, ncol = 4)
  if (trait == "x") {
    cm[1, 3] <- 1
    cm[2, 4] <- 1
    cm[3, 1] <- 1
    cm[4, 2] <- 1
  } else {
    cm[2, 1] <- 1
    cm[1, 2] <- 1
    cm[3, 4] <- 1
    cm[4, 3] <- 1
  }
  rownames(cm) <- c("00", "01", "10", "11")
  colnames(cm) <- c("00", "01", "10", "11")
  names(root) <- c("00", "01", "10", "11")

  states <- m$states
  for (i in seq_len(length(states))) {
    from <- m$observed[[states[i]]]
    if (root[from] == 0) {
      pi[i, ] <- 0
      next
    }
    for (j in seq_len(length(states))) {
      if (i == j) {
        # skip diagonal
        next
      }

      # in a rates model we should be in the same rate
      if (m$rates > 1) {
        row_state <- strsplit(states[i], split = "[", fixed = TRUE)[[1]]
        col_state <- strsplit(states[j], split = "[", fixed = TRUE)[[1]]
        if (row_state[2] != col_state[2]) {
          pi[i, j] <- 0
          next
        }
      }
      to <- m$observed[[states[j]]]
      if (cm[from, to] != 0) {
        next
      }
      pi[i, j] <- 0
    }
  }

  return(pi)
}

active_ancestor_vector <- function(sc) {
  if (sc == "12") {
    return(c(1, 1, 0, 0))
  }
  if (sc == "13") {
    return(c(1, 0, 1, 0))
  }
  if (sc == "24") {
    return(c(0, 1, 0, 1))
  }
  if (sc == "34") {
    return(c(0, 0, 1, 1))
  }
  stop(sprintf("invalid scenario: '%s'", sc))
}

default_pi_vector <- function(states) {
  pi <- rep(0, length(states))
  for (i in seq_len(length(states))) {
    s <- strsplit(states[i], split = "", fixed = TRUE)[[1]]
    if (length(s) < 2) {
      next
    }
    if (s[2] == "1") {
      pi[i] <- 1
    }
  }
  return(pi)
}

set_pi_matrix <- function(m_PI, v_pi) {
  m <- matrix(0, nrow = nrow(m_PI), ncol = ncol(m_PI))
  for (i in seq_len(nrow(m))) {
    for (j in seq_len(ncol(m))) {
      if (m_PI[i, j] == 0) {
        next
      }
      m[i, j] <- v_pi[m_PI[i, j]]
    }
  }
  return(m)
}
