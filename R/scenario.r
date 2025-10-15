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
# for the given trait birth nodes.
birth_events <- function(t, youngest) {
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
  return(ev)
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
