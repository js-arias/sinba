# youngest_birth_event returns the youngest node
# for the birth of a given root state.
youngest_birth_event <- function(t, enc, root) {
  youngest <- list()

  # single trait
  if (root == "0" || root == "1") {
    state <- "1"
    if (root == "1") {
      state <- "0"
    }
    youngest[[1]] <- birth_node_states(t, enc, state)
    return(youngest)
  }

  # two traits
  for (i in c(1:2)) {
    root_state <- substr(root, i, i)
    state <- "1"
    if (root_state == "1") {
      state <- "0"
    }
    et <- list()
    for (j in seq_len(length(t$tip))) {
      tp <- t$tip[j]
      et[[tp]] <- list(
        name = tp,
        state = substr(enc[[tp]]$state, i, i)
      )
    }
    youngest[[i]] <- birth_node_states(t, et, state)
  }
  return(youngest)
}

birth_node_states <- function(t, enc, state) {
  ns <- rep(FALSE, length(t$parent))
  for (j in seq_len(length(t$tip))) {
    obs <- enc[[t$tip[j]]]$state
    if (obs == state || obs == "p") {
      ns[j] <- TRUE
    }
  }
  return(dollo(t, ns))
}

# dollo detects the youngest node
# in which an state can be born.
dollo <- function(t, ns) {
  # dollo down pass
  for (i in seq(from = length(t$edge), to = 1)) {
    n <- t$edge[i]
    if (n == t$root_id) {
      # skip root
      next
    }
    # update ancestor
    a <- t$parent[n]
    ns[a] <- ns[a] | ns[n]
  }

  dn <- dollo_uppass(t, t$root_id, ns)
  if (dn == 0) {
    dn <- t$root_id
  }
  return(dn)
}

# dollo_uppass returns the first node
# in which the state born.
dollo_uppass <- function(t, n, ns) {
  # the node does not have the state
  if (!ns[n]) {
    return(0)
  }

  # the node is a tip
  if (n < t$root_id) {
    return(0)
  }

  desc <- which(t$parent == n)
  ws <- 0
  for (i in desc) {
    if (ns[i]) {
      ws <- ws + 1
    }
  }
  if (ws > 1) {
    # more than one descendant has the state
    # so the state in the node is present
    return(n)
  }

  for (i in desc) {
    x <- dollo_uppass(t, i, ns)
    if (x > 0) {
      return(x)
    }
  }
  return(0)
}
