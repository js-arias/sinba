# youngest_birth_event returns the youngest node
# for the birth of each state.
youngest_birth_event <- function(t, cond) {
  youngest <- list()

  # a single trait
  if (ncol(cond) == 2) {
    nodes <- c()
    for (i in seq_len(ncol(cond))) {
      ns <- rep(FALSE, nrow(cond))
      for (n in seq_len(length(t$tip))) {
        # set tips
        if (cond[n, i] == 0) {
          ns[n] <- TRUE
        }
      }
      nodes <- c(nodes, dollo(t, ns))
    }
    youngest[[1]] <- nodes
    return(youngest)
  }

  # two traits
  # trait 1
  nodes <- c()
  # state 0
  ns <- rep(FALSE, nrow(cond))
  for (n in seq_len(length(t$tip))) {
    # set tips
    if (cond[n, 1] == 0 || cond[n, 2] == 0) {
      ns[n] <- TRUE
    }
  }
  nodes <- c(nodes, dollo(t, ns))
  # state 1
  ns <- rep(FALSE, nrow(cond))
  for (n in seq_len(length(t$tip))) {
    # set tips
    if (cond[n, 3] == 0 || cond[n, 4] == 0) {
      ns[n] <- TRUE
    }
  }
  nodes <- c(nodes, dollo(t, ns))
  youngest[[1]] <- nodes

  # trait 2
  nodes <- c()
  # state 0
  ns <- rep(FALSE, nrow(cond))
  for (n in seq_len(length(t$tip))) {
    # set tips
    if (cond[n, 1] == 0 || cond[n, 3] == 0) {
      ns[n] <- TRUE
    }
  }
  nodes <- c(nodes, dollo(t, ns))
  # state 1
  ns <- rep(FALSE, nrow(cond))
  for (n in seq_len(length(t$tip))) {
    # set tips
    if (cond[n, 2] == 0 || cond[n, 4] == 0) {
      ns[n] <- TRUE
    }
  }
  nodes <- c(nodes, dollo(t, ns))
  youngest[[2]] <- nodes

  return(youngest)
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

  return(dollo_uppass(t, t$root_id, ns))
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
