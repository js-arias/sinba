# phylo_node_age returns the age of a node
# in an ape tree
# (an object of class "phylo").
phylo_node_age <- function(phy, n) {
  s <- 0
  while (TRUE) {
    e <- which(phy$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    s <- s + phy$edge.length[e]
    n <- phy$edge[e, 1]
  }
  return(s)
}

# phylo_to_sinba transforms an ape tree
# (an object of class "phylo")
# into a tree used for most internal functions
# in sinba.
# The reason is speed.
phylo_to_sinba <- function(tree) {
  parent <- rep(0, length(tree$tip.label) + tree$Nnode)
  br_len <- rep(0, length(tree$tip.label) + tree$Nnode)
  ages <- rep(0, length(tree$tip.label) + tree$Nnode)
  for (i in seq_len(nrow(tree$edge))) {
    n <- tree$edge[i, 2]
    parent[n] <- tree$edge[i, 1]
    br_len[n] <- tree$edge.length[i]
    ages[n] <- ages[parent[n]] + br_len[n]
  }
  return(list(
    root_id = length(tree$tip.label) + 1,
    parent = parent,
    edge = tree$edge[, 2],
    br_len = br_len,
    ages = ages,
    tip = tree$tip.label
  ))
}

# tree_to_cpp transforms a tree into a set of vectors
# to be used in a C++ function.
tree_to_cpp <- function(t) {
  return(list(
    parent = as.integer(t$parent - 1),
    nodes = as.integer(t$edge - 1),
    branch = t$br_len
  ))
}

# is_parent returns true if node p is parent of n.
is_parent <- function(t, p, n) {
  if (p == n) {
    # a node cannot be parent of itself
    return(FALSE)
  }
  while (n > 0) {
    a <- t$parent[n]
    if (a == p) {
      return(TRUE)
    }
    n <- a
  }
  return(FALSE)
}

# node_size returns the number of terminals
# for a given node.
node_size <- function(t, n) {
  if (n < t$root_id) {
    # a tip
    return(1)
  }
  sum <- 0
  for (x in seq_len(length(t$tip))) {
    if (is_parent(t, n, x)) {
      sum <- sum + 1
    }
  }
  return(sum)
}

# get_node_by_len returns a node
# at a given "age" l (length from the root),
# in the path between n and the root.
get_node_by_len <- function(t, l, n) {
  if (l > t$age[n]) {
    return(0)
  }
  while (n > 0) {
    a <- t$age[n] - t$br_len[n]
    if (a < l) {
      return(n)
    }
    n <- t$parent[n]
  }
  return(t$root_id)
}

# path_to_node returns the path from the root
# to a given node.
path_to_node <- function(t, n) {
  x <- c()
  while (TRUE) {
    x <- c(n, x)
    if (n == t$root_id) {
      break()
    }
    n <- t$parent[n]
  }
  return(x)
}


# prob_birth returns the probability of a birth event.
prob_birth <- function(t) {
  return(max(t$ages) / sum(t$br_len))
}

# active_status returns a vector with the active status
# of each node.
active_status <- function(t, n1, n2) {
  st <- rep(0, length(t$parent))
  for (i in seq_len(length(st))) {
    if (i == n2) {
      st[i] <- 2
      next
    }
    if (i == n1) {
      st[i] <- 1
      next
    }
    if (is_parent(t, n2, i)) {
      st[i] <- 2
      next
    }
    if (is_parent(t, n1, i)) {
      st[i] <- 1
    }
  }
  return(st)
}
