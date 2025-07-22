# phylo_to_sinba transforms an ape tree
# (an object of class "phylo")
# into a tree used for most internal functions
# in sinba.
# The reason is speed.
phylo_to_sinba <- function(tree) {
  parent <- rep(0, length(tree$tip.label) + tree$Nnode)
  br_len <- rep(0, length(tree$tip.label) + tree$Nnode)
  for (i in seq_len(nrow(tree$edge))) {
    parent[tree$edge[i, 2]] <- tree$edge[i, 1]
    br_len[tree$edge[i, 2]] <- tree$edge.length[i]
  }
  return(list(
    root_id = length(tree$tip.label) + 1,
    parent = parent,
    edge = tree$edge[, 2],
    br_len = br_len,
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

# length_to_root returns the total length from a node
# towards the root.
length_to_root <- function(t, n) {
  s <- 0
  while (n > 0) {
    s <- s + t$br_len[n]
    n <- t$parent[n]
  }
  return(s)
}

# prob_birth returns the probability of a birth event.
prob_birth <- function(t) {
  max <- 0
  for (i in seq_len(length(t$tip))) {
    l <- length_to_root(t, i)
    if (l > max) {
      max <- l
    }
  }
  return(max / sum(t$br_len))
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
