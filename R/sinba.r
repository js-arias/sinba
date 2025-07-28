# node_age returns the elapsed time
# since the root time.
node_age <- function(t, n) {
  s <- 0
  while (TRUE) {
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    s <- s + t$edge.length[e]
    n <- t$edge[e, 1]
  }
  return(s)
}

# path_to_node returns the path from the root
# to a given node.
path_to_node <- function(t, n) {
  x <- c()
  while (TRUE) {
    x <- c(n, x)
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    n <- t$edge[e, 1]
  }
  return(x)
}

# is_parent_in_tree returns true if node p is parent of n.
is_parent_in_tree <- function(t, p, n) {
  if (p == n) {
    # a node cannot be parent of itself
    return(FALSE)
  }
  while (TRUE) {
    e <- which(t$edge[, 2] == n)
    if (length(e) == 0) {
      break
    }
    a <- t$edge[e, 1]
    if (a == p) {
      return(TRUE)
    }
    n <- a
  }
  return(FALSE)
}

# node_length returns a list with the total distance
# from the node to the root.
node_lengths <- function(t) {
  nl <- c()
  for (i in seq_len(length(t$tip.label) + t$Nnode)) {
    l <- length_to_root(t, i)
    nl <- c(nl, l)
  }
  return(nl)
}

# get_node_by_age returns a node in a path
# using the age a of the node
get_node_by_age <- function(p, ages, a) {
  for (i in p) {
    if (ages[i] >= a) {
      return(i)
    }
  }
  return(-1)
}

# node_size returns the number of terminals
# for a given node.
node_size <- function(t, n) {
  if (n <= length(t$tip.label)) {
    return(1)
  }
  sum <- 0
  for (x in seq_len(length(t$tip.label))) {
    if (is_parent_in_tree(t, n, x)) {
      sum <- sum + 1
    }
  }
  return(sum)
}
