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
