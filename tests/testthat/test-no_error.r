# This file test that external functions of the package
# are working
# (i.e., do not test behavior)

library(phytools)

test_that("fit_sinba works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(fit_sinba(tree, data_xz, new_model("IND")))
})

test_that("fit_sinba works with single birth", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(fit_sinba(tree, data_xz, new_model("IND"),
    single_birth = TRUE
  ))
})

test_that("fit_sinba works with hidden model", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)
  hm <- new_hidden_model(
    list(
      list(
        trait = "x",
        state = 1,
        hidden = c("a", "b", "c")
      )
    )
  )

  expect_no_error(fit_sinba(tree, data_xz, hm))
})

test_that("pagel works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(fit_pagel(tree, data_xz, new_model("IND")))
})

test_that("pagel works with hidden model", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  hm <- new_hidden_model(
    list(
      list(
        trait = "x",
        state = 1,
        hidden = c("a", "b", "c")
      )
    )
  )

  expect_no_error(fit_pagel(tree, data_xz, hm))
})

test_that("fit_fixed_matrix works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  fs <- fit_sinba(tree, data_xz, new_model("IND"))

  expect_no_error(fit_fixed_matrix(
    tree, data_xz, fs$Q,
    model = new_model("IND")
  ))
})

test_that("fit_fixed_matrix works with hidden model", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)
  hm <- new_hidden_model(
    list(
      list(
        trait = "x",
        state = 1,
        hidden = c("a", "b", "c")
      )
    )
  )

  fs <- fit_sinba(tree, data_xz, hm)

  expect_no_error(fit_fixed_matrix(
    tree, data_xz, fs$Q,
    model = hm
  ))
})

test_that("fixed_sinba works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  fp <- fit_pagel(tree, data_xz, new_model("IND"))
  root_id <- length(tree$tip.label) + 1
  births <- list()
  for (i in 1:2) {
    b <- list(
      node = root_id,
      age = 0
    )
    births[[i]] <- b
  }

  expect_no_error(fixed_sinba(
    tree, data_xz, fp$Q, births
  ))
})

test_that("fixed_sinba works with hidden models", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)
  hm <- new_hidden_model(
    list(
      list(
        trait = "x",
        state = 1,
        hidden = c("a", "b", "c")
      )
    )
  )

  fp <- fit_pagel(tree, data_xz, model = hm)
  root_id <- length(tree$tip.label) + 1
  births <- list()
  for (i in 1:2) {
    b <- list(
      node = root_id,
      age = 0
    )
    births[[i]] <- b
  }

  expect_no_error(fixed_sinba(
    tree, data_xz, fp$Q, births,
    model = hm
  ))
})

test_that("fit_sinba_single works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_z <- data.frame(LETTERS, z)

  expect_no_error(fit_sinba_single(tree, data_z, model = new_model(traits = 1)))
})

test_that("fit_sinba_single works with hidden model", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_z <- data.frame(LETTERS, z)
  hm <- new_hidden_model(
    list(
      list(
        trait = "x",
        state = 1,
        hidden = c("a", "b", "c")
      )
    ),
    traits = 1
  )

  expect_no_error(fit_sinba_single(tree, data_z, model = hm))
})

test_that("fixed_sinba_single works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_z <- data.frame(LETTERS, z)

  fz <- fit_sinba_single(tree, data_z)

  expect_no_error(fixed_sinba_single(
    tree, data_z, fz$Q, fz$birth
  ))
})

test_that("fixed_sinba works with hidden models", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_z <- data.frame(LETTERS, z)
  hm <- new_hidden_model(
    list(
      list(
        trait = "x",
        state = 1,
        hidden = c("a", "b", "c")
      )
    ),
    traits = 1
  )

  fz <- fit_sinba_single(tree, data_z, model = hm)
  expect_no_error(fixed_sinba_single(
    tree, data_z, fz$Q, fz$birth,
    model = hm
  ))
})

test_that("map_sinba works", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  fz <- fit_sinba(tree, data_xz, new_model("IND"))

  expect_no_error(map_sinba(fz, n = 2))
})

test_that("map_sinba from raw data works", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(map_sinba(
    tree = tree, data = data_xz,
    model = new_model("IND"), n = 2
  ))
})
