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

  fs <- fit_sinba(tree, data_xz, hm)

  expect_no_error(fit_fixed_matrix(
    tree, data_xz, fs$Q,
    model = hm
  ))
})
