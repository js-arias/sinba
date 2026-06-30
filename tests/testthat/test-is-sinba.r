library(phytools)

test_that("is_full_sinba", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  y <- c(0, 0, rep(1, 5), rep(0, 19))
  z <- c(rep(0, 20), rep(1, 6))
  w <- c(1, rep(1, 24), 1)

  data_xz <- data.frame(LETTERS, x, z)

  # single trait
  expect_equal(
    is_full_sinba(tree, data.frame(LETTERS, x), new_model(traits = 1)),
    TRUE
  )
  expect_equal(
    is_full_sinba(tree, data.frame(LETTERS, w), new_model(traits = 1)),
    FALSE
  )

  # two traits
  expect_equal(
    is_full_sinba(tree, data.frame(LETTERS, x, x)),
    TRUE
  )
  expect_equal(
    is_full_sinba(tree, data.frame(LETTERS, x, y)),
    TRUE
  )
  expect_equal(
    is_full_sinba(tree, data.frame(LETTERS, x, w)),
    FALSE
  )
  expect_equal(
    is_full_sinba(tree, data.frame(LETTERS, x, z)),
    FALSE
  )
})
