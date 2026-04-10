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

test_that("fit_sinba works with a simultaneous birth", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(fit_simultaneous(tree, data_xz, new_model("IND")))
})

test_that("fit_sinba works with hidden model", {
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

  expect_no_error(fit_sinba(tree, data_xz, model_as(hm, "IND")))
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

  expect_no_error(fit_pagel(tree, data_xz, model_as(hm, "IND")))
})

test_that("fit_mixed works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(fit_mixed(
    tree, data_xz,
    trait = "x"
  ))
})

test_that("fit_fixed_births works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  births <- list(
    list(
      node = 31,
      age = 0
    ),
    list(
      node = 28,
      age = 0
    )
  )

  expect_no_error(fit_fixed_births(
    tree, data_xz, births
  ))
})

test_that("fit_fixed_matrix works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  q <- matrix(c(
    0,    1e-6, 2, 0,
    3,    0,    0, 2,
    1e-6, 0,    0, 1e-6,
    0,    1e-6, 3, 0
  ), byrow = TRUE, nrow = 4)

  expect_no_error(fit_fixed_matrix(
    tree, data_xz, q,
    model = new_model("IND")
  ))
})

test_that("fit_fixed_matrix works with hidden model", {
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
  hq <- matrix(c(
    0,    3,    3,    0,   3,    0,   3,   0,
    1e-6, 0,    0,    3,   0,    3,   0,   3,
    0.01, 0,    0,    3,   0.5,  0,   0.5, 0,
    0,    0.01, 1e-6, 0,   0,    0.5, 0,   0.5,
    0.01, 0,    0.5,  0,   0,    3,   0.5, 0,
    0,    0.01, 0,    0.5, 1e-6, 0,   0,   0.5,
    0.01, 0,    0.5,  0,   0.5,  0,   0,   3,
    0,    0.01, 0,    0.5, 0,    0.5, 0,   0
  ), byrow = TRUE, nrow = 8)

  expect_no_error(fit_fixed_matrix(
    tree, data_xz, hq,
    model = model_as(hm, "IND")
  ))
})

test_that("fixed_sinba works", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  births <- list(
    list(
      node = 31,
      age = 0
    ),
    list(
      node = 28,
      age = 0
    )
  )

  q <- matrix(c(
    0,    1e-6, 2, 0,
    3,    0,    0, 2,
    1e-6, 0,    0, 1e-6,
    0,    1e-6, 3, 0
  ), byrow = TRUE, nrow = 4)

  expect_no_error(fixed_sinba(
    tree, data_xz, q, births
  ))
})

test_that("fixed_sinba works with hidden models", {
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

  hq <- matrix(c(
    0,    3,    3,    0,   3,    0,   3,   0,
    1e-6, 0,    0,    3,   0,    3,   0,   3,
    0.01, 0,    0,    3,   0.5,  0,   0.5, 0,
    0,    0.01, 1e-6, 0,   0,    0.5, 0,   0.5,
    0.01, 0,    0.5,  0,   0,    3,   0.5, 0,
    0,    0.01, 0,    0.5, 1e-6, 0,   0,   0.5,
    0.01, 0,    0.5,  0,   0.5,  0,   0,   3,
    0,    0.01, 0,    0.5, 0,    0.5, 0,   0
  ), byrow = TRUE, nrow = 8)

  births <- list(
    list(
      node = 31,
      age = 0
    ),
    list(
      node = 28,
      age = 0
    )
  )

  expect_no_error(fixed_sinba(
    tree, data_xz, hq, births,
    model = model_as(hm, "IND")
  ))
})

test_that("fit_sinba works on a single trait", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_z <- data.frame(LETTERS, z)

  expect_no_error(fit_sinba(tree, data_z, model = new_model(traits = 1)))
})

test_that("fit_sinba_single works with hidden model", {
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

  expect_no_error(fit_sinba(tree, data_z, model = hm))
})

test_that("fixed_sinba works on a single trait", {
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_z <- data.frame(LETTERS, z)

  birth <- list(
    list(
      node = 28,
      age = 0
    )
  )
  q <- matrix(c(
    0, 1,
    2, 0
  ), byrow = TRUE, nrow = 2)

  expect_no_error(fixed_sinba(
    tree, data_z, q, birth,
    model = new_model(traits = 1)
  ))
})

test_that("fixed_sinba works a single trait and hidden models", {
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
  hq <- matrix(c(
    0,    1.7, 0.07, 13,
    10,   0,   1.1,  1e-6,
    0.03, 5,   0,    15,
    25,   35,  10,   0
  ), byrow = TRUE, nrow = 4)

  birth <- list(
    list(
      node = 28,
      age = 0
    )
  )

  expect_no_error(fixed_sinba(
    tree, data_z, hq, birth,
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

test_that("map_pagel works", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  fz <- fit_pagel(tree, data_xz, new_model("IND"))

  expect_no_error(map_pagel(fz, n = 2))
})

test_that("map_pagel from raw data works", {
  skip_on_cran()
  set.seed(6)

  tree <- pbtree(n = 26, tip.label = LETTERS)
  x <- c(0, 0, rep(1, 12), rep(0, 12))
  z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))
  data_xz <- data.frame(LETTERS, x, z)

  expect_no_error(map_pagel(
    tree = tree, data = data_xz,
    model = new_model("IND"), n = 2
  ))
})
