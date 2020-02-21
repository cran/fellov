context("Overlap functions")

e1_1 <- list("c" = c(1), "P" = matrix(5, ncol = 1), "r" = 0.8)
e1_2 <- list("c" = c(0), "P" = matrix(2, ncol = 1), "r" = 1.3)
e2_1 <- list("c" = c(1,1), "P" = matrix(c(1,0,0,1), ncol = 2), "r" = 0.9)
e2_2 <- list("c" = c(0,1), "P" = matrix(c(3,1,1,2), ncol = 2), "r" = 1.3)
e2_3 <- list("c" = c(1,0), "P" = matrix(c(2,0.5,0.5,1), ncol = 2), "r" = 1.2)

test_that("is_feasible_point", {
  expect_is(is_feasible_point(c(1), ell = e1_1), "is_feasible_point")
  expect_is(is_feasible_point(c(1,1), ell = e2_1), "is_feasible_point")
})

test_that("feasible_point", {
  expect_is(feasible_point(list(e1_1, e1_2)), "feasible_point")
  expect_is(feasible_point(list(e2_1, e2_2, e2_3)), "feasible_point")
})

test_that("feasible_overlap", {
  expect_is(feasible_overlap(list(e1_1, e1_2)), "feasible_overlap")
  expect_is(feasible_overlap(list(e2_1, e2_2, e2_3)), "feasible_overlap")
})

test_that("marginal_overlap", {
  expect_is(marginal_overlap(list(e1_1, e1_2)), "marginal_overlap")
  expect_is(marginal_overlap(list(e2_1, e2_2, e2_3)), "marginal_overlap")
})

test_that("pairwise_overlap", {
  expect_is(pairwise_overlap(list(e1_1, e1_2)), "pairwise_overlap")
  expect_is(pairwise_overlap(list(e2_1, e2_2, e2_3)), "pairwise_overlap")
})
