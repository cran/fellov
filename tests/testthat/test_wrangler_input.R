context("Wrangeler Input")

P <- matrix(c(2,0,1,0,1.5,-0.2,1,-0.2,1), ncol = 3)
e1 <- list("c" = c(1,2,3), "P" = P, "r" = 2)
e2 <- wrangle_ellipse(e1, out_params = c("S"))
e3 <- wrangle_ellipse(e1, out_params = c("L"))
e4 <- wrangle_ellipse(e1, out_params = c("e", "U"))
e5 <- wrangle_ellipse(e1, out_params = c("D", "U"))
e6 <- wrangle_ellipse(e3, out_params = c("S"))
e7 <- wrangle_ellipse(e4, out_params = c("S")) #
e8 <- wrangle_ellipse(e5, out_params = c("S")) #
e9 <- wrangle_ellipse(e2, out_params = c("L"))
e10 <- wrangle_ellipse(e4, out_params = c("L"))
e11 <- wrangle_ellipse(e5, out_params = c("L"))
e12 <- wrangle_ellipse(e2, out_params = c("e", "U"))
e13 <- wrangle_ellipse(e3, out_params = c("e", "U"))
e14 <- wrangle_ellipse(e5, out_params = c("e", "U"))
e15 <- wrangle_ellipse(e2, out_params = c("D", "U"))
e16 <- wrangle_ellipse(e3, out_params = c("D", "U"))
e17 <- wrangle_ellipse(e4, out_params = c("D", "U"))

tol <- 1e-10

test_that("Wrangles correctly", {
  expect_true(all(wrangle_ellipse(e1, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e2, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e3, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e4, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e5, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e6, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e7, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e8, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e9, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e10, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e11, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e12, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e13, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e14, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e15, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e16, "P")[[1]] - P < tol))
  expect_true(all(wrangle_ellipse(e17, "P")[[1]] - P < tol))
})


test_that("ellipoids", {
  expect_error(wrangle_ellipse(ell = NULL, out_params = NULL),
               "'ell' must be a named list or a list of named lists.")
  expect_error(wrangle_ellipse(ell = list("c" = c(1),
                                          "P" = matrix(c(1,2,2,4), ncol = 2)),
                                 out_params = c("c", "P")),
               "The different parameters imply different dimensions.")
})

test_that("out parameters", {
  expect_error(wrangle_ellipse(ell = list(), out_params = c("H")),
               "One or more characters from 'out_prams' are not recognised.")
  expect_error(wrangle_ellipse(ell = list(), out_params = c("H", "B")),
               "One or more characters from 'out_prams' are not recognised.")
})
