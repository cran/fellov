#' Feasibility of all Pairwise Ellipse Overlaps
#'
#' Determin if pairs of ellipses intersect non-emptily..
#'
#' The \code{pairwise_overlap} functions goes through all pairs of ellipses from
#' \code{ell} and checks if their intersection is non-empty.
#'
#' Note that if all pairs of ellipses intersect this does not mean that the
#' intersection of all the ellipses is non-empty. The example below is
#' constructed to illustrate this.
#'
#' @param ell a list of at least two (non degenerate) ellipses; see
#'   \code{\link{wrangle_ellipse}}.
#' @param ... additional arguments to be passed to the low level funcitons.
#'
#' @return
#' The \code{pairwise_overlap} function returns an object of
#' \code{\link[base]{class}} "\code{pairwise_overlap}" with the following
#' components:
#' \item{intersection }{ a data.frame where the two first columns specify the
#'   two ellipses intersected and the last coloumn indicate if they have a
#'   non-empty intersection.}
#' \item{call }{ the matched call.}
#'
#' @seealso \code{\link{wrangle_ellipse}} for detailed on ellipse
#' parameterization.
#'
#' @examples
#' ## three different two dimensional ellipses
#' e1 <- list(c = c(0, 0.7), P = matrix(c(0.2, 0, 0, 3), ncol = 2), r = 0.5)
#' e2 <- list(c = c(0, 1), P = matrix(c(3, -1.5, -1.5, 1), ncol = 2), r = 1)
#' e3 <- list(c = c(1.5, 1), P = matrix(c(3, 1.2, 1.2, 1), ncol = 2), r = 1.2)
#' # Note: These ellipses have been chosen so all pairs intersect,
#' #       but the intersection of all three is empty.
#'
#' # test pairwise overlaps
#' pairwise_overlap(list(e1, e2, e3))
#'
#' ## regression example
#' # generate data
#' n <- 100
#' E <- rbinom(n, 2, 0.5)
#' X <- rnorm(n, 3 * E, 1)
#' Y <- rnorm(n, 2 + 1.5 * E, 1)
#' m0 <- lm(Y ~ X, data = data.frame(Y,X), subset = (E == 0))
#' m1 <- lm(Y ~ X, data = data.frame(Y,X), subset = (E == 1))
#' m2 <- lm(Y ~ X, data = data.frame(Y,X), subset = (E == 2))
#'
#' # create 95% confidence ellipses and check pairwise intersection
#' q <- qchisq(0.95, 2) # df = 2, as there are two covariates (1, X)
#' E0 <- list(c = coefficients(m0), S = vcov(m0), r = q)
#' E1 <- list(c = coefficients(m1), S = vcov(m1), r = q)
#' E2 <- list(c = coefficients(m2), S = vcov(m2), r = q)
#' pairwise_overlap(list("model 0" = E0, "model 1" = E1, "model 2" = E2))
#'
#' @export
pairwise_overlap <- function(ell, ...) {
  if (!is.list(ell)) {
    stop("'ell' must be a list of ellipses")
  }
  if (length(ell) < 2) {
    stop("'ell' must be a list of length at least 2.")
  }

  new_ell <- lapply(ell, function(e) {
    wrangle_ellipse(e, out_params = c("n","c", "e", "U", "r"))
  })

  nam <- if (is.null(names(ell))) {
    paste0("ell ", seq_along(ell))
  }  else {
    names(ell)
  }
  if (any(nam == "")) {nam[nam == ""] <- paste0("ell ", which(nam == ""))}
  test_list <- utils::combn(seq_along(ell), 2)
  res <- sapply(seq_len(ncol(test_list)), function(i) {
    t <- test_list[,i]
    e1 <- new_ell[[t[1]]]
    e2 <- new_ell[[t[2]]]
    ell_nam <- c(nam[t[1]], nam[t[2]])

    diff <- e2[["c"]] - e1[["c"]]

    if (e1[["n"]] != e2[["n"]]) {
      stop("All ellipses must have equal dimension.")
    }
    if (e1[["n"]] == 1) {
      S2 <- e2[["U"]] %*% (1 / e2[["e"]]) %*% t(e2[["U"]])
      S1 <- e1[["U"]] %*% (1 / e1[["e"]]) %*% t(e1[["U"]])
      res <- diff ^ 2 <= S1 * e1[["r"]] + S2 * e2[["r"]]
      return(c(ell_nam, res))
    }
    if (t(diff)%*%e1[["U"]]%*%diag(e1[["e"]])%*%t(e1[["U"]])%*%diff<=e1[["r"]]) {

      return(c(ell_nam, TRUE))
    }
    if (t(diff)%*%e2[["U"]]%*%diag(e2[["e"]])%*%t(e2[["U"]])%*%diff<=e2[["r"]]) {
      return(c(ell_nam, TRUE))
    }

    q1 <- sqrt(e1[["r"]])
    D1 <- diag(sqrt(e1[["e"]]))
    S2 <- e2[["U"]] %*% diag(1 / e2[["e"]]) %*% t(e2[["U"]])

    # Transform
    c <- (1 / q1) * D1 %*% t(e1[["U"]]) %*% diff
    S <- D1 %*% t(e1[["U"]]) %*% S2 %*% e1[["U"]] %*% D1
    vals <- eigen(S)
    d <- sqrt(make_positive(vals$values, silent = TRUE))
    U <- vals$vectors

    y <- -t(U) %*% c
    y <- abs(y)

    # f goes monotone, quadratically to -1, so sure and fast convergence
    f <- function(t) sum(( d * y / (t + d^2))^2) - 1
    df <- function(t) - 2 * sum((y * d)^2 / (t + d^2)^3)

    t0 <- 0
    ft0 <- f(t0)

    while(ft0 > 1e-4){
      t0 <- t0 - f(t0) / df(t0)
      ft0 <- f(t0)
    }

    x0 <- y * d^2 / (d^2 + t0) # projection of y onto (c,C)
    dist <- sqrt(sum((y - x0)^2))
    return(c(ell_nam, dist < 1))
  })

  inter <- data.frame("ellipse1" = res[1,],
                      "ellipse2" = res[2,],
                      "intersection" = as.logical(res[3,]),
                      stringsAsFactors = FALSE)

  res <- structure(list("intersections" = inter, call = match.call()),
                   class = "pairwise_overlap")
  return(res)
}




# (very) adhoc way of dealing with negative eigenvalues
make_positive <- function (v, silent = TRUE){
  w <- which(v < 10^(-14))
  if (length(w) > 0 & !silent)
    warning("Some eigenvalues are below 10^(-14) and will therefore be regularized")
  for (ww in w) {
    if (v[ww] < 0) {
      v[ww] <- v[ww] - 2 * v[ww] + 10^(-14)
    }
    else {
      v[ww] <- v[ww] + 10^(-14)
    }
  }
  return(v)
}





#' @export
print.pairwise_overlap <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("All Pairs of Ellipses Intersect:\n  ")
  if (all(x$intersections[,3])) {cat("TRUE\n\n")} else {cat("FALSE\n\n")}

  cat("Pair Specific Intersection Info:\n")
  names(x$intersections) <- c(" Ellipse 1", " Ellipse 2", " Intersect")
  print(x$intersections)

  cat("\n")
  invisible(x)
}
