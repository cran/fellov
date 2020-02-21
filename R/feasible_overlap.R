#' Find Smallest Feasible Ellipse Overlap
#'
#' \code{feasible_overlap} will find the smallest radius such that the ellipses
#' have a non-empty intersection.
#'
#' Given a list of ellipses \code{ell} the function \code{feasible_overlap} will
#' find the smallest radius such that the ellipses from \code{ell} overlap. This
#' is done by solving the following quadratically constrained problem
#'
#' \tabular{cclr}{
#' \ifelse{html}{\out{&nbsp;&nbsp;&nbsp;&nbsp;&nbsp}}{\eqn{\qquad}}\tab
#' \ifelse{html}{\out{&nbsp min<sub>(x,s)</sub> &nbsp}}{\eqn{\min_{(x,s)}\quad}}
#'  \tab s \tab \cr
#' \tab s.t. \tab
#' \ifelse{html}{\out{(x - c<sub>i</sub>)<sup>T</sup> P<sub>i</sub> (x - c<sub>i</sub>) - r<sub>i</sub> &#8804 s}}{\eqn{(x - c_i)^T P_i (x - c_i) - r_i \leq s}}
#' \tab \ifelse{html}{\out{&nbsp;&nbsp}}{\eqn{\quad}} for all i = 1, ..., d\cr
#' }
#' To solve this convex problem the logarithmic barrier method is used.
#'
#' Note that it is not necessary to specify ellipse radii in \code{ell}.
#'
#'
#' @param ell a list of at least two (non degenerate) ellipses; see
#'   \code{\link{wrangle_ellipse}}.
#' @param ... additional arguments to be passed to internal functions.
#'
#' @return \code{feasible_overlap} returns an object of \code{\link[base]{class}}
#' "\code{feasible_overlap}". This object is a list with the following entries:
#'   \item{radii}{the smallest ellipse radii resulting in a non-empty
#'     intersection.}
#'   \item{x}{The limiting common point.}
#'   \item{distance}{The ellipse specific distances.}
#'   \item{call}{The matched call.}
#'
#' @seealso \code{\link{wrangle_ellipse}} for detailed on ellipse
#' parameterization.
#'
#' @examples
#' ## two dimensional ellipses
#' e1 <- list("c" = c(1,1), "P" = matrix(c(2,0,0,0.5), ncol = 2))
#' e2 <- list("c" = c(0,0), "S" = matrix(c(1, 0.2, 0.2, 2), ncol = 2), "r" = 1)
#' # note: it is not necessary to specify an ellipse radius "r"
#'
#' feasible_overlap(list(e1, e2))
#'
#'
#' ## regression example
#' # generate data
#' n <- 100
#' E <- rbinom(n, 1, 0.5)
#' X <- rnorm(n, E * 3, 1)
#' Y <- rnorm(n, 2 + 1.5 * X, 1)
#'
#' # create confidence region ellipses
#' m0 <- lm(Y ~ X, data = data.frame(Y, X), subset = (E == 0))
#' m1 <- lm(Y ~ X, data = data.frame(Y, X), subset = (E == 1))
#' ConfRegion0 <- list(c = coefficients(m0), S = vcov(m0))
#' ConfRegion1 <- list(c = coefficients(m0), S = vcov(m0))
#'
#' # find smallest radius
#' res <- feasible_overlap(list(ConfRegion0, ConfRegion1))
#' # this radius now corresponds to the chisq quantile at which
#' # the confidence regions intersect non-emptily.
#' # In other words the (1 - alpha)-confidence intervals intersect for alpha:
#' alpha <- pchisq(res$radii, 2)
#'
#' @export

feasible_overlap <- function(ell, ...) {
  if (!is.list(ell)) {stop("'ell' must be a list of ellipses")}
  if (length(ell) < 2) {stop("'ell' must be a list of length at least 2.")}
  if (! all(sapply(ell, function(e) {is.list(e)}))) {
    stop("'ell' must be a list of ellipses.")
  }
  new_ell <- wrangle_ellipse(ell, out_params = c("n", "c", "P", "q"))
  new_ell <- lapply(new_ell, function(e) {c(e, "r" = 0)})
  nn <- sapply(new_ell, function(e) {e$n})
  if (! all(nn[1] == nn)) {
    stop("The ellipses do not have the same dimensions.")
  }
  res <- find_overlap(new_ell, ...)
  dist <- data.frame("distance" = distance_to_center(res$x, new_ell),
                     row.names = names(ell))


  res <- list(radii = res$s, x = res$x, distance = dist, call = match.call())
  return(structure(res, class = "feasible_overlap"))
}

#' @export
print.feasible_overlap <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Limiting Intersection Point:\n  ", x$x, "\n\n")
  cat("Ellipse Specific Distance:\n")
  pp <- seq_len(nrow(x$distance))
  if (is.null(row.names(x$distance))) {
    row.names(x$distance) <- paste0("ellipse ", pp, ": ")
  }
  if (all(row.names(x$distance) == pp)) {
    row.names(x$distance) <- paste0("ellipse ", pp, ": ")
  }

  print(x$distance)
  cat("\n")
  invisible(x)
}

