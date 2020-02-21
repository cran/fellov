#' Find Feasible Point in Ellipse Overlap
#'
#' \code{feasible_point} will find a point in the interior of the intersection
#' of two or more fully specified ellipses. If the intersections is empty
#' \code{NA} is returned.
#'
#' \code{feasible_point} will find a point in the interior of the intersection
#' of two or more fully specified ellipses \code{ell}. If the intersections is
#' empty \code{NA} is returned.
#'
#' @param ell a list of at least two (non degenerate) ellipses; see
#'   \code{\link{wrangle_ellipse}}.
#' @param ... additional arguments to be passed to internal functions.
#'
#' @return \code{feasible_point} returns an object of \code{\link[base]{class}}
#'   "\code{feasible_point}" with the following entries
#'   \item{x }{ An interior point.}
#'   \item{distance }{ A data.frame with the ellipse specific distances.}
#'   \item{optim }{ The final internal optimization value.}
#'   \item{call }{ The matched call.}
#'
#' @seealso \code{\link{wrangle_ellipse}} for detailed on ellipse
#' parameterization.
#'
#' @examples
#' # two dimensional ellipses
#' e1 <- list("c" = c(1,2), "P" = matrix(c(2,0,0,1), ncol = 2), "r" = 3)
#' e2 <- list("c" = c(0,0), "S" = matrix(c(1, 0.2, 0.2, 2), ncol = 2), "r" = 1)
#'
#' # find point in intersection
#' feasible_point(list(e1, e2))
#'
#' # make new ellipse
#' e3 <- list("c" = c(2,2), "P" = matrix(c(1,0,0,1), ncol = 2), "r" = 0.5)
#'
#' # now there is no overlap
#' feasible_point(list(e1, e2, e3))
#'
#' @export

feasible_point <- function(ell, ...) {
  if (!is.list(ell)) {stop("'ell' must be a list of ellipses")}
  if (length(ell) < 2) {stop("'ell' must be a list of length at least 2.")}
  if (! all(sapply(ell, function(e) {is.list(e)}))) {
    stop("'ell' must be a list of ellipses.")
  }
  new_ell <- wrangle_ellipse(ell, out_params = c("n", "c", "P", "q", "r"))
  nn <- sapply(new_ell, function(e) {e$n})
  if (! all(nn[1] == nn)) {
    stop("The ellipses do not have the same dimensions.")
  }
  res <- find_overlap(new_ell, ...)
  dist <- distance_to_center(res$x, new_ell)
  radii <- sapply(new_ell, function(e) {e$r})
  int <- dist <= radii
  out_dist <- data.frame("distance" = dist, "radius" = radii, "interior" = int,
                         row.names = names(ell))
  if (res$s <= 0) {xx <- res$x} else {xx <- NA}
  res <- list(x = xx, distance = out_dist, optim = res$s, call = match.call())
  return(structure(res, class = "feasible_point"))
}

#' @export
print.feasible_point <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Interior Point:\n  ")
  if (x$optim >= 0) {
    cat("The ellipses do not intersect !")
  } else {
    cat(x$x, "\n\nEllipse Specific Distances:\n")
    if (is.null(row.names(x$distance))) {
      row.names(x$distance) <- paste0("ellipse ", seq_len(nrow(x$distance)), ": ")
    }
    print(x$distance)
  }
  cat("\n")
  invisible(x)
}





#' Determine If A Point Is In Ellipse Overlap
#'
#' \code{is_feasible_point} will determine if a given point is in the interior
#' of the intersection of one or more fully specified ellipes.
#'
#' Given a point \code{is_feasible_point} will check if this point is in the
#' intersection of the list of ellipses \code{ell}. Note that this function will
#' not check if the intersection is non-empty.
#'
#' @param point a numeric of length equal to the dimensions of the ellipses in
#'   \code{ell}.
#' @param ell a list of at least one ellipse; see \code{\link{wrangle_ellipse}}.
#'
#' @return \code{is_feasible_point} returns an object of
#'   \code{\link[base]{class}} "\code{is_feasible_point}". This object is a list
#'   containing the following components:
#'   \item{point }{ the inputted \code{point}.}
#'   \item{fasible }{ logical; is \code{TRUE} when the point \code{x} is in the
#'     interior of all ellipses.}
#'   \item{distance }{ a data.frame with the distance from \code{x} to the center
#'     of each ellipse, the radius of each ellipse and a logical indicator,
#'     which is \code{TRUE} when \code{x} is an element in the ellipse.}
#'   \item{call }{ the match call.}
#'
#' @seealso \code{\link{wrangle_ellipse}} for detailed on ellipse
#' parameterization.
#'
#' @examples
#' e1 <- list("c" = c(1,1), "P" = matrix(c(3,1,1,2), ncol = 2), "r" = 2)
#' e2 <- list("c" = c(0,2), "S" = matrix(c(4,1,1,1), ncol = 2), "r" = 3)
#'
#' is_feasible_point(c(1.1,0.9), e1)
#' is_feasible_point(c(1,0), list(e1, e2))
#'
#' @export

is_feasible_point <- function(point, ell) {
  if (!is.numeric(point)) {
    stop("'point' must be a numeric vector or matrix.")
  }
  if (!is.list(ell)) {
    stop("'ell' must be a named list or a list of named lists.")
  } else if (any(sapply(ell, function(e) {! is.list(e)}))) {
    ne <- 1
    ell <- list(ell)
  } else {
    ne <- 2
  }
  new_ell <- wrangle_ellipse(ell, out_params = c("n", "c", "P", "r"))
  nn <- sapply(new_ell, function(e) {e$n})
  if (! all(length(point) == nn)) {
    stop("One or more ellipse does not have the same dimensions as 'point'.")
  }
  dist <- distance_to_center(point, new_ell)
  radii <- sapply(new_ell, function(e) {e$r})
  int <- dist <= radii
  if (is.null(names(new_ell))) {
    nam <- paste0("ellipse ", seq_along(new_ell))
  } else {
    nam <- names(new_ell)
    if (any(nam == "")) {
      which <- which(nam == "")
      nam[which] <- paste0("ellipse ", which)
    }
  }
  out_dist <- data.frame("distance" = dist, "radius" = radii, "interior" = int,
                         row.names = nam)

  res <- list(point = point, feasible = all(int), distance = out_dist,
              call = match.call())
  return(structure(res, class = "is_feasible_point"))
}

#' @export
print.is_feasible_point <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  txt <- paste(deparse(x$call$point))
  nn <- nchar(txt)
  if (nn >= 12) {
    txt <- paste0(substring(txt, 1, 6),"...", substring(txt, nn - 3, nn))
  }
  cat("\nThe point", txt, "is in the ellipse intersection:\n  ", x$feasible)
  cat("\n\nDistances by Ellipse:\n")
  print(x$distance)
  cat("\n")
  invisible(x)
}


distance_to_center <- function(x, ell) {
  sapply(ell, function(e) {
    diff <- matrix(x - e$c, ncol = 1)
    t(diff) %*% e$P %*% diff
  })
}
