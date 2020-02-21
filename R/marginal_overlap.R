#' Feasibility of all Marginal Ellipse Overlaps
#'
#' Determin if the projections of ellipses onto each margin overlap.
#'
#' The ellipses are projected onto the specified \code{margins}. For each margin
#' the intersection of the projected ellipses is found. The  \code{Lower} and
#' \code{Upper} endpoints of the intersection interval is reported. If the
#' intersection along a margin is empty then \code{Lower} and \code{Upper}
#' is reported as \code{NA}.
#'
#' Note that if the ellipses overlap when projected onto each margin this does
#' not imply that the ellipses themselves intersect non-emptily. The example
#' below is constructed to illustrate this.
#'
#' @param ell a list of at least two (non degenerate) ellipses; see
#'   \code{\link{wrangle_ellipse}}.
#' @param margins either "all" or a vector indicating the margins to project the
#'   ellipses onto and take intersections.
#'
#' @return \code{marginal_overlap} returns an object of \code{\link[base]{class}}
#' "marginal_overlap" which contains a data.frame where the coloumns descibe the
#' following
#'   \item{Margin}{Intputted \code{margins}.}
#'     \item{Overlap}{Whether the ellipses overlap when projected onto
#'       corresponding margin.}
#'     \item{Lower}{Lower endpoint of intersection interval. \code{NA} if the
#'       intersection is empty.}
#'     \item{Upper}{Upper endpoint of intersection interval. \code{NA} if the
#'       intersection is empty.}
#'
#' @seealso \code{\link{wrangle_ellipse}} for detailed on ellipse
#' parameterization.
#'
#' @examples
#' ## two dimensional ellipses
#' e1 <- list(c = c(0.1, 0), P = matrix(c(3, 0, 0, 1), ncol = 2), r = 1)
#' e2 <- list(c = c(1, 1), P = matrix(c(3, 1.2, 1.2, 1), ncol = 2), r = 0.8)
#' e3 <- list(c = c(2, 1.5), P = matrix(c(1, 0.6, 0.6, 1), ncol = 2), r = 0.4)
#' # Note: These three ellipses have been chosen so (some of) the marginal
#' #       projections intersect, but the actual ellipses do not intersect.
#'
#' # Ellipses e1 and e2 overlap when projected onto margin 1 and 2 respectivly.
#' marginal_overlap(list(e1, e2))
#'
#' # Adding ellipse e3:
#' # Then there is no overlap when projecting onto margin 1
#' marginal_overlap(list(e1, e2, e3), margins = c(1))
#'
#'
#' ## regression example
#' n <- 100
#' E <- rbinom(n, 1, 0.5)
#' X <- rnorm(n, E * 3, 1)
#' Y <- rnorm(n, 2 + X, 1)
#' lm_E0 <- lm(Y ~ X, data = data.frame(Y, X), subset = (E == 0))
#' lm_E1 <- lm(Y ~ X, data = data.frame(Y, X), subset = (E == 1))
#'
#' # create 95% confidence ellipses and check marginal overlap
#' q <- qchisq(0.95, 2) # df = 2, as there are two covariates (1, X)
#' ell0 <- list(c = coefficients(lm_E0), S = vcov(lm_E0), r = q)
#' ell1 <- list(c = coefficients(lm_E1), S = vcov(lm_E1), r = q)
#' marginal_overlap(list(ell0, ell1))
#' @export
marginal_overlap <- function(ell, margins = "all") {
  if (!is.list(ell)) {
    stop("'ell' must be a list of ellipses")
  }
  if (length(ell) < 2) {
    stop("'ell' must be a list of length at least 2.")
  }

  new_ell <- lapply(ell, function(e) {
    wrangle_ellipse(e, out_params = c("n", "c", "S", "r"))
  })

  dim <- new_ell[[1]]$n
  if (identical(margins,"all")) {
    n <- dim
    mar <- seq_len(n)
  } else if (is.numeric(margins)) {
    n <- length(margins)
    mar <- margins
    if (n > dim) {
      stop("Too many 'margins' specified.")
    }
    if (any(mar > dim)) {
      stop("'margins' do not correspond to ellipse dimension.")
    }
    if (length(unique(mar)) < length(mar)) {
      stop("'margins' contains duplicates.")
    }
  } else {
    stop("'margins' must either be 'all' or a vector of integers indicating the margins.")
  }

  end <- sapply(new_ell, function(e) {
    q <- sqrt(e$r * diag(e$S)[mar])
    L <- e$c[mar] - q
    U <- e$c[mar] + q
    return(c(L, U))
  })

  if (n == 1) {
    L_max <- max(end[1,])
    U_min <- min(end[2,])
  } else {
    L_max <- apply(end[seq_len(n), ], 1, max)
    U_min <- apply(end[seq_len(n) + n, ], 1, min)
  }
  overlap <- L_max <= U_min

  intersection <- data.frame("Margin" = mar, "Overlap" = overlap,
                             "Lower" = L_max, "Upper" = U_min)
  intersection[!overlap, c(3,4)] <- c(NA, NA)

  return(structure(list("intersection" = intersection, "call" = match.call()),
                   class = "marginal_overlap"))
}

#' @export
print.marginal_overlap <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("The Ellipses Overlap Along All Specified Margins:\n  ")
  if (all(x$intersection[,2]) == 1) {cat("TRUE\n\n")} else {cat("FALSE\n\n")}

  cat("Margin Specific Overlap:\n")
  interval <- apply(x$intersection, 1, function(d) {
    if (is.na(d["Lower"])) {
      "  Empty Intersection"
    } else {
      paste0("  [", round(d["Lower"], digits = 3), "; ",
             round(d["Upper"], digits = 3) ,"]")
    }
  })
  row_names <- if (all(row.names(x$intersection) == paste0(seq_len(nrow(x$intersection))))) {
    paste0("Margin ", x$intersection[,1])
  } else {
    row.names(x$intersection)
  }
  dat <- data.frame("Overlap" = as.logical(x$intersection[,2]),
                    interval,
                    row.names = row_names,
                    fix.empty.names = FALSE)
  print(dat)

  cat("\n")
  invisible(x)
}
