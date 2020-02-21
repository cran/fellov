# QCQP Solver for Feasibility Problems
#
# \code{find_overlap} is used to find the smallest radius of two or more
# ellipses that results in the existance of a feasible point, i.e. the ellipse
# intersect.
#
# The function \code{find_overlap} finds the smallest possible radius resulting
# in a non-empty intersection of the ellipses \code{ell}. This is done by
# solving the following quadratically constrained problem
#
# \tabular{cclr}{
# \ifelse{html}{\out{&nbsp;&nbsp;&nbsp;&nbsp;&nbsp}}{\eqn{\qquad}}\tab
# \ifelse{html}{\out{&nbsp min<sub>(x,s)</sub> &nbsp}}{\eqn{\min_{(x,s)}\quad}}
#  \tab s \tab \cr
# \tab s.t. \tab
# \ifelse{html}{\out{(x - c<sub>i</sub>)<sup>T</sup> P<sub>i</sub> (x - c<sub>i</sub>) - r<sub>i</sub> &#8804 s}}{\eqn{(x - c_i)^T P_i (x - c_i) - r_i \leq s}}
# \tab \ifelse{html}{\out{&nbsp;&nbsp}}{\eqn{\quad}} for all i = 1, ..., d\cr
# }
# To solve this convex problem the logarithmic barrier method is used.
#
#
# @param ell a list of length at least 2 describing the ellipses.
# @param initial a numeric vectors or matrix with initial values for the
#   parameters to be optimized over.
# @param centering_steps number of centering steps, must be atleast 1.
# @param tol_newton tolerance for the newton method in each centering step.
# @param mu tuning parameter for centering method.
# @param alpha tuning parameter for backtracking lineseach in newton  steps.
# @param beta tuning parameter for backtracking lineseach in newton  steps.
# @param tol_step meachnism for jumping to the next centoring step
#
# @seealso \code{\link{wrangle_ellipse}} for detailed on ellipse
# parameterization. The functions \code{\link{feasible_overlap}} and
# \code{\link{feasible_point}} both use \code{find_overlap} internally.
#
# @return
# The function \code{find_overlap} returns an object of \code{\link[base]{class}}
# "\code{find_overlap}" with the following components:
# \item{x }{ the global minimum; a point in the smallest possible intersection.}
# \item{s }{ the optimal value at the global minimum.}
# \item{call }{ the matched call.}
#
#
# @keywords internal
#
# @export
find_overlap <- function(ell, initial, centering_steps = 3L, tol_newton = 0.1, mu = 20, alpha = 0.3, beta = 0.8, tol_step = 1e-5) {
  if (! is.numeric(centering_steps)) {stop("'centering_steps' must be an integer larget than 1.")}
  centering_steps <- as.integer(centering_steps)
  if (centering_steps < 1L) {stop("'centering_steps' must be an integer larger than 1.")}
  if (!is.numeric(tol_newton)) {stop("'tol_newton' must be a small positive number.")}
  if (tol_newton <= 1e-10) {stop("'tol_newton' must be larger then 1e-10.")} # IS THIS TOO RESTRICTIVE ???
  if (! is.numeric(mu)) {stop("'mu' must be a number strictly larger than 1.")}
  if (mu <= 1) {stop("'mu' must be a number strictly larger than 1.")}
  if (! is.numeric(alpha)) {stop("'alpha' must be a number between 0 and 0.5.")}
  if (alpha <= 0 | alpha >= 0.5) {stop("'alpha' must be a number between 0 and 0.5.")}
  if (! is.numeric(beta)) {stop("'beta' must be a number between 0 and 1.")}
  if (beta <= 0 | beta >= 1) {stop("'beta' must be a number between 0 and 1.")}
  if (!is.numeric(tol_step)) {stop("'tol_step' must be a small positive number.")}
  if (tol_step <= 1e-10) {stop("'tol_step' must be larger then 1e-10.")} # IS THIS TOO RESTRICTIVE ???

  if (missing(initial)) {
    tmp <- sapply(ell, function(e) {e$c})
    if (is.null(dim(tmp))) { x <- mean(tmp) } else { x <- rowMeans(tmp) }
  } else if (length(initial) == length(ell[[1]]$c)) {
    if (is.numeric(initial)) {
      x <- initial
    } else {
      stop("'initial' must be a numeric vector of length equal to ellipse dimension.")
    }
  } else {
    stop("'initial' must be a numeric vector of length equal to ellipse dimension.")
  }

  # -- Define relevamt functions
  pp <- seq_along(x)
  distance_to_edge <- function(x) {
    sapply(ell, function(e) {
      diff <- x - e$c
      t(diff) %*% e$P %*% diff - e$r
    })
  }
  speed <- function(x) {
    lapply(ell, function(e) {
      2 * e$P %*% x + e$q
    })
  }
  acceleration <- function(x) {
    lapply(ell, function(e) {
      2 * e$P
    })
  }
  target <- function(s, tune, dist) {
    tune * s - sum(log(s - dist))
  }
  target_gradient <- function(s, tune, dist, spee) {
    top <- tune - sum(1 / (s - dist))
    v <- lapply(pp, function(i) {
      (1 / (s - dist[[i]])) * spee[[i]]
    })
    bottom <- Reduce("+", v)
    return(matrix(c(top, bottom), ncol = 1))
  }
  target_hessian <- function(s, tune, dist, spee, acc) {
    dsds <- sum(1 / (s - dist) ^ 2)
    v <- lapply(pp, function(i) {
      (1 / (s - dist[[i]]) ^ 2) * spee[[i]]
    })
    dsdx <- - Reduce("+", v)
    m <- lapply(pp, function(i) {
      (1 / (s - dist[[i]]) ^ 2) * spee[[i]] %*% t(spee[[i]]) +
        (1 / (s - dist[[i]])) * acc[[i]]
    })
    dxdx <- Reduce("+", m)
    return(cbind(rbind(dsds, dsdx), rbind(t(dsdx), dxdx)))
  }
  linesearch_test <- function(min_s, s, x, dist, decrement, tune, delta, step) {
    s_delta <- s + step * delta[1, ]
    if (s_delta <= min_s) {return(TRUE)}
    dist_delta <- distance_to_edge(x + step * delta[-1, ])
    if (max(dist_delta) >= s_delta){return(TRUE)}
    if(target(s_delta, tune, dist_delta) >
       (target(s, tune, dist) + alpha * step * decrement)) {
      return(TRUE)
    }
    return(FALSE)
  }



  #  -- Centering
  min_s <- sum(sapply(ell, function(e) {-e$r}))
  dist <- distance_to_edge(x)
  s <- max(dist) + 1e-5
  tune <- (sum(1 / (s - dist)) * s + 1) / s

  for (i in seq_len(centering_steps)) {

    # -- centering :  Newtons method
    lambda <- tol_newton * 5
    ss <- TRUE # allowes jumps to next centering if s becomes smaller than 0
    while ( (lambda / 2 > tol_newton) & ss) {
      # .. calculate direction and decrement
      dist <- distance_to_edge(x)
      spee <- speed(x)
      acc <- acceleration(x)
      gradient <- target_gradient(s, tune, dist, spee)
      hessian <- target_hessian(s, tune, dist, spee, acc)
      if (rcond(hessian) < 2e-15) {
        add <- diag(diag(hessian) * 1e-3)
        inv_hessian <- solve(hessian + add)
      } else {
        inv_hessian <- solve(hessian)
      }
      delta <- - inv_hessian %*% gradient
      lambda <- - t(gradient) %*% delta

      # .. linesearch
      step <- 1
      while (linesearch_test(min_s, s, x, dist, lambda, tune, delta, step)) {
        if (step < tol_step) {
          step <- 0
          ss <- FALSE
        } else {
          step <- beta * step
        }
      }

      # .. update (s,x)
      s <- s + step * delta[1, ]
      x <- x + step * delta[-1, ]
      if (s < 0) {
        ss <- FALSE
        }
    }
    tune <- tune * mu
  }

  return(structure(list("x" = x, "s" = s, "call" = match.call()),
                   class = "find_overlap"))
}
