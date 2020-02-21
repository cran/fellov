#' Ellipse Wrangeler
#'
#' \code{wrangle_ellipse} is used to wrangle one or more ellipses from one
#' parametrization to another.
#'
#' Takes ellipse parameters and and calculates the wanted \code{out_params}. A
#' parameterization is a named list, where each named entry is a parameter. The
#' following parameters are accepted both input and output:
#'
#' \itemize{
#'   \item{\code{n} : }{ dimension of ellipse; an integer.}
#'   \item{\code{c} : }{ center of the ellipse; a vector.}
#'   \item{\code{P} : }{ precision matrix - inverse of \code{S}; a positive
#'     definit, symmetric matrix.}
#'   \item{\code{S} : }{ deviation matrix - inverse of \code{P}; a positive
#'     definit, symmetric matrix.}
#'  \item{\code{r} : }{ radius; a positive number.}
#'  \item{\code{q} : }{ cross term \code{-Pc}; a vector.}
#'  \item{\code{L} : }{ Cholesky decomposition of \code{P}}
#'  \item{\code{e} : }{ eigen values of \code{P}; a vector of eigenvalues.}
#'  \item{\code{U} : }{ eigen vectors of \code{P}; a matrix, where each column
#'    is an eigen vector.}
#'  \item{\code{D} : }{ diagnonal matrix with \code{sqrt(e)} as diagonal entries.}
#' }
#'
#' An ellipse \code{E} may be fully parameterized using the above parameters in
#' the following ways:
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  (x - c)<sup>T</sup> P (x - c) &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: (x-c)^T P (x-c) \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  (x - c)<sup>T</sup> S<sup>-1</sup> (x - c) &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: (x-c)^T S^{-1} (x-c) \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  (x - c)<sup>T</sup> LL<sup>T</sup> (x - c) &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: (x-c)^T LL^T (x-c) \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  (x - c)<sup>T</sup> UD<sup>2</sup>U<sup>T</sup> (x - c) &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: (x-c)^T UD^2U^T (x-c) \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  || L<sup>T</sup>(x - c) ||<sub>2</sub> &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: ||L^T(x-c)||_2 \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  || DU<sup>T</sup>(x - c) ||<sub>2</sub> &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: ||DU^T(x-c)||_2 \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  c + L<sup>-T</sup>w, ||w|| &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: c + L^{-T}w, \ ||w|| \leq r\}}}
#' \ifelse{html}{\out{<center>E = { x &#8712 R<sup>p</sup> :  c + UD<sup>-1</sup>w, ||w|| &#8804 r }</center>}}{\deqn{E = \{x \in \R^p: c + UD^{-1}w, \ ||w|| \leq r\}}}
#'
#' To ensure that all of the above parameters can be calculated it is advised
#' (but in some cases not needed) that the input ellipses are fully
#' parameterized.
#'
#'
#' @param ell a list of (non degenerate) ellipses to be wrangled. An ellipse is
#'   a named list and each entry corresponds to a parameter. To ensure all
#'   \code{out_params} can be calculated one of the parametrizations listed
#'   below in the description must be specified. Some \code{out_params} do
#'   not require a fully parametrized ellipse and so partially specified
#'   ellipses can be used.
#'
#' @param out_params a vector of names of the output parameters. A list of
#'   possible parameters is given below in the details.
#'
#' @return A list of wrangled ellipses. The wrangled ellipses are now given by
#'   the \code{out_params}.
#'
#'
#' @examples
#' # two dimensional unite ball
#' e2d <- list(c = c(0,0), S = matrix(c(1,0,0,1), ncol = 2), r = 1)
#'
#' # three dimensional ellipse
#' e3d <- list(c = c(3,2,1), P = matrix(c(3,1,2,1,5,0,2,0,2), ncol = 3))
#'
#' f1 <- wrangle_ellipse(e2d) # (c,P,r) parameterization
#' f2 <- wrangle_ellipse(e2d, out_params = c("c", "e", "U", "r"))
#' f3 <- wrangle_ellipse(list("ellipse1" = e2d, "ellipse2" = e3d),
#'                       c("n", "c", "U", "D"))
#'
#' @export
wrangle_ellipse <- function(ell, out_params = c("c", "P", "r")) {
  pos <- c("n", "c", "P", "S", "q", "r", "L", "e", "D", "U") # these need to correspond to find_x functions !
  if (!is.list(ell)) {
    stop("'ell' must be a named list or a list of named lists.")
  } else if (any(sapply(ell, function(e) {! is.list(e)}))) {
      ne <- 1
      ell <- list(ell)
  } else {
    ne <- 2
  }
  if (!is.character(out_params)) {
    stop("'out_params' must be a vector of characters")
  } else if (any(! out_params %in% pos)) {
    stop("One or more characters from 'out_prams' are not recognised.")
  }

  if (! "n" %in% out_params) {
    sapply(ell, function(e) {find_n(e)})
  }

  new_ell <- lapply(ell, function(e) {
    output <- lapply(out_params, function(i) {
      if (is.null(e[[i]])) { i } else { e[[i]] }
    })
    names(output) <- out_params
    lapply(output, function(o) {
      if (!is.character(o)) { o } else {
        funname <- paste0("find_", o)
        foo <- get(funname, mode = "function")
        return(foo(e))
      }
    })
  })

  if (ne == 1) { # if input not list of lists
    new_ell <- unlist(new_ell, recursive = FALSE)
  }
  return(new_ell)
}



find_n <- function(e) {
  nn <- numeric(0)
  if (!is.null(e[["c"]])) { nn <- c(nn, length(e[["c"]])) }
  if (!is.null(e[["q"]])) { nn <- c(nn, length(e[["q"]])) }
  if (!is.null(e[["e"]])) { nn <- c(nn, length(e[["e"]])) }
  if (!is.null(e[["D"]])) { nn <- c(nn, dim(e[["D"]])) }
  if (!is.null(e[["U"]])) { nn <- c(nn, dim(e[["U"]])) }
  if (!is.null(e[["P"]])) { nn <- c(nn, dim(e[["P"]])) }
  if (!is.null(e[["S"]])) { nn <- c(nn, dim(e[["S"]])) }
  if (!is.null(e[["L"]])) { nn <- c(nn, dim(e[["L"]])) }
  if (all(nn[1] == nn)) {
    return(nn[1])
  } else {
    stop("The different parameters imply different dimensions.")
  }
}

find_c <- function(e) {
  if (! is.null(e[["c"]])) {
    return(e[["c"]])
  }
  stop("It is not possible to determin 'P' from the given ellipse")
}

find_P <- function(e) {
  if (! is.null(e[["P"]])) {
    return(e[["P"]])
  } else if (! is.null(e[["S"]])) {
    return(solve(e[["S"]]))
  } else if (! is.null(e[["L"]])) {
    return(e[["L"]] %*% t(e[["L"]]))
  } else if (! (is.null(e[["U"]]) | is.null(e[["D"]]))) {
    return(e[["U"]] %*% e[["D"]] %*% e[["D"]] %*% t(e[["U"]]))
  } else if (! (is.null(e[["e"]]) | is.null(e[["U"]]))) {
    return(e[["U"]] %*% diag(e[["e"]]) %*% t(e[["U"]]))
  }
  stop("It is not possible to determin 'P' from the given ellipse")
}

find_S <- function(e) {
  if (! is.null(e[["S"]])) {
    return(e[["S"]])
  } else if (! is.null(e[["P"]])) {
    return(solve(e[["P"]]))
  } else if (! is.null(e[["L"]])) {
    return(solve(e[["L"]] %*% t(e[["L"]])))
  } else if (! (is.null(e[["e"]]) | is.null(e[["U"]]))) {
    return(e[["U"]] %*% diag(1 / e[["e"]]) %*% t(e[["U"]]))
  } else if (! (is.null(e[["D"]]) | is.null(e[["U"]]))) {
    return(e[["U"]] %*% diag(1 / (diag(e[["D"]])^2)) %*%t(e[["U"]]))
  }
  stop("It is not possible to determin 'S' from the given ellipse")
}

find_r <- function(e) {
  if (!is.null(e[["r"]])) {
    return(e[["r"]])
  }
  stop("It is not possible to determin 'r' from the given ellipse")
}

find_e <- function(e) {
  if (! is.null(e[["e"]])) {
    return(e[["e"]])
  } else if (!is.null(e[["P"]])) {
    return(eigen(e[["P"]], only.values = TRUE)$values)
  } else if (! is.null(e[["S"]])) {
    return(1 / eigen(e[["S"]], only.values = TRUE)$values)
  } else if (! is.null(e[["D"]])) {
    return(diag(e[["D"]])^2)
  } else if (! is.null(e[["L"]])) {
    return(eigen(e[["L"]]%*%t(e[["L"]]), only.values = TRUE)$values)
  }
  stop("It is not possible to determin 'e' from the given ellipse")
}


find_D <- function(e) {
  if (!is.null(e[["D"]])) {
    return(e[["D"]])
  } else if (!is.null(e[["e"]])) {
    return(sqrt(diag(e[["e"]])))
  } else if (!is.null(e[["P"]])) {
    return(diag(sqrt(eigen(e[["P"]], only.values = TRUE)$values)))
  } else if (!is.null(e[["S"]])) {
    return(return(diag(1/sqrt(eigen(e[["S"]], only.values = TRUE)$values))))
  } else if (! is.null(e[["L"]])) {
    return(diag(sqrt(eigen(e[["L"]]%*%t(e[["L"]]), only.values = TRUE)$values)))
  }
  stop("It is not possible to determin 'D' from the given ellipse")
}


find_U <- function(e) {
  if (!is.null(e[["U"]])) {
    return(e[["U"]])
  } else if (!is.null(e[["P"]])) {
    return(eigen(e[["P"]])$vectors)
  } else if (!is.null(e[["S"]])) {
    return(eigen(e[["S"]])$vectors)
  } else if (! is.null(e[["L"]])) {
    return(eigen(e[["L"]]%*%t(e[["L"]]))$vectors)
  }
  stop("It is not possible to determin 'U' from the given ellipse")
}




find_q <- function(e) {
  if (is.null(e[["c"]])) {
    stop("It is not possible to determin 'q' from the given ellipse")
  }
  if (! is.null(e[["P"]])) {
    return(- 2 * (e[["P"]] %*% e[["c"]]))
  } else if (! is.null(e[["S"]])) {
    return(- 2 * (solve(e[["S"]]) %*% e[["c"]]))
  } else if (! (is.null(e[["D"]]) | is.null(e[["U"]]))) {
    return(- 2 * e[["U"]]%*%e[["D"]]%*%e[["D"]]%*%t(e[["U"]]) %*% e[["c"]])
  } else if (! (is.null(e[["e"]]) | is.null(e[["U"]]))) {
    return(- 2 * e[["U"]]%*%diag(e[["e"]])%*%t(e[["u"]]) %*% e[["c"]])
  } else if (! is.null(e[["L"]])) {
    return(- 2 * e[["L"]]%*%t(e[["L"]]) %*% e[["c"]])
  }
  stop("It is not possible to determin 'q' from the given ellipse")
}



find_L <- function(e) {
  if (!is.null(e[["L"]])) {
    return(e[["L"]])
  } else if (!is.null(e[["P"]])) {
    return(t(chol(e[["P"]])))
  } else if (!is.null(e[["S"]])) {
    return(t(chol(solve(e[["S"]]))))
  } else if (! (is.null(e[["D"]]) | is.null(e[["U"]]))) {
    return(t(chol(e[["U"]]%*%e[["D"]]%*%e[["D"]]%*%t(e[["U"]]))))
  } else if (! (is.null(e[["e"]]) | is.null(e[["U"]]))) {
    return(t(chol(e[["U"]]%*%diag(e[["e"]])%*%t(e[["U"]]))))
  }
  stop("It is not possible to determin 'L' from the given ellipse")
}

