#' Auxiliary function for controlling the GAM copula fitting
#'
#' @param formula a formula object, starting with the \code{~} operator and
#' followed by the covariates, separated by \code{+} operators.
#' Default: \code{~ 1}, i.e. only intercept.
#' @param family_set a vector of bivariate copula families as implemented in
#' \code{\link{gamCopula}}, i.e. \code{1} Gaussian, \code{2} Student-t,
#' \code{301} Double Clayton type I (standard and rotated 90 degrees),
#' \code{302} Double Clayton type II (standard and rotated 270 degrees),
#' \code{303} Double Clayton type III (survival and rotated 90 degrees),
#' \code{304} Double Clayton type IV (survival and rotated 270 degrees),
#' \code{401} Double Gumbel type I (standard and rotated 90 degrees),
#' \code{402} Double Gumbel type II (standard and rotated 270 degrees),
#' \code{403} Double Gumbel type III (survival and rotated 90 degrees),
#' \code{404} Double Gumbel type IV (survival and rotated 270 degrees).
#' All available copula families are used per default.
#' @param tau logical; specifies whether Kendall's tau (default: \code{TRUE}) or the
#' copula parameter (\code{FALSE}) should be parametrized by the \code{formula}.
#' @param method character; Fisher-scoring (default: \code{"FS"}) or Newton-Raphson
#' (default: \code{"NR"}) for the fitting.
#' @param tol.rel double; relative tolerance for \code{"FS"/"NR"} algorithm.
#' Default: \code{tol.rel = 0.001}.
#' @param n.iters integer; maximal number of iterations for \code{"FS"/"NR"} algorithm.
#' Default: \code{n.iters = 10}.
#' @param ... unused.
#'
#' @return a list with components named as the arguments.
#'
#' @export
gam_control <- function(formula = ~ 1,
                        family_set = c(1, 2, 301:304, 401:404),
                        tau = TRUE,
                        method = "FS",
                        tol.rel = 0.001,
                        n.iters = 10,
                        ...) {
  list(formula = formula,
       family_set = family_set,
       tau = tau,
       method = method,
       tol.rel = tol.rel,
       n.iters = n.iters)
}


