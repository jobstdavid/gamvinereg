#' Predict quantiles from D-vine GAM copula regression model
#'
#' @param object an object of class \code{gamvinereg}.
#' @param newdata data frame of variables for which to predict the quantiles.
#' @param alpha vector of quantile levels.
#' @param cores integer; the number of cores used for computations.
#' Default setting is \code{cores = 1}.
#' @param ... unused.
#'
#' @return a data frame of quantiles where each column corresponds to one value of alpha.
#'
#' @examples
#' # load data for station Hannover
#' data(station)
#'
#' # formula for gamvinereg
#' formula <- obs ~ t2m_mean + t2m_sd
#'
#' # data transformation to uniform scale
#' ustation <- as.data.frame(VineCopula::pobs(station[, all.vars(formula)]))
#'
#' ## generate margin distributions
#' # generate left-truncated at 0 normal distribution
#' gen.trun(par = 0, family = "NO", type = "left")
#' # generate left-censored normal distribution
#' gen.cens(family = "NO", type = "left")
#' # generate log-normal distribution
#' gen.Family(family = "NO", type = "log")
#'
#' ## margin specifications to select from
#' margins <- list(
#'
#'  # margin for variable obs (normal distribution)
#'  list(list(family = NO(),
#'            mu.formula = obs ~ sin1 + cos1,
#'            sigma.formula = ~ sin1 + cos1)),
#'
#'  # margin for variable t2m_mean (normal distribution)
#'  list(list(family = NO(),
#'            mu.formula = t2m_mean ~ sin1 + cos1,
#'            sigma.formula = ~ sin1 + cos1)),
#'
#'  # margins for variable t2m_sd (log-normal, left-truncated/censored at 0 normal distribution)
#'  list(list(family = logNO(),
#'            mu.formula = t2m_sd ~ sin1 + cos1,
#'            sigma.formula = ~ sin1 + cos1),
#'       list(family = NOtr(),
#'            mu.formula = t2m_sd ~ sin1 + cos1,
#'            sigma.formula = ~ sin1 + cos1),
#'       list(family = NOlc(),
#'            mu.formula = Surv(t2m_sd, t2m_sd > 0, type = "left") ~ sin1 + cos1,
#'            sigma.formula = ~ sin1 + cos1))
#'
#')
#'
#'
#' # fit gamvinereg with time-dependent linear correlation model
#' (object1 <- gamvinereg(formula = formula,
#'                        data = station,
#'                        control = gam_control(formula = ~ sin1 + cos1),
#'                        margins = margins))
#'
#' # fit gamvinereg with time-dependent spline based correlation model
#' (object2 <- gamvinereg(formula = formula,
#'                        data = station,
#'                        control = gam_control(formula = ~ s(doy, bs = "cc")),
#'                        margins = margins))
#'
#' # fit gamvinereg with on uniform data with constant correlation model
#' # fixed variable order t2m_sd - t2m_mean
#' (object3 <- gamvinereg(formula = formula,
#'                        data = ustation,
#'                        order = c("t2m_sd", "t2m_mean"),
#'                        uscale = TRUE))
#'
#' # predict median
#' med <- predict(object1, newdata = station, alpha = 0.5)
#'
#' # predict deciles
#' dec <- predict(object2, newdata = station, alpha = 1:9/10)
#'
#' # predict quartiles
#' q <- predict(object3, newdata = ustation, alpha = c(0.25, 0.75))
#'
#'
#'
#' @importFrom stats complete.cases
#' @export
predict.gamvinereg <- function(object, newdata, alpha = 0.5, cores = 1, ...) {

  # criterion for checking if data is already in [0,1]
  uscale <- is.na(object$stats$var_cll[1])

  if (!uscale) {

    margins_formulas <- unlist(object$margins$specifications, recursive = FALSE)
    margins_predictors <- unlist(as.vector(sapply(1:length(margins_formulas), function(k) c(all.vars(margins_formulas[[k]]$mu.formula)[-1],
                                                                                            all.vars(margins_formulas[[k]]$sigma.formula),
                                                                                            all.vars(margins_formulas[[k]]$nu.formula),
                                                                                            all.vars(margins_formulas[[k]]$tau.formula)))))
    control_predictors <- all.vars(object$gam_formula)
    additional_predictors <- c(margins_predictors, control_predictors)
    dummy_formula <- update(object$formula, reformulate(c(".", additional_predictors)))

  } else {

    control_predictors <- all.vars(object$gam_formula)
    dummy_formula <- update(object$formula, reformulate(c(".", control_predictors)))

  }


  # create output data.frame
  output <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = length(alpha), dimnames = list(c(), alpha)))
  # subset data frame
  if (!missing(newdata)) {
    mf <- model.frame(dummy_formula, newdata, na.action = NULL)
    cc <- complete.cases(mf[, -1])
    mf <- mf[cc, ]
  } else {
    mf <- model.frame(dummy_formula, parent.frame(), na.action = NULL)
  }

  # GAM copula variables
  x <- mf[, all.vars(object$gam_formula), drop = FALSE]

  # to u-scale
  u <- get_pits(olddata = object$model_frame,
                newdata = mf,
                margins = object$margins$models[object$order],
                uscale = uscale,
                cores = cores)

  # calculate quantile on uniform (u)-scale
  q <- qdvine(u = u,
              x = x,
              vine = object$vine,
              alpha = alpha,
              cores = cores)

  if (!uscale) {
    # transform to original y-scale
    q <- to_yscale(model = object$margins$models[[1]],
                   olddata = object$model_frame,
                   newdata = mf,
                   q = q)
  }

  # fill in available quantiles
  output[cc, ] <- q

  return(output)
}
#'
#'
#' Following function are helper functions for predict.gamvinereg:
#'
#' calculate PIT values
#' @noRd
#' @importFrom gamlss predictAll
get_pits <- function(olddata, newdata, margins, uscale, cores = 1) {

  d <- length(margins)
  vars <- names(margins)

  if (!uscale) {

    u <- mclapply(1:d, function(k) {

      m <- margins[[k]]
      y <- as.numeric(newdata[, vars[k]])
      # estimate in-sample parameters
      if(class(m$y) == "Surv") {
        par <- lapply(1:length(m$parameters), function(k) predict(object = m,
                                                                  newdata = newdata,
                                                                  data = olddata,
                                                                  what = m$parameters[k],
                                                                  type = "response"))
        assign(vars[k], y)
        y <- eval(m$mu.formula[[2]])
        par <- c(list(y), par)
        names(par) <- c("q", m$parameters)

      } else {
        par <- predictAll(object = m,
                          newdata = newdata,
                          data = olddata,
                          type = "response")[m$parameters]
        par <- c(list(y), par)
        names(par) <- c("q", m$parameters)
      }
      # get in-sample u-data of y
      u <- as.numeric(do.call(paste0("p", m$family[1]), par))
      u

    }, mc.cores = cores)

    u <- matrix(unlist(u), ncol = d, byrow = FALSE)

  } else {

    # check if data in [0,1]
    u <- newdata[, vars]
    if ( any(u > 1) | any(u < 0) ) {
      stop("Data has to be in the interval [0,1]!")
    }

  }

  u <- matrix(as.matrix(truncate_u(u)), ncol = d)
  colnames(u) <- vars
  u <- as.data.frame(u)

  return(u)

}
#'
#' calculate quantiles from conditional D-vine GAM copula
#' @noRd
#' @importFrom VineCopula BiCopHinv2
qdvine <- function(u, x, vine, alpha, cores = 1) {
  d <- ncol(vine$structure)
  if (ncol(u) != d - 1)
    stop("Dimensions of u and vine are not compatible")

  ## obtain diagonal entries in V matrix
  n <- nrow(u)
  V <- array(NA, dim = c(d, d, n))
  V[d, -1, ] <- t(u)
  V2 <- V
  if (d > 2) {
    for (j in (d - 1):2) {
      for (k in (d - 1):j) {

        # temporär data
        u_e <- data.frame(u = cbind(V2[k + 1, j, ], V[k + 1, j + 1, ]), x)

        # bivariate copula info
        bc <- vine$pair_copulas[[d - k]][[j]]$res
        # family
        family <- bc@family
        # parameters
        par <- as.numeric(gamBiCopPredict(bc, newdata = u_e, target = "par")$par)
        if(family == 2) {
          par2 <- bc@par2
        } else {
          par2 <- 0
        }

        # pseudo observations for next tree
        if (family %in% c(1:2, 5)) {
          # evaluation works component-wise, i.e. u_e[k, 1], u_e[k, 2], par[k], par2[k] (vectorized version)
          V[k, j, ] <- BiCopHfunc1(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
          V2[k, j, ] <- BiCopHfunc2(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
        } else if (family %in% c(301:304, 401:404)) {
          # update family
          family <- sapply(1:length(par), function(k) getFams(family, par[k]))
          # evaluation works component-wise, i.e. u_e[k, 1], u_e[k, 2], par[k], family[k] (vectorized version)
          V[k, j, ] <- BiCopHfunc1(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
          V2[k, j, ] <- BiCopHfunc2(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
        }

      }
    }
    tmp <- t(V[-1, 2, ])
  } else {
    tmp <- as.numeric(V[-1, 2, ])
  }

  # return quantiles
  out <- mclapply(1:length(alpha), function(j) {

    U <- cbind(alpha[j], tmp)
    for (i in 1:(d-1)) {

      # temporär data
      u_e <- data.frame(u = U[, i:(i+1)], x)

      # bivariate copula info
      bc <- vine$pair_copulas[[d - i]][[1]]$res
      # family
      family <- bc@family
      # parameters
      par <- as.numeric(gamBiCopPredict(bc, newdata = u_e, target = "par")$par)
      if(family == 2) {
        par2 <- bc@par2
      } else {
        par2 <- 0
      }

      # pseudo observations for next tree
      if (family %in% c(1:2, 5)) {
        U[, i+1] <- BiCopHinv2(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
      } else if (family %in% c(301:304, 401:404)) {
        # update family
        family <- sapply(1:length(par), function(k) getFams(family, par[k]))
        U[, i+1] <- BiCopHinv2(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
      }

    }
    U[, d]

  }, mc.cores = cores)

  out <- matrix(truncate_u(unlist(out)), ncol = length(alpha), byrow = FALSE)

  return(out)

}
#'
#' transform uniform data back to response scale
#' @noRd
to_yscale <- function(model, olddata, newdata, q) {
  if(class(model$y) == "Surv") {
    par <- lapply(1:length(model$parameters), function(k) predict(object = model,
                                                                  newdata = newdata,
                                                                  data = olddata,
                                                                  what = model$parameters[k],
                                                                  type = "response"))
    names(par) <- model$parameters
  } else {
    par <- predictAll(object = model,
                      newdata = newdata,
                      data = olddata)[model$parameters]
  }
  par <- lapply(1:ncol(q), function(k) c(p = list(q[, k]), par))
  y <- sapply(1:ncol(q), function(k) do.call(paste0("q", model$family[1]), par[[k]]))

  return(y)
}
#'
