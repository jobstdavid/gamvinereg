#' Conditional probability integral transform (PIT)
#'
#' Calculates the conditional distribution of the response given the covariates.
#'
#' @param object an object of class \code{gamvinereg}.
#' @param newdata data frame of response variable and its covariates.
#' @param cores integer; the number of cores used for computations.
#' Default setting is \code{cores = 1}.
#' @param ... unused.
#'
#' @return a vector of PIT values.
#'
#' @examples
#' # load data for station Hannover
#' data(station)
#'
#' # formula for gamvinereg
#' formula <- obs ~ t2m_mean + t2m_sd
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
#' (object <- gamvinereg(formula = formula,
#'                       data = station,
#'                       control = gam_control(formula = ~ sin1 + cos1),
#'                       margins = margins))
#'
#' # calculate PIT values
#' u <- cpit(object)
#' # should be approximately uniform
#' hist(u, prob = TRUE)
#'
#'
#' @export
cpit <- function(object, newdata, cores = 1, ...) {

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

  if (missing(newdata)) {
    mf <- model.frame(dummy_formula, object$model_frame, na.action = na.pass)
  } else {
    mf <- model.frame(dummy_formula, newdata, na.action = na.pass)
  }
  # create output vector
  output <- rep(NA, nrow(mf))
  # subset data frame
  cc <- complete.cases(mf[, -1])
  mf <- mf[cc, ]

  # GAM copula variables
  x <- mf[, all.vars(object$gam_formula), drop = FALSE]

  # to u-scale
  u <- get_pits(olddata = object$model_frame,
                newdata = mf,
                margins = object$margins$models[c(names(object$margins$models)[1], object$order)],
                uscale = uscale,
                cores = cores)
  psobs_old <-  list(direct = array(u[, 1], dim = c(1, 1, nrow(u))),
                     indirect = array(NA, dim = c(1, 1, nrow(u))))

  # conditional distribution
  cond_dist <- function(new_var, pair_copulas, psobs_old) {
    d <- dim(psobs_old$direct)[1] + 1
    n <- length(new_var)

    psobs <- list(
      direct = array(NA, dim = c(d, d, n)),
      indirect = array(NA, dim = c(d, d, n))
    )

    psobs$direct[-1, -d, ] <- psobs_old$direct
    psobs$indirect[-1, -d, ] <- psobs_old$indirect
    psobs$direct[d, d, ] <- new_var

    for (i in rev(seq_len(d - 1))) {
      # get data for current edge
      u_e <- matrix(NA, n, 2)
      u_e[, 1] <- psobs$direct[i + 1, i, ]
      u_e[, 2] <- if (i == d - 1) {
        psobs$direct[i + 1, i + 1, ]
      } else {
        psobs$indirect[i + 1, i + 1, ]
      }
      u_e <- data.frame(u = u_e, x)

      pc_fit <- pair_copulas[[d - i]][[i]]

      # bivariate copula info
      # family
      family <- pc_fit$res@family
      # parameters
      par <- as.numeric(gamBiCopPredict(pc_fit$res, newdata = u_e, target = "par")$par)
      if(family == 2) {
        par2 <- pc_fit$res@par2
      } else {
        par2 <- 0
      }

      # pseudo observations for next tree
      if (family %in% c(1:2, 5)) {

        # evaluation works component-wise, i.e. u_e[k, 1], u_e[k, 2], par[k], par2[k] (vectorized version)
        psobs$direct[i, i, ] <- BiCopHfunc2(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
        psobs$indirect[i, i, ] <- BiCopHfunc1(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)

      } else if (family %in% c(301:304, 401:404)) {

        # update family
        family <- sapply(1:length(par), function(k) getFams(family, par[k]))
        # evaluation works component-wise, i.e. u_e[k, 1], u_e[k, 2], par[k], family[k] (vectorized version)
        psobs$direct[i, i, ] <- BiCopHfunc2(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)
        psobs$indirect[i, i, ] <- BiCopHfunc1(u1 = u_e[, 1], u2 = u_e[, 2], family = family, par = par, par2 = par2)

      }

    }

    return(psobs)

  }

  for (var in object$order) {
    psobs_old <- cond_dist(new_var = u[, var],
                           pair_copulas = object$vine$pair_copulas,
                           psobs_old = psobs_old)
  }

  # get PIT values
  u <- psobs_old$direct[1, 1, ]

  # fill in available PIT values
  output[cc] <- u

  return(output)

}
