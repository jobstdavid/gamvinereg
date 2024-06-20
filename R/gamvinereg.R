#' D-vine GAM copula based quantile regression
#'
#' D-vine GAM copula based quantile regression as described by Jobst, Möller and Groß (2023).
#'
#' @param formula a formula object, with the response on the left of an
#' \code{~} operator, and the covariates, separated by \code{+} operators,
#' on the right for the D-vine.
#' @param data data frame containing all variables for the D-vine GAM copula.
#' @param margins list of margin specification(s) for each variable in the same order as in the formula argument.
#' Each margin distribution is specified as list by the arguments \code{"family"}, \code{"mu.formula"}, \code{"sigma.formula"}
#' and if necessary \code{"nu.formula"} and \code{"tau.formula"} according to \code{\link{gamlss}}.
#' @param selcrit selection criterion for the margins and GAM copulas.
#' \code{"AIC"} (default), \code{"BIC"} or \code{"logLik"} are the possible choices.
#' @param control control parameters for the GAM copulas. The \code{\link{gam_control}} function
#' provides the default setting.
#' @param order the order of covariates in the D-vine, provided as vector of
#' variable names.
#' @param uscale logical; if \code{TRUE}, then \code{gamvinereg} assumes that the variables for the D-vine
#' are already in [0,1] and the \code{margins} argument will be ignored.
#' Default setting is \code{uscale = FALSE}.
#' @param cores integer; the number of cores used for computations.
#' Default setting is \code{cores = 1}.
#' @param ... unused.
#'
#' @return An object of class \code{gamvinereg}. It is a list containing the elements
#' \describe{
#'   \item{formula}{the formula used for the D-vine fit.}
#'   \item{gam_formula}{the formula used for the GAM copula fit.}
#'   \item{selcrit}{criterion used for margin and variable selection.}
#'   \item{model_frame}{the data used to fit the regression model.}
#'   \item{margins}{list consisting of the selected models, compared margin
#'   specifications and its corresponding statistics.}
#'   \item{vine}{a list consisting of the selected pair copulas, the D-vine
#'   strucuture and the variable names.}
#'   \item{stats}{fit statistics such as conditional
#'   log-likelihood/AIC/BIC and p-values for each variable's contribution.}
#'   \item{order}{order of the covariates chosen by the variable selection
#'   algorithm.}
#'   \item{selected_vars}{indices of selected variables.}}
#' Use the \code{\link{predict.gamvinereg}} function to predict conditional quantiles.
#' Furthermore, apply the \code{\link{cpit}} function to calculate the conditional
#' probability integral transform (PIT) values. Additionally, use the
#' \code{\link{cpdf}} function to calculate the conditional density values.
#' The \code{summary} function of a \code{gamvinereg} object shows the contribution of
#' each selected variable with the associated p-value derived from a likelihood ratio test.
#'
#' @references
#'
#' Kraus, D. and Czado, C. (2017). D-vine copula based quantile regression,
#' Computational Statistics and Data Analysis, 110, 1-18.doi: 10.1016/j.csda.2016.12.009.
#'
#' Jobst, D., Möller, A., and Groß, J. (2023). D-Vine GAM Copula based
#' Quantile Regression with Application to Ensemble Postprocessing. doi: 10.48550/arXiv.1510.04161.
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
#' med <- predict(object1, alpha = 0.5)
#'
#' # predict deciles
#' dec <- predict(object2, alpha = 1:9/10)
#'
#' # predict quartiles
#' q <- predict(object3, newdata = ustation, alpha = c(0.25, 0.75))
#'
#' @importFrom stats model.frame na.omit reformulate update
#' @export
gamvinereg <- function(formula, data, margins, selcrit = "AIC", control = gam_control(),
                       order = NA, uscale = FALSE, cores = 1, ...) {

  ## pre-processing
  if (!uscale) {

    # remove unused variables
    # create dummy formula for subsetting data frame for needed variables
    margins_formulas <- unlist(margins, recursive = FALSE)
    margins_predictors <- unlist(as.vector(sapply(1:length(margins_formulas), function(k) c(all.vars(margins_formulas[[k]]$mu.formula)[-1],
                                                                                            all.vars(margins_formulas[[k]]$sigma.formula),
                                                                                            all.vars(margins_formulas[[k]]$nu.formula),
                                                                                            all.vars(margins_formulas[[k]]$tau.formula)))))
    control_predictors <- all.vars(control$formula)
    additional_predictors <- c(margins_predictors, control_predictors)
    dummy_formula <- update(formula, reformulate(c(".", additional_predictors)))

  } else {

    # create dummy formula for subsetting data frame for needed variables
    control_predictors <- all.vars(control$formula)
    dummy_formula <- update(formula, reformulate(c(".", control_predictors)))
    margins <- NULL

  }

  # subset data frame
  if (!missing(data)) {
    mf <- model.frame(dummy_formula, data, na.action = na.omit)
  } else {
    mf <- model.frame(dummy_formula, parent.frame(), na.action = na.omit)
  }

  # marginal variables
  vars <- all.vars(formula)
  d <- length(vars)
  # GAM copula variables
  x <- mf[, all.vars(control$formula), drop = FALSE]

  # estimation of the margins and transformation to copula data
  margins <- fit_margins(data = mf,
                         vars = vars,
                         margins = margins,
                         selcrit = selcrit,
                         uscale = uscale,
                         cores = cores)
  u <- to_uscale(data = mf,
                 margins = margins,
                 uscale = uscale,
                 cores = cores)

  ## initialization
  current_fit <- initialize_fit(u = u)
  status <- initialize_status(margins = margins,
                              selcrit = selcrit,
                              uscale = uscale)

  ## estimation
  # no predetermined order
  if (any(is.na(order))) {

    # sequential forward variable selection algorithm
    for (i in seq_len(d - 1)) {
      new_fits <- mclapply(status$remaining_vars + 1, function(k) {
        xtnd_vine(new_var = u[, k],
                  old_fit = current_fit,
                  control = control,
                  selcrit = selcrit,
                  x = x,
                  cores = cores)
      }, mc.cores = cores)

      # update selection criteria and copulas
      status <- update_status(status = status,
                              new_vines = new_fits)

      # check if optimum is found
      if (status$optimum_found) {
        break
      }

      # select 'best' fit
      current_fit <- new_fits[[status$best_ind]]
    }

    # update selected variables
    names(status$selected_vars) <- vars[status$selected_vars + 1]

  } else { # predetermined order

    # check that only available variables are in the 'order'
    check_order(order = order,
                var_nms = vars)

    # to guarantee that all desired variables in 'order' will be included
    status$selcrit <- "keep_all"

    # conditional D-vine copula estimation using fixed order
    for (var in order) {
      current_fit <- xtnd_vine(u[, var], current_fit, control, selcrit, x, cores = cores)
      status <- update_status(status, list(current_fit))
    }

    # update selected variables
    status$selected_vars <- sapply(order, function(nm) which(vars == nm) - 1)
    status$remaining_vars <- numeric(0)

    # for finalize_vinereg_object correct copula selcrit
    status$selcrit <- selcrit

  }

  # prepare final gamvinereg object
  finalize_vinereg_object(
    formula = formula,
    control_formula = control$formula,
    model_frame = mf,
    margins = margins,
    vine = current_fit$vine,
    status = status,
    var_nms = vars
  )

}


