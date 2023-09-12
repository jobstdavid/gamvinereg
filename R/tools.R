utils::globalVariables(c("k", "RS"))
#' command above helps to suppress annoying R-check notes
#' for missing global variables
#'
#' ensures that data is not to close to 0 and 1
#' @noRd
truncate_u <- function(u) {
  pmin(pmax(u, 1e-10), 1 - 1e-10)
}
#'
#' estimates margin distributions via the R-package gamlss
#' @noRd
#' @import gamlss.cens gamlss.dist gamlss.tr
#' @importFrom stats logLik
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom parallel mclapply
fit_margins <- function(data, vars, margins, selcrit, uscale, cores = 1) {

  if (!uscale) {

    margins_tmp <- margins
    d <- length(margins)
    l <- sapply(1:length(margins), function(k) length(margins[[k]]))

    # check if margin specifications of each variable appear in same order as in the formula argument
    vars_check <- c()
    for (k in 1:d) {
      for (j in 1:l[k]) {
        vars_check <- c(vars_check, vars[k] == all.vars(margins[[k]][[j]]$mu.formula)[1])
      }
    }
    if (!all(vars_check)) {
      stop("Variables in the margin specifications are not in the same order as in the formula argument!")
    }

    margins <- unlist(margins, recursive = FALSE)
    n <- length(margins)
    data_tmp <- data

    fit <- mclapply(1:n, function(k) {

      # get GAMLSS settings
      formula <- margins[[k]]$mu.formula
      sigma.formula <- margins[[k]]$sigma.formula
      nu.formula <- margins[[k]]$nu.formula
      tau.formula <- margins[[k]]$tau.formula
      family <- margins[[k]]$family

      # GAMLSS
      gamlss(formula = formula,
             sigma.formula = sigma.formula,
             nu.formula = nu.formula,
             tau.formula = tau.formula,
             family = family,
             data = data_tmp,
             method = RS(),
             control = gamlss.control(trace = FALSE))

    }, mc.cores = cores)

    # save selcrit value, logLik and edf according to margin model structure
    ll <- sapply(1:length(margins), function(k) as.numeric(logLik(fit[[k]])))
    edf <- sapply(1:length(margins), function(k) attr(logLik(fit[[k]]), "df"))
    crit.val <- sapply(1:length(margins), function(k) do.call(selcrit, list(fit[[k]])))
    ll.sel <- edf.sel <- crit.val.sel <- list()
    i <- 1
    for (k in 1:d) {
      ll.sel[[k]] <- ll[i:(i-1+l[k])]
      edf.sel[[k]] <- edf[i:(i-1+l[k])]
      crit.val.sel[[k]] <- crit.val[i:(i-1+l[k])]
      i <- i+l[k]
    }

    # get index for marginal model for selcrit
    if (selcrit == "logLik") {
      sel.index <- sapply(1:d, function(k) which.max(crit.val.sel[[k]]))
    } else {
      sel.index <- sapply(1:d, function(k) which.min(crit.val.sel[[k]]))
    }

    i <- 1
    index <- c()
    for (k in 1:d) {
      s <- i:(i-1+l[k])
      index <- c(index, s[sel.index[k]])
      i <- i+l[k]
    }

    fit <- fit[index]
    names(fit) <- vars

    out <- list(models = fit,
                specifications = margins_tmp,
                stats = list(ll = ll.sel,
                             edf = edf.sel,
                             selcrit = crit.val.sel))
    names(out$stats)[3] <- selcrit

    gc()

  } else {

    # create empty dummy margins object
    d <- length(vars)
    l <- lapply(1:d, function(k) NA)
    out <- list(models = l,
                specifications = l,
                stats = list(ll = l,
                             edf = l,
                             selcrit = l))
    names(out$models) <- vars
    names(out$stats)[3] <- selcrit

  }

  return(out)

}
#'
#' transforms raw data to uniform data in [0,1]
#' @noRd
#' @importFrom stats predict
to_uscale <- function(data, margins, uscale, cores = 1) {

  if (!uscale) {

    d <- length(margins$models)

    u <- mclapply(1:d, function(k) {

      m <- margins$models[[k]]
      y <- m$y
      # estimate in-sample parameters
      par <- lapply(1:length(m$parameters), function(k) predict(object = m,
                                                                what = m$parameters[k],
                                                                type = "response"))
      par <- c(list(y), par)
      names(par) <- c("q", m$parameters)
      # get in-sample u-data of y
      u <- as.numeric(do.call(paste0("p", m$family[1]), par))
      u

    }, mc.cores = cores)

    u <- matrix(unlist(u), ncol = d, byrow = FALSE)
    # ensure that observations are not to close to 0 and 1
    u <- truncate_u(u)
    colnames(u) <- names(margins$models)

  } else {

    # check if data is in [0,1]
    vars <- names(margins$models)
    u <- data[, vars]
    if ( any(u > 1) | any(u < 0) ) {
      stop("Data has to be in the interval [0,1]!")
    }

    # ensure that observations are not to close to 0 and 1
    u <- truncate_u(u)
    colnames(u) <- vars

  }

  u <- as.data.frame(u)

  return(u)

}
#'
#' initializes estimation start for the D-vine GAM copula
#' @noRd
initialize_fit <- function(u) {
  list(
    # 1-dimensional (= empty) vine
    vine = list(pair_copulas = list(list()), matrix = as.matrix(1)),
    # array for storing pseudo-observations
    psobs = list(
      direct = array(u[, 1], dim = c(1, 1, nrow(u))),
      indirect = array(NA, dim = c(1, 1, nrow(u)))
    ),
    # conditional log-likelihood of the model
    cll = 0
  )
}
#'
#' initializes estimation status for the D-vine GAM copula
#' @noRd
initialize_status <- function(margins, selcrit, uscale) {

  if (!uscale) {
    ll <- logLik(margins$models[[1]])
    list(
      # remaining variable indices to select from
      remaining_vars = seq_len(length(margins$models) - 1),
      # variables indices included in the model
      selected_vars = NULL,
      # selection criterion
      selcrit = selcrit,
      # conditional log-liklihood (only unconditional margin so far)
      clls = as.numeric(ll),
      # number of parameters in current model
      edf = attr(ll, "df"),
      # TRUE when no improvement is possible
      optimum_found = FALSE
    )
  } else {
    list(
      # remaining variable indices to select from
      remaining_vars = seq_len(length(margins$models) - 1),
      # variables indices included in the model
      selected_vars = NULL,
      # selection criterion
      selcrit = selcrit,
      # conditional log-liklihood (only unconditional margin so far)
      clls = NA,
      # number of parameters in current model
      edf = NA,
      # TRUE when no improvement is possible
      optimum_found = FALSE
    )
  }
}
#'
#' updates estimation status for the D-vine GAM copula
#' @noRd
update_status <- function(status, new_vines) {
  p <- length(new_vines)
  clls <- sapply(1:p, function(k) new_vines[[k]]$cll)
  edf <- sapply(1:p, function(k) new_vines[[k]]$edf)
  n <- dim(new_vines[[1]]$psobs$direct)[3]
  crits <- calculate_crits(clls, edf, n, status$selcrit)

  if (max(crits) <= 0) {
    # optimum found, keep old fit
    status$optimum_found <- TRUE
  } else {
    status$best_ind <- which.max(crits)
    status$selected_vars <- c(
      status$selected_vars,
      status$remaining_vars[status$best_ind]
    )
    status$remaining_vars <- setdiff(
      status$remaining_vars,
      status$selected_vars
    )
    status$clls = c(status$clls, clls[status$best_ind])
    status$edf = c(status$edf, edf[status$best_ind])
  }

  return(status)
}
#'
#' calculates selection criterion for sequential forward estimation
#' @noRd
calculate_crits <- function(clls, edf, n, selcrit) {
  clls - switch(
    selcrit,
    "logLik" = 0,
    "AIC" = edf,
    "BIC" = edf * log(n) / 2,
    "keep_all" = -Inf
  )
}
#'
#' extends D-vine GAM copula via sequential forward estimation
#' @noRd
#' @importFrom VineCopula BiCopHfunc1 BiCopHfunc2
#' @importFrom gamCopula gamBiCopFit gamBiCopPredict
#' @importFrom stats4 logLik
xtnd_vine <- function(new_var, old_fit, control, selcrit = "AIC", x, cores = 1) {
  d <- dim(old_fit$psobs$direct)[1] + 1
  n <- length(new_var)

  psobs <- list(
    direct = array(NA, dim = c(d, d, n)),
    indirect = array(NA, dim = c(d, d, n))
  )
  psobs$direct[-1, -d, ] <- old_fit$psobs$direct
  psobs$indirect[-1, -d, ] <- old_fit$psobs$indirect
  psobs$direct[d, d, ] <- new_var
  old_fit$vine$pair_copulas[[d - 1]] <- list()

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

    # estimate bivariate GAM copula
    pc_fit <- mclapply(1:length(control$family_set), function(k)  {

      gamBiCopFit(data = u_e,
                  formula = control$formula,
                  family = control$family_set[k],
                  tau = control$tau,
                  method = control$method,
                  tol.rel = control$tol.rel,
                  n.iters = control$n.iters,
                  verbose = FALSE)

    }, mc.cores = cores)

    val.selcrit <- sapply(1:length(pc_fit), function(k) do.call(selcrit, list(pc_fit[[k]]$res)))
    if (selcrit == "logLik") {
      index <- which.max(val.selcrit)
    }  else {
      index <- which.min(val.selcrit)
    }
    pc_fit <- pc_fit[[index]]


    old_fit$vine$pair_copulas[[d - i]][[i]] <- pc_fit

    # bivariate copula info
    # family
    family <- pc_fit$res@family
    # parameters
    par <- as.numeric(gamBiCopPredict(pc_fit$res, target = "par")$par)
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

  # some info about last selected bivariate gam copula
  cll <- logLik(pc_fit$res)
  edf <- attr(cll, "df")
  # some info about vine copula extension
  vine <- list(pair_copulas = old_fit$vine$pair_copulas,
               structure = gen_dvine_mat(d))

  # output
  list(vine = vine, psobs = psobs, cll = cll, edf = edf)

}
#'
#' selects appropriate copula according to the R-package VineCopula
#' @noRd
getFams <- function(family, par) {

  fam <- family

  if (par > 0 && fam %in% c(301, 302)) {
    return(3)
  } else if (par > 0 && fam %in% c(303, 304)) {
    return(13)
  } else if (par < 0 && fam %in% c(301, 303)) {
    return(23)
  } else if (par < 0 && fam %in% c(302, 304)) {
    return(33)
  } else if (par > 0 && fam %in% c(401, 402)) {
    return(4)
  } else if (par > 0 && fam %in% c(403, 404)) {
    return(14)
  } else if (par < 0 && fam %in% c(401, 403)) {
    return(24)
  } else if (par < 0 && fam %in% c(402, 404)) {
    return(34)
  } else {
    return(fam)
  }

  return(fam)
}
#'
#' creates D-vine structure according to the R-package VineCopula
#' @noRd
gen_dvine_mat <- function(d) {

  mat <- diag(d:1)
  for (j in 1:(d-1)) {
    mat[(j+1):d, j] <- 1:(d-j)
  }

  return(mat)

}
#'
#' checks fixed order in the D-vine
#' @noRd
check_order <- function(order, var_nms) {
  if (!all(order %in% var_nms))
    stop("unknown variable name in 'order'; ",
         "allowed values: '", paste(var_nms[-1], collapse = "', '"), "'.")
  if (any(order == var_nms[1]))
    stop("response variable '", var_nms[1],
         "' must not appear in 'order'.")
}
#'
#' creates final gamvinereg (S3) object
#' @noRd
#' @importFrom stats pchisq
finalize_vinereg_object <- function(formula, control_formula, model_frame, margins, vine, status, var_nms) {

  # adjust model matrix and names
  reorder <- status$selected_vars
  reorder[order(reorder)] <- seq_along(status$selected_vars)
  vine$names <- c(var_nms[1], names(reorder)[reorder])

  # compute fit statistics
  nobs <- nrow(model_frame)
  var_edf <- status$edf
  var_cll <- status$cll
  var_caic <- -2 * var_cll + 2 * var_edf
  var_cbic <- -2 * var_cll + log(nobs) * var_edf
  var_p_value <- pchisq(2 * var_cll, var_edf, lower.tail = FALSE)
  var_p_value[1] <- NA
  cll <- sum(var_cll)
  edf <- sum(var_edf)
  caic <- sum(var_caic)
  cbic <- sum(var_cbic)

  stats <- list(
    nobs = nobs,
    edf = edf,
    cll = cll,
    caic = caic,
    cbic = cbic,
    var_edf = var_edf,
    var_cll = var_cll,
    var_caic = var_caic,
    var_cbic = var_cbic,
    var_p_value = var_p_value
  )

  # return results as S3 object
  out <- list(
    formula = formula,
    gam_formula = control_formula,
    selcrit = status$selcrit,
    model_frame = model_frame,
    margins = margins,
    vine = vine,
    stats = stats,
    order = var_nms[status$selected_vars + 1],
    selected_vars = as.numeric(status$selected_vars)
  )
  class(out) <- "gamvinereg"

  return(out)

}
#'
#'

