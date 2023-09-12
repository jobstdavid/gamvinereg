#' gives overview over fitted D-vine GAM copula regression model
#' @noRd
#' @export
print.gamvinereg <- function(x, ...) {
  cat("D-vine regression model: ")
  n_predictors <- length(x$order)
  if (n_predictors <= 10) {
    predictors <- paste(x$order, collapse = ", ")
  } else {
    predictors <- paste(x$order[1:10], collapse = ", ")
    predictors <- paste0(predictors, ", ... (", n_predictors - 10, " more)")
  }
  cat(names(x$model_frame)[1], "|", predictors, "\n")
  cat("GAM copula regression model: ")
  if (x$vine$pair_copulas[[1]][[1]]$res@tau) {
    cat(paste0("tau | ", as.character(x$gam_formula)[2], sep = ""), "\n")
  } else {
    cat(paste0("par | ", as.character(x$gam_formula)[2], sep = ""), "\n")
  }
  stats <- unlist(x$stats[1:5])
  stats <- paste(names(stats), round(stats, 2), sep = " = ")
  cat(paste(stats, collapse = ", "), "\n")
  invisible(x)
}
#'
#' gives summary over fitted D-vine GAM copula regression model
#' @noRd
#' @export
summary.gamvinereg <- function(object, ...) {
  data.frame(
    var = c(names(object$model_frame)[1], object$order),
    edf = object$stats$var_edf,
    cll = object$stats$var_cll,
    caic = object$stats$var_caic,
    cbic = object$stats$var_cbic,
    p_value = object$stats$var_p_value
  )
}
#'
#'
