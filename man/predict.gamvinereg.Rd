% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.gamvinereg.R
\name{predict.gamvinereg}
\alias{predict.gamvinereg}
\title{Predict quantiles from D-vine GAM copula regression model}
\usage{
\method{predict}{gamvinereg}(object, newdata, alpha = 0.5, cores = 1, ...)
}
\arguments{
\item{object}{an object of class \code{gamvinereg}.}

\item{newdata}{data frame of variables for which to predict the quantiles.}

\item{alpha}{vector of quantile levels.}

\item{cores}{integer; the number of cores used for computations.
Default setting is \code{cores = 1}.}

\item{...}{unused.}
}
\value{
a data frame of quantiles where each column corresponds to one value of alpha.
}
\description{
Predict quantiles from D-vine GAM copula regression model
}
\examples{
# load data for station Hannover
data(station)

# formula for gamvinereg
formula <- obs ~ t2m_mean + t2m_sd

# data transformation to uniform scale
ustation <- as.data.frame(VineCopula::pobs(station[, all.vars(formula)]))

## generate margin distributions
# generate left-truncated at 0 normal distribution
gen.trun(par = 0, family = "NO", type = "left")
# generate left-censored normal distribution
gen.cens(family = "NO", type = "left")
# generate log-normal distribution
gen.Family(family = "NO", type = "log")

## margin specifications to select from
margins <- list(

 # margin for variable obs (normal distribution)
 list(list(family = NO(),
           mu.formula = obs ~ sin1 + cos1,
           sigma.formula = ~ sin1 + cos1)),

 # margin for variable t2m_mean (normal distribution)
 list(list(family = NO(),
           mu.formula = t2m_mean ~ sin1 + cos1,
           sigma.formula = ~ sin1 + cos1)),

 # margins for variable t2m_sd (log-normal, left-truncated/censored at 0 normal distribution)
 list(list(family = logNO(),
           mu.formula = t2m_sd ~ sin1 + cos1,
           sigma.formula = ~ sin1 + cos1),
      list(family = NOtr(),
           mu.formula = t2m_sd ~ sin1 + cos1,
           sigma.formula = ~ sin1 + cos1),
      list(family = NOlc(),
           mu.formula = Surv(t2m_sd, t2m_sd > 0, type = "left") ~ sin1 + cos1,
           sigma.formula = ~ sin1 + cos1))

)


# fit gamvinereg with time-dependent linear correlation model
(object1 <- gamvinereg(formula = formula,
                       data = station,
                       control = gam_control(formula = ~ sin1 + cos1),
                       margins = margins))

# fit gamvinereg with time-dependent spline based correlation model
(object2 <- gamvinereg(formula = formula,
                       data = station,
                       control = gam_control(formula = ~ s(doy, bs = "cc")),
                       margins = margins))

# fit gamvinereg with on uniform data with constant correlation model
# fixed variable order t2m_sd - t2m_mean
(object3 <- gamvinereg(formula = formula,
                       data = ustation,
                       order = c("t2m_sd", "t2m_mean"),
                       uscale = TRUE))

# predict median
med <- predict(object1, alpha = 0.5)

# predict deciles
dec <- predict(object2, alpha = 1:9/10)

# predict quartiles
q <- predict(object3, newdata = ustation, alpha = c(0.25, 0.75))



}
