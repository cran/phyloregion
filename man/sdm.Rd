% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdm.R
\name{sdm}
\alias{sdm}
\title{Species distribution models}
\usage{
sdm(
  x,
  predictors = NULL,
  pol = NULL,
  thin = TRUE,
  thin.size = 500,
  algorithm = "all",
  size = 50,
  width = 50000,
  mask = FALSE
)
}
\arguments{
\item{x}{A dataframe containing the species occurrences
and geographic coordinates. Column 1 labeled as "species", column 2 "lon",
column 3 "lat".}

\item{predictors}{A \code{SpatRaster} to extract values from the
locations in x on which the models will be projected.}

\item{pol}{A vector polygon specifying the boundary to restrict the
prediction. If \code{NULL}, the extent of input points is used.}

\item{thin}{Whether to thin occurrences}

\item{thin.size}{The size of the thin occurrences.}

\item{algorithm}{Character. The choice of algorithm to run the species
distribution model. Available algorithms include:
\itemize{
\item \dQuote{all}: Calls all available algorithms: GLM, and MAXENT.
\item \dQuote{GLM}: Calls only Generalized linear model.
\item \dQuote{MAXENT}: Calls only Maximum entropy.
}}

\item{size}{Minimum number of points required to successfully run
a species distribution model especially for species with few occurrences.}

\item{width}{Width of buffer in meter if x is in longitude/latitude CRS.}

\item{mask}{logical. Should y be used to mask? Only used if pol is a SpatVector}
}
\value{
A list with the following objects:
\itemize{
  \item \code{ensemble_raster} The ensembled raster that predicts
  the potential species distribution based on the algorithms selected.
  \item \code{data} The dataframe of occurrences used to implement the model.
  \item \code{polygon} Map polygons of the predicted distributions
  analogous to extent-of-occurrence range polygon.
  \item \code{indiv_models} Raster layers for the separate models that
  predict the potential species distribution.
}
}
\description{
This function computes species distribution models using
two modelling algorithms: generalized linear models,
and maximum entropy (only if \code{rJava} is available).
Note: this is an experimental function, and may change in the future.
}
\examples{
\donttest{
# get predictor variables
library(predicts)
f <- system.file("ex/bio.tif", package="predicts")
preds <- rast(f)
#plot(preds)

# get species occurrences
b <- file.path(system.file(package="predicts"), "ex/bradypus.csv")
d <- read.csv(b)

# fit ensemble model for four algorithms
m <- sdm(d, predictors = preds, algorithm = "all")
# plot(m$ensemble_raster)
# plot(m$polygon, add=TRUE)
}
}
\references{
Zurell, D., Franklin, J., König, C., Bouchet, P.J., Dormann, C.F., Elith, J.,
Fandos, G., Feng, X., Guillera‐Arroita, G., Guisan, A., Lahoz‐Monfort, J.J.,
Leitão, P.J., Park, D.S., Peterson, A.T., Rapacciuolo, G., Schmatz, D.R.,
Schröder, B., Serra‐Diaz, J.M., Thuiller, W., Yates, K.L., Zimmermann, N.E.
and Merow, C. (2020), A standard protocol for reporting species distribution
models. \emph{Ecography}, \strong{43}: 1261-1277.
}
