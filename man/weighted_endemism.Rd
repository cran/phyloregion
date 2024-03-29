% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WeightedEndemism.R
\name{weighted_endemism}
\alias{weighted_endemism}
\title{Measure the distribution of narrow-ranged or endemic species.}
\usage{
weighted_endemism(x)
}
\arguments{
\item{x}{A (sparse) community matrix.}
}
\value{
A data frame of species traits by site.
}
\description{
\code{weighted_endemism} is species richness inversely weighted
by species ranges.
}
\examples{
library(terra)
data(africa)
p <- vect(system.file("ex/sa.json", package = "phyloregion"))
Endm <- weighted_endemism(africa$comm)
m <- merge(p, data.frame(grids=names(Endm), WE=Endm), by="grids")
m <- m[!is.na(m$WE),]

plot(m, "WE", col = hcl.colors(20), type="continuous")

}
\references{
Crisp, M.D., Laffan, S., Linder, H.P. & Monro, A. (2001) Endemism in the
Australian flora. \emph{Journal of Biogeography} \strong{28}: 183–198.

Daru, B.H., Farooq, H., Antonelli, A. & Faurby, S. (2020) Endemism
patterns are scale dependent. \emph{Nature Communications} \strong{11}
: 2115.
}
\keyword{bioregion}
