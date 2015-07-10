\name{sample.grid}
\docType{sample.grid}
\alias{sample.grid}
\alias{sample.grid,SpatialPoints,SpatialPoints-method}
\alias{sample.grid,SpatialPointsDataFrame,SpatialPoints-method}
\title{sample spatial points by grids}
\description{Get a subset of a object of class \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"} avoiding spatial clustered distribution. Spatial points are overlayed with spatial grids with a specified cell size and then get a subset from each grid with a specified number at most. If one grid has less points than the specified number, all the points are taken. If one grid has more points than the specified number, only this number of points are taken by \code{\link{sample}}.
}
\usage{
\S4method{merge}{SpatialPoints,SpatialPointsDataFrame}(obj, cell.size, n, bbox, \dots)
}
\arguments{
  \item{obj}{Spatial-Points,-PointsDataFrame object}
  \item{cell.size}{the cell size of the overlayed \code{"SpatialGridDataFrame"}}in the form of c(x,y)
  \item{n}{the maximum number points in each grid}
  \item{bbox}{the bounding box of output \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"}; it is set the same as the \code{obj} if missing}  
  \item{\dots}{other optional arguments that can be passed to \code{over}}
}
\note{This function can be used when there are too much point observations to be handled, especially for spatially clustered observations. The total number of sampled points are determined by \code{cell.size} and \code{n} together. You will get fewer the sampled points when \code{cell.size} is larger, or/and when \code{n} is smaller. Similar sample sizes can be achieved by differen combination of \code{cell.size} and \code{n}.}
\value{
Returns an object of type \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"} that contains a subset of the {obj}.
}
\author{Wei Shangguan}
\examples{
library(sp)
data(isis)
profs <- isis[["sites"]]
coordinates(profs) <- ~  LONWGS84 + LATWGS84
proj4string(profs) <- CRS("+proj=longlat +datum=WGS84")
#sample SpatialPointsDataFrame: 
prof1 <- sample.grid(profs, cell.size = c(1,2), n = 5)
str(prof1)
#sample SpatialPoints:
profs <- as(profs, "SpatialPoints")
prof2 <- sample.grid(profs, cell.size = c(1,1), n = 2)
str(prof2)
}
\keyword{methods}