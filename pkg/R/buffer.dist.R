# Purpose        : Derive distances to points;
# Maintainer     : Tomislav Hengl (tom.hengl@isric.org)
# Contributions  : ;
# Dev Status     : Pre-Alpha
# Note           : not recommended for big data;


setMethod("buffer.dist", signature(observations = "SpatialPointsDataFrame", predictionDomain = "SpatialPixelsDataFrame"), function(observations, predictionDomain, classes, width, ...){
  if(missing(width)){ width <- sqrt(areaSpatialGrid(predictionDomain)) }
  if(!length(classes)==length(observations)){ stop("Length of 'observations' and 'classes' does not match.") }
  s <- list(NULL)
  for(i in 1:length(levels(classes))){
    pnts <- observations[which(classes==levels(classes)[i]),1]@coords
    if(!is.null(pnts)&nrow(pnts)>0){
      r <- rasterize(pnts, y=raster(predictionDomain))
      s[[i]] <- distance(r, width=width, ...)
    }
  }
  s <- s[sapply(s, function(x){!is.null(x)})]
  s <- brick(s)
  s <- as(s, "SpatialPixelsDataFrame")
  s <- s[predictionDomain@grid.index,]
  return(s)
})