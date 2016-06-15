# Purpose        : Automate prediction of soil properties / soil types;
# Maintainer     : Tomislav Hengl (tom.hengl@isric.org)
# Contributions  : ;
# Dev Status     : Pre-Alpha
# Note           : ;

setMethod("autopredict", signature(target = "SpatialPointsDataFrame", covariates = "SpatialPixelsDataFrame"), function(target, covariates, auto.plot=TRUE, ...){
  
  ## parent call:
  parent_call <- as.list(substitute(list(...)))[-1]
  ## TH: TO-DO estimate processing time
  
  ## predictive components:
  spc.fm <- as.formula(paste("~", paste(names(covariates), collapse = "+")))
  covariates <- spc(covariates, spc.fm)
  ## Model:
  fm <- as.formula(paste(names(target)[1], "~", paste(names(covariates@predicted), collapse = "+")))
  if(is.factor(target@data[,1])){
    m <- spmultinom(fm, target, covariates@predicted, ...)
    if(auto.plot==TRUE){ plotKML(m@predicted, folder.name=names(target)[1], file.name=paste0(names(target)[1], "_predicted.kml")) }
    return(m)
  }
  if(is.numeric(target@data[,1])){
    if(!any(names(parent_call) %in% "method")){
      method <- "ranger"
    }
    ## TH: TO-DO add ensemble predictions
    m <- fit.gstatModel(target, fm, covariates@predicted, method=method, ...)
    ## predict:
    p <- predict(m, covariates@predicted)
    if(auto.plot==TRUE){ plotKML(p, folder.name=names(target)[1], file.name=paste0(names(target)[1], "_predicted.kml")) }
    return(p)
  }

})


makePixels <- function(x, y, factors, pixel.size=1e2, t.dens=.2, min.dim=50, max.dim=2000, gdalwarp, sigma=1e3, show.progress=TRUE, area.poly=NULL, remove.files=TRUE){
  
  if(requireNamespace("spatstat", quietly = TRUE)&requireNamespace("maptools", quietly = TRUE)){
    
    if(is.na(proj4string(x))){ stop("Object 'x' missing SRS / proj4string") }
    if(missing(gdalwarp)){
      gdalwarp <- .programPath(utility="gdalwarp")
    }
    if(!class(x)=="SpatialPoints"){ stop("object of class 'SpatialPoints' expected") }
    if(any(!file.exists(y))){ stop("Missing files in 'y'") }
    if(missing(factors)){ factors = rep(FALSE, length(y)) }
    
    ## determine spatial domain:
    ncol = round(abs(diff(x@bbox[1,]))/pixel.size); nrow = round(abs(diff(x@bbox[2,]))/pixel.size)
    if(nrow>max.dim|ncol>max.dim){ stop("Grid too large. Consider increasing 'pixel.size'") }
    if(nrow<min.dim|ncol<min.dim){ stop("Grid too small. Consider reducing 'pixel.size'") }
    grid.owin <- owin(x@bbox[1,], x@bbox[2,], mask=matrix(TRUE,nrow,ncol))
    x.ppp <- spatstat::ppp(x=x@coords[,1], y=x@coords[,2], window=grid.owin)
    dens <- spatstat::density.ppp(x.ppp, sigma=sigma)
    dens <- maptools::as.SpatialGridDataFrame.im(dens)
    if(is.null(area.poly)){
      dx <- quantile(dens@data[,1], t.dens, na.rm=TRUE)
      dens@data[,1] <- ifelse(dens@data[,1]<dx, NA, dens@data[,1])
    } else {
      if(!class(area.poly)=="SpatialPolygons"){ stop("object of class 'SpatialPolygons' expected") }
      area.poly <- as(rasterize(area.poly, raster(dens), field=1), "SpatialGridDataFrame")
      dens@data[,1] <- ifelse(is.na(area.poly@data[,1]), NA, dens@data[,1])
    }
    ## generate mask map:
    dens <- as(dens[1], "SpatialPixelsDataFrame")
    if(length(dens)>1e6){ warning("Grid contains more than 1e6 elements. Consider reducing 'pixel.size'.") }
    ## copy grid parameters
    proj4string(dens) = proj4string(x)
    attr(dens@bbox, "dimnames")[[1]] = attr(x@bbox, "dimnames")[[1]]
    attr(dens@coords, "dimnames")[[2]] = attr(x@coords, "dimnames")[[2]]
    
    ## resample all layers to target grid:
    message(paste("Resampling values from", length(y), "rasters..."))
    if (show.progress) { pb <- txtProgressBar(min=0, max=length(y), style=3) }
    out.fname <- list(NULL)
    for(i in 1:length(y)){
      if(.Platform$OS.type == "windows"){
        fname <- normalizePath(y[i], winslash="/")
      } else {
        fname <- normalizePath(y[i], sep="/")
      }
      ## decide about whether to aggregate or downscale
      gi <- GDALinfo(fname)
      if(gi[["res.x"]]<pixel.size & !factors[i]){ 
        resample.method = "average"
      } else {
        if(gi[["res.x"]]>pixel.size & !factors[i]){ 
          resample.method = "cubicspline"
        } else {
          resample.method = "near"
        }
      } 
      out.fname[[i]] = set.file.extension(paste0("r_", basename(fname)), ".tif")
      if(!file.exists(out.fname[[i]])){ system(paste0(gdalwarp, ' ', fname, ' ', out.fname[[i]], ' -t_srs \"', proj4string(dens), '\" -co \"COMPRESS=DEFLATE\" -r \"', resample.method, '\" -tr ', pixel.size, ' ', pixel.size,' -te ', paste(as.vector(dens@bbox), collapse=" "))) }
      if (show.progress) { setTxtProgressBar(pb, i) }
    }
    if (show.progress) { close(pb) }
    ## round up numbers:
    dens@data[,1] <- (dens@data[,1]/max(dens@data[,1], na.rm=TRUE))*1e2
    for(i in 1:length(y)){
      dens@data[,i+1] <- rgdal::readGDAL(out.fname[[i]])$band1[dens@grid.index]
    }
    names(dens) = c("dens", unlist(out.fname))
    ## Assign factors where necessary
    for(i in 1:length(y)){
      if(factors[i]){ dens@data[,i+1] = as.factor(dens@data[,i+1]) }
    }
    ## Remove all layers without any variation:
    rm.v <- !sapply(dens@data, function(x){var(as.numeric(x), na.rm=TRUE)})>0
    if(any(rm.v)){
      dens <- dens[!rm.v]
    }
    if(remove.files==TRUE){
      unlink(unlist(out.fname))
    }
    return(dens)
  } else {
    stop("Missing packages 'maptools' and/or 'spatstat'")
  }
}
