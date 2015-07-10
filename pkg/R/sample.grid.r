# Purpose        :  get a subset of spatial points with at least n points in each cell;
# Maintainer     : Wei Shangguan (shanggv@hotmail.com);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl);
# Status         : tested
# Note           :
#


sample.grid.SpatialPointsDataFrame <- function(obj, cell.size, n, bbox, ...){
    #to avoid the gid in obj
    if ("gid" %in% names(obj)){
        gid <- obj$gid
    }else{
        gid <- NULL
    }
    obj$gid <- 1:length(obj)
    if(missing(bbox)) { bbox <- obj@bbox }
    if(missing(cell.size)){
        # print warning:
        if(length(obj)>10000){
        warning("You need to assgn cell.size", immediate. = TRUE)
        }
    }
    require(raster)
    x <- GridTopology(cellcentre.offset=bbox[,1],
            cellsize=cell.size,
            cells.dim=c(ceiling(abs(diff(bbox[1,])/cell.size[1])),
            ncols=ceiling(abs(diff(bbox[2,])/cell.size[2]))))
    r.sp <- SpatialGrid(x, proj4string = obj@proj4string)
    r <- raster(r.sp)
    in.r <- rasterize(obj, r, field = "gid")
    grd <- as(in.r, "SpatialGridDataFrame")
    ov <- over(obj, grd, ...)
    pnts <- tapply(obj$gid, ov, function(x, n){
        if(length(x) <= n) x
        else sample(x,n)
        }, n=n)
    pnts <- unlist(pnts, use.names = FALSE)
    ret <- obj[obj$gid %in% pnts, ]
    if(is.null(gid)){
        ret$gid <- NULL
    }else{
        ret$gid <- gid
    }
    return(ret)
}
sample.grid.SpatialPoints <- function(obj, cell.size, n, bbox, ...){
    obj <- SpatialPointsDataFrame(obj, data.frame(id = 1:length(obj)))
    ret <- sample.grid.SpatialPointsDataFrame(obj, cell.size, n, bbox, ...)
    ret <- as(ret, "SpatialPoints")
    return(ret)
}
sample.grid <- function(obj, cell.size, n, bbox, ...) UseMethod("sample.grid")
setMethod("sample.grid", "SpatialPoints", sample.grid.SpatialPoints)
setMethod("sample.grid", "SpatialPointsDataFrame", sample.grid.SpatialPointsDataFrame)
