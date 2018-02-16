## mapping functions modified from stackoverflow qustion 10620862
## as answered by Josh O'Brien

library(sp)
library(maps)
library(maptools)
library(rgdal)
library(rgeos)

## Convert a "maps" map to a "SpatialLines" map
makeSLmap <- function() {
    llCRS <- CRS("+proj=longlat +ellps=WGS84")
    wrld <- map("world", interior = FALSE, plot=FALSE,
                xlim = c(-179, 179), ylim = c(-89, 89))
    wrld_p <- pruneMap(wrld, xlim = c(-179, 179))
    
    map2SpatialLines(wrld_p, proj4string = llCRS)
}

## Clip SpatialLines neatly along the antipodal meridian
sliceAtAntipodes <- function(SLmap, lon_0) {
    ## Preliminaries
    long_180 <- (lon_0 %% 360) - 180
    llCRS  <- CRS("+proj=longlat +ellps=WGS84")  ## CRS of 'maps' objects
    eqcCRS <- CRS("+proj=eqc")
    ## Reproject the map into Equidistant Cylindrical/Plate Caree projection 
    SLmap <- spTransform(SLmap, eqcCRS)
    
    ## Make a narrow SpatialPolygon along the meridian opposite lon_0
    L  <- Lines(Line(cbind(long_180, c(-89, 89))), ID="cutter")
    SL <- SpatialLines(list(L), proj4string = llCRS)
    SP <- gBuffer(spTransform(SL, eqcCRS), 10, byid = TRUE)
    
    ## Use it to clip any SpatialLines segments that it crosses
    ii <- which(gIntersects(SLmap, SP, byid=TRUE))
    
    # Replace offending lines with split versions
    # (but skip when there are no intersections (as, e.g., when lon_0 = 0))
    if(length(ii) > 0) { 
        SPii <- gDifference(SLmap[ii], SP, byid=TRUE)
        SLmap <- rbind(SLmap[-ii], SPii)  
    }
    
    return(SLmap)
}

## re-center, and clean up remaining streaks
recenterAndClean <- function(SLmap, lon_0) {
    llCRS <- CRS("+proj=longlat +ellps=WGS84")  ## map package's CRS
    newCRS <- CRS(paste("+proj=eqc +lon_0=", lon_0, sep=""))
    
    ## Recenter 
    SLmap <- spTransform(SLmap, newCRS)
    
    return(SLmap)
}

## Put it all together:
Recenter <- function(lon_0 = -100, grid=FALSE, ...) {                        
    SLmap <- makeSLmap()
    SLmap2 <- sliceAtAntipodes(SLmap, lon_0)
    
    out <- recenterAndClean(SLmap2, lon_0)
    return(out)
}

x <- readOGR('dodoB.kml', 'dodoB')
x@data <- data.frame(x@data[, 1, drop = FALSE],
                     project = c('amazon', 'atlanticForest', 'us-china', 
                                 'us-china', 'hawaii', 'palau'),
                     area = c(3931680, 719858, 1530108, 
                              1450601, 140406, 67469),
                     col = hsv(c(0, 0.09, 0.15, 0.15, 0.7, 0.8)), 
                     stringsAsFactors = FALSE)

m <- Recenter(170)
x <- spTransform(x, CRS(proj4string(m)))

pdf('fig_map.pdf', width = 10, height = 4)
layout(matrix(1:2, ncol = 2), widths = c(2, 1))

## the actual map
par(mar = rep(0.5, 4))
plot(m, xlim = c(-1e+07, 2e+07), col = 'gray', lwd = 0.75)
plot(x, col = colAlpha(x$col, 0.7), border = x$col, add = TRUE, lwd = 2)
polygon(c(2.1e+07, 1.5e+07, 2.1e+07), c(-3e+06, 4e+06, 7.5e+06), col = 'white', 
        border = 'white')

## the scale plot

aScaleMax <- tapply(x$area, x$project, mean)
aScaleMin <- c(10^5, 10^4, 10^2, 10, 10^5)
aResoMax <- aScaleMax / 10^3
aResoMin <- c(100, 100, 1, 1, 50)
aCol <- tapply(x$col, x$project, head, n = 1)

smoothPoly <- function(x, y, ...) {
    xpr <- exp(spline(1:length(x), log(x), n = 50)[[2]])
    ypr <- exp(spline(1:length(y), log(y), n = 50)[[2]])
    
    polygon(xpr, ypr, ...)
}

par(mar = c(4, 2.5, 2, 1.5) + 0.5, mgp = c(2, 0.5, 0))

plot(1, log = 'xy', axes = FALSE, type = 'n',
     xlim = c(10^1, 10^7), ylim = c(0.5, 10^6), 
     xlab = 'Spatial scale', ylab = 'Spatial resolution')

for(i in 1:length(aScaleMax)) {
    # polygon(c(aScaleMin[i], aScaleMax[i], aScaleMax[i]), 
    #         c(aResoMin[i], aResoMax[i], aScaleMax[i] / 10), 
    #         col = colAlpha(aCol[i], 0.5), border = aCol[i], lwd = 3)
    smoothPoly(c(aScaleMin[i], aScaleMax[i], aScaleMax[i], aScaleMin[i]), 
               c(aResoMin[i], aResoMax[i], aScaleMax[i] / 10, aResoMin[i]), 
               col = colAlpha(aCol[i], 0.4), border = aCol[i], lwd = 1)
}


for(i in 1:length(aScaleMax)) {
    # polygon(c(aScaleMin[i], aScaleMax[i], aScaleMax[i]), 
    #         c(aResoMin[i], aResoMax[i], aScaleMax[i] / 10), 
    #         col = colAlpha(aCol[i], 0.5), border = aCol[i], lwd = 3)
    smoothPoly(c(aScaleMin[i], aScaleMax[i], aScaleMax[i], aScaleMin[i]), 
               c(aResoMin[i], aResoMax[i], aScaleMax[i] / 10, aResoMin[i]), 
               border = aCol[i], lwd = 2)
}

axis(1, at = 10^c(1, 4, 7), labels = c(expression(10~km^2), 
                                       expression(10^4~km^2),
                                       expression(10^7~km^2)))
axis(2, at = 10^c(0, 3, 6), labels = c(expression(1~km^2), 
                                       expression(10^3~km^2),
                                       expression(10^6~km^2)))
box()

dev.off()
