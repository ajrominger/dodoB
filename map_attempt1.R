library(maps)
library(maptools)
library(rgeos)
library(sp)

foo <- map('world', fill = TRUE, plot = FALSE)
foo2 <- map2SpatialPolygons(foo, foo$names)
foo2 <- gBuffer(foo2, byid = TRUE, width = 0)
m <- gUnionCascaded(foo2)

m <- SpatialPolygons(list(
    Polygons(m@polygons[[1]]@Polygons[sapply(m@polygons[[1]]@Polygons, 
                                                   function(p) !p@hole & p@area > 0.1)
                                            ], ID=1)))
proj4string(m) <- '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

m2 <- spTransform(m, CRS('+proj=mill +lon_0=-150 +lat_0=0 +R=6371000 +units=m +no_defs'))

pdf('foo.pdf')
# map('world')
plot(m2)
dev.off()
