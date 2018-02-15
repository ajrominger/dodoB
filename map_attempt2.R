library(rgdal)
library(maps)
library(mapdata)

x <- readOGR('dodoB.kml', 'dodoB')
x@data <- data.frame(x@data[, 1, drop = FALSE],
                     project = c('amazon', 'atlanticForest', 'us-china', 
                                 'us-china', 'hawaii', 'palau'),
                     area = c(3931680, 719858, 1530108, 
                              1450601, 140406, 12753),
                     col = hsv(c(0, 0.09, 0.15, 0.15, 0.7, 0.8)), 
                     stringsAsFactors = FALSE)

par(mfrow = c(1, 2), mar = rep(0.1, 4))
plot(x, col = 'red')
box()
newlonglat <- '+proj=longlat +lon_0=-180 +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
x <- spTransform(x, CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=180 +o_lat_p=90 +lon_0=-10 +ellps=sphere +no_defs"),
                 use_ob_tran=TRUE)

plot(x, col = 'red')
box()

par(fg = 'transparent')
map('world', col = 'gray', fill = TRUE)
plot(x, add = TRUE, col = x$col)
axis(1)
