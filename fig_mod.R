library(ape)
library(viridis)
library(socorro)
library(pika)
library(jpeg)
library(raster)
library(rgeos)

## function to read in image and make polygon output
arthOutline <- function(f, maxCol = 100) {
    a <- readJPEG(f)
    
    r <- as.matrix(as.raster(a))
    cols <- unique(as.vector(r))
    
    r <- raster(matrix((r %in% cols[colSums(col2rgb(cols)) < maxCol]) * 1, nrow = nrow(r)))
    
    l <- rasterToContour(r)
    x <- l@lines[[1]]@Lines[[1]]@coords
    
    return(x)
}

## funciton to add arthropod to a plot by name
addArth <- function(arth, x, y, width, col, alpha = 1) {
    d <- '~/Dropbox/Research/grants/macrosystems/figs'
    a <- switch(arth,
                'fly' = arthOutline(file.path(d, 'fly.jpg')),
                'moth' = arthOutline(file.path(d, 'moth.jpg'), 400),
                'beetle' = arthOutline(file.path(d, 'beetle.jpg')))
    
    yasp <- diff(range(par('usr')[3:4])) / diff(range(par('usr')[1:2]))
    if(arth %in% c('fly', 'beetle')) yasp <- yasp * 1.5
    
    a[, 2] <- a[, 2] * yasp
    
    asp <- width / diff(range(a[, 1]))
    a <- a * asp
    
    ai <- a
    for(i in 1:length(x)) {
        ai[, 1] <- a[, 1] - mean(a[, 1]) + x[i]
        ai[, 2] <- a[, 2] - mean(a[, 2]) + y[i]
        
        polygon(ai, col = rgb(t(col2rgb(col)), maxColorValue = 255, alpha = alpha*255), 
                border = ifelse(alpha == 1, 'black', col))
    }
}


## regional phylo
set.seed(5)
tre <- rphylo(5, 0.9, 0.8)
set.seed(2)
trt <- rTraitCont(tre)
# trt[tre$tip.label == 't2'] <- 2 * trt[tre$tip.label == 't2']
trt['t5'] <- -0.05
trt['t3'] <- -0.02


pdf('fig_mod_phylo.pdf', width = 3, height = 3)
par(mar = rep(1, 4))
plot(tre, show.tip.label = FALSE, edge.width = 3)

xTip <- max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)
ordTip <- c('t3', 't4', 't1', 't5', 't2')
colTrt <- quantCol(trt[ordTip], pal = viridis(10))

par(xpd = NA)
arth <- c('beetle', 'beetle', 'fly', 'fly', 'moth')
for(i in 1:length(arth)) {
    addArth(arth[i], x = 1.05 * xTip, y = i, width = ifelse(arth[i] == 'moth', 0.6, 0.4), col = colTrt[i])
}

dev.off()

## local comm
locPlot <- function(r, x0, y0, n0, trt, env, ...) {
    ## xy for dots
    coord <- seq(-r, r, length.out = ceiling(1.1 * sqrt(n0 / (pi*0.5^2))))
    xy <- expand.grid(coord, coord)
    xy <- xy[sqrt(xy[, 1]^2 + xy[, 2]^2) <= 1.05 * r, ]
    xy[, 1] <- xy[, 1] + x0
    xy[, 2] <- xy[, 2] + y0
    
    ## spp ID
    n1 <- sample(trt, size = nrow(xy), 
                 prob = dnorm(trt, mean(trt) + env * diff(range(trt)), 0.1), 
                 replace = TRUE)
    
    ## plotting
    rOut <- diff(range(xy[, 1])) / 2 + diff(xy[1:2, 1]) / 1.5
    
    polygon(x0 + 1.1*rOut * cos(seq(0, 2*pi, length.out = 100)), 
            y0 + 1.1*rOut * sin(seq(0, 2*pi, length.out = 100)), 
            col = colAlpha(quantCol(env, viridis(10), xlim = range(trt)), 0.2))
    
    n1 <- n1[(xy[, 1] - x0)^2 + (xy[, 2] - y0)^2 < rOut^2]
    
    ## jerry rig
    n1 <- n1[c(1, 4, 5, 2, 3)] 
    if(n1[1] == trt['t1']) {
        n1[1] <-trt['t3']
        names(n1)[1] <- 't3'
    }
    
    xy <- xy[(xy[, 1] - x0)^2 + (xy[, 2] - y0)^2 < rOut^2, ]
    
    arth1 <- match(names(n1), ordTip)
    col1 <- quantCol(n1, viridis(10), xlim = range(trt))
    
    for(i in 1:length(n1)) {
        addArth(arth[arth1[i]], x = xy[i, 1], y = xy[i, 2], width = ifelse(arth[arth1[i]] == 'moth', 0.55, 0.3), col = col1[i])
    }
    
    return(arth1)
}

pdf('fig_mod_comm.pdf', width = 3, height = 3)
par(mar = rep(0.1, 4))
plot(1, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, type = 'n', axes = FALSE)
set.seed(1)
abund1 <- locPlot(1.25, -0.5, 1, 15, trt, 0.9 * max(trt), cex = 1)
set.seed(2)
abund2 <- locPlot(1.25, 0.5, -1, 15, trt, 0.5 * min(trt), cex = 1)
dev.off()

## trait evol through time

rwalk <- function(n, x0, x1, sd) {
    X <- cumsum(c(0, rnorm(n - 1, 0, sd)))
    
    newX <- X - ((X[n] - X[1]) / (n - 1) * (1:n - 1) + X[1]) + 
        (x1 - x0) / (n - 1) * (1:n - 1) + x0
    
    newX[newX > max(x0, x1)] <- max(x0, x1)
    newX[newX < min(x0, x1)] <- min(x0, x1)
    
    return(newX)
}

set.seed(5)
x1 <- rwalk(100, 0.07, trt[tre$tip.label == 't3'], 0.005)
set.seed(1)
x2 <- rwalk(50, x1[51], trt[tre$tip.label == 't4'], 0.005)


pdf('fig_mod_trt.pdf', width = 3, height = 2)
par(mar = c(2, 2, 0, 0) + 0.5, mgp = c(1, 0, 0))

plot(1:100, x1, ylim = range(x1, x2), type = 'n', axes = FALSE, 
     xlab = 'Millions of years', ylab = 'Trait values')
box()
axisArrows(1, length = 0.1, lwd = 2)

segments(x0 = 1:99, x1 = 2:100, y0 = x1[-100], y1 = x1[-1], 
         col = quantCol(x1[-101], pal = viridis(10), xlim = range(trt)), 
         lwd = 2)

segments(x0 = 51:99, x1 = 52:100, y0 = x2[-50], y1 = x2[-1], 
         col = quantCol(x2[-50], pal = viridis(10), xlim = range(trt)), 
         lwd = 2)

dev.off()




## demography through time
set.seed(123)
y1 <- cumsum(c(10, rnorm(100)))

pdf('fig_mod_pop.pdf', width = 3.5, height = 1.8)
par(mar = c(2, 2, 0, 0) + 0.5, mgp = c(1, 0, 0))
plot(cumsum(c(0, rexp(y1[-1]))), y1, type = 's', lwd = 2, axes = FALSE,
     xlab = 'Hundreds of years', ylab = '', col = quantCol(trt['t3'], viridis(5), xlim = range(trt)))
mtext('Abundance', side = 2, line = 0.25)
box()
axisArrows(1, length = 0.1, lwd = 2)
dev.off()


## coalesent
set.seed(12)
coal <- rphylo(10, 1, 0.7)

pdf('fig_mod_coal.pdf', width = 3.5, height = 1.8)
par(mar = rep(0.1, 4))
plot(coal, show.tip.label = FALSE, edge.width = 3,
     edge.col = quantCol(trt['t3'], viridis(5), xlim = range(trt)))
dev.off()


## abundances

pdf('fig_mod_abund.pdf', width = 1.5, height = 3.5)
par(mfcol = c(2, 1), mar = c(3, 1.3, 0, 0) + 0.2, mgp = c(1.5, 1, 0))

comm1 <- table(abund1)
barplot(sort(comm1, decreasing = TRUE), col = colTrt[c(2, 1, 4, 3)], 
        ylab = '', yaxt = 'n', names.arg = NA)
# axisArrows(2, length = 0.1)
box()
mtext('Abundance', side = 2, line = 0.25)

comm2 <- table(abund2)
barplot(sort(comm2, decreasing = TRUE), col = colTrt[as.numeric(names(sort(comm2, TRUE)))], 
        ylab = '', yaxt = 'n', names.arg = NA)
# axisArrows(2, length = 0.1)
box()
mtext('Abundance', side = 2, line = 0.25)

dev.off()




## pop gen

pdf('fig_mod_genDiv.pdf', width = 1.5, height = 3.5)
par(mfcol = c(2, 1), mar = c(3, 1.3, 0, 0) + 0.2, mgp = c(1.5, 1, 0))

comm1 <- table(abund1)
barplot(c(4, 2, 1, 0.25), col = colTrt[c(2, 1, 4, 3)], 
        ylab = '', yaxt = 'n', names.arg = NA)
# axisArrows(2, length = 0.1)
box()
mtext(expression('Genetic div'~(pi)), side = 2, line = 0.25)

comm2 <- table(abund2)
barplot(c(1, 2, 1.7, 3), col = colTrt[as.numeric(names(sort(comm2, TRUE)))], 
        ylab = '', yaxt = 'n', names.arg = NA)
# axisArrows(2, length = 0.1)
box()
mtext(expression('Genetic div'~(pi)), side = 2, line = 0.25)

dev.off()


# 
# col1 <- colAlpha(quantCol(0.9 * max(trt), viridis(10), xlim = range(trt)), 0.2)
# col2 <- colAlpha(quantCol(0.5 * min(trt), viridis(10), xlim = range(trt)), 0.2)
# 
# pdf('fig_mod_abund.pdf', width = 3, height = 6)
# par(mfcol = c(2, 1), oma = c(3, 3, 0, 0) + 0.1, mar = rep(0.2, 4), mgp = c(1, 1, 0))
# 
# set.seed(1)
# plot(sort(rtnegb(100, 0.1, mu = 6), TRUE), log = 'y', 
#      axes = FALSE, frame.plot = TRUE, 
#      panel.first = rect(par('usr')[1], 10^par('usr')[3], par('usr')[2], 10^par('usr')[4], 
#                         col = col1, border = NA))
# 
# set.seed(12)
# plot(sort(rtnegb(100, 4, mu = 6), TRUE), log = 'y', 
#      axes = FALSE, frame.plot = TRUE, 
#      panel.first = rect(par('usr')[1], 10^par('usr')[3], par('usr')[2], 10^par('usr')[4], 
#                         col = col2, border = NA))
# 
# mtext('Species rank', side = 1, line = 0.5, outer = TRUE, cex = 1.5)
# mtext('Abundance', side = 2, line = 0.5, outer = TRUE, cex = 1.5)
# 
# dev.off()
# 
# 
# ## pi
# 
# pdf('fig_mod_genDiv.pdf', width = 3, height = 6)
# par(mfcol = c(2, 1), oma = c(3, 3, 0, 0) + 0.1, mar = rep(0.2, 4), mgp = c(1, 1, 0))
# 
# set.seed(1)
# plot(sort(rnbinom(100, 1, mu = 20), TRUE), 
#      axes = FALSE, frame.plot = TRUE, 
#      panel.first = rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], 
#                         col = col1, border = NA))
# 
# set.seed(12)
# plot(sort(rnbinom(100, 4, mu = 20), TRUE), 
#      axes = FALSE, frame.plot = TRUE, 
#      panel.first = rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], 
#                         col = col2, border = NA))
# 
# mtext('Population rank', side = 1, line = 0.5, outer = TRUE, cex = 1.5)
# mtext(expression('Genetic diversity'~(pi)), side = 2, line = 0.5, outer = TRUE, cex = 1.5)
# 
# dev.off()
# 



# ## complex plotting
# 
# layout(matrix(c(1, 1, 7, 7, 1, 1, 7, 7, 2, 2, 8, 9, 2, 2, 8, 9, 
#                 3, 4, 0, 0, 5, 6, 0, 0), 
#               nrow = 4))
# 
# ### tree
# par(mar = rep(0.1, 4))
# plot(tre, show.tip.label = FALSE)
# tiplabels(pch = 21, bg = quantCol(trt, pal = viridis(10)))
# 
# ### local
# par(mar = rep(0.1, 4))
# plot(1, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, type = 'n', axes = FALSE)
# locPlot(0.75, -0.5, 1.5, 20, trt, 0.9 * max(trt), cex = 1)
# locPlot(0.75, 0.5, -1.5, 20, trt, 0.5 * min(trt), cex = 1)
# 
# ### abundances
# 
# plot(sort(rtnegb(100, 0.1, mu = 6), TRUE), log = 'y')
# plot(sort(rtnegb(100, 1, mu = 6), TRUE), log = 'y')
# 
# ### pi
# plot(sort(rnbinom(100, 1, mu = 20), TRUE))
# plot(sort(rnbinom(100, 4, mu = 20), TRUE))
# 
# 
# ### traits
# 
# par(mar = c(2, 2, 3, 0) + 0.5, mgp = c(1, 0, 0))
# plot(1:100, x1, ylim = range(x1, x2), type = 'n', axes = FALSE, 
#      xlab = 'Millions of years', ylab = 'Trait values')
# box()
# axisArrows(1, length = 0.1, lwd = 2)
# 
# segments(x0 = 1:99, x1 = 2:100, y0 = x1[-100], y1 = x1[-1], 
#          col = quantCol(x1[-101], pal = viridis(10), xlim = range(trt)), 
#          lwd = 2)
# 
# segments(x0 = 51:99, x1 = 52:100, y0 = x2[-50], y1 = x2[-1], 
#          col = quantCol(x2[-50], pal = viridis(10), xlim = range(trt)), 
#          lwd = 2)