library(ape)
library(viridis)
library(socorro)
library(pika)

## regional phylo
set.seed(5)
tre <- rphylo(20, 0.9, 0.8)
set.seed(2)
trt <- rTraitCont(tre)
trt[tre$tip.label == 't2'] <- 2 * trt[tre$tip.label == 't2']

pdf('fig_mod_phylo.pdf', width = 3, height = 3)
par(mar = rep(0.1, 4))
plot(tre, show.tip.label = FALSE)
tiplabels(pch = 21, bg = quantCol(trt, pal = viridis(10)))
dev.off()

## local comm
locPlot <- function(r, x0, y0, n0, trt, env, ...) {
    ## xy for dots
    coord <- seq(-r, r, length.out = ceiling(sqrt(n0 / (pi*0.5^2))))
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
    
    polygon(x0 + rOut * cos(seq(0, 2*pi, length.out = 100)), 
            y0 + rOut * sin(seq(0, 2*pi, length.out = 100)), 
            col = colAlpha(quantCol(env, viridis(10), xlim = range(trt)), 0.2))
    
    points(xy, bg = quantCol(n1, viridis(10), xlim = range(trt)), 
           pch = 21, ...)
    
}

pdf('fig_mod_comm.pdf', width = 3, height = 3)
par(mar = rep(0.1, 4))
plot(1, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, type = 'n', axes = FALSE)
locPlot(0.75, -0.5, 1, 20, trt, 0.9 * max(trt), cex = 1)
locPlot(0.75, 0.5, -1, 20, trt, 0.5 * min(trt), cex = 1)
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

set.seed(123)
x1 <- rwalk(100, -0.1, trt[tre$tip.label == 't11'], 0.01)
x2 <- rwalk(50, x1[51], trt[tre$tip.label == 't2'], 0.01)


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

pdf('fig_mod_pop.pdf', width = 4, height = 2)
par(mar = c(2, 2, 0, 0) + 0.5, mgp = c(1, 0, 0))
plot(cumsum(c(0, rexp(y1[-1]))), y1, type = 's', lwd = 2, axes = FALSE,
     xlab = 'Hundreds of years', ylab = 'Abundance')
box()
axisArrows(1, length = 0.1, lwd = 2)
dev.off()


## coalesent
set.seed(12)
coal <- rphylo(10, 1, 0.7)

pdf('fig_mod_coal.pdf', width = 4, height = 2)
par(mar = rep(0.1, 4))
plot(coal, show.tip.label = FALSE)
dev.off()

foo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
s <- diff(c(foo$xx[-(1:Ntip(coal))], max(foo$xx)))
m <- s * 2:Ntip(coal)
mut <- sample(Ntip(coal) - 1, 12, prob = m, replace = TRUE)


## abundances

col1 <- colAlpha(quantCol(0.9 * max(trt), viridis(10), xlim = range(trt)), 0.2)
col2 <- colAlpha(quantCol(0.5 * min(trt), viridis(10), xlim = range(trt)), 0.2)

pdf('fig_mod_abund.pdf', width = 3, height = 6)
par(mfcol = c(2, 1), oma = c(3, 3, 0, 0) + 0.1, mar = rep(0.2, 4), mgp = c(1, 1, 0))

set.seed(1)
plot(sort(rtnegb(100, 0.1, mu = 6), TRUE), log = 'y', 
     axes = FALSE, frame.plot = TRUE, 
     panel.first = rect(par('usr')[1], 10^par('usr')[3], par('usr')[2], 10^par('usr')[4], 
                        col = col1, border = NA))

set.seed(12)
plot(sort(rtnegb(100, 4, mu = 6), TRUE), log = 'y', 
     axes = FALSE, frame.plot = TRUE, 
     panel.first = rect(par('usr')[1], 10^par('usr')[3], par('usr')[2], 10^par('usr')[4], 
                        col = col2, border = NA))

mtext('Species rank', side = 1, line = 0.5, outer = TRUE, cex = 1.5)
mtext('Abundance', side = 2, line = 0.5, outer = TRUE, cex = 1.5)

dev.off()


## pi

pdf('fig_mod_genDiv.pdf', width = 3, height = 6)
par(mfcol = c(2, 1), oma = c(3, 3, 0, 0) + 0.1, mar = rep(0.2, 4), mgp = c(1, 1, 0))

set.seed(1)
plot(sort(rnbinom(100, 1, mu = 20), TRUE), 
     axes = FALSE, frame.plot = TRUE, 
     panel.first = rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], 
                        col = col1, border = NA))

set.seed(12)
plot(sort(rnbinom(100, 4, mu = 20), TRUE), 
     axes = FALSE, frame.plot = TRUE, 
     panel.first = rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], 
                        col = col2, border = NA))

mtext('Population rank', side = 1, line = 0.5, outer = TRUE, cex = 1.5)
mtext(expression('Genetic diversity'~(pi)), side = 2, line = 0.5, outer = TRUE, cex = 1.5)

dev.off()




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