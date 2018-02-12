library(ape)
library(viridis)
library(socorro)

## regional phylo
set.seed(1)
tre <- rphylo(20, 0.9, 0.8)
set.seed(1)
trt <- rTraitCont(tre)

plot(tre, show.tip.label = FALSE)
tiplabels(pch = 21, bg = quantCol(trt, pal = viridis(10)))

## local comm
n0 <- 20

coord <- seq(-1, 1, length.out = ceiling(sqrt(n0 / (pi*0.5^2))))
xy <- expand.grid(coord, coord)
xy <- xy[sqrt(xy[, 1]^2 + xy[, 2]^2) <= 1.05, ]

n1 <- sample(trt, size = nrow(xy), 
             prob = dnorm(trt, mean(trt) + 0.25 * diff(range(trt)), 0.1), 
             replace = TRUE)
n2 <- sample(trt, size = nrow(xy), 
             prob = dnorm(trt, mean(trt) - 0.25 * diff(range(trt)), 0.1), 
             replace = TRUE)


plot(xy, asp = 1, bg = quantCol(n1, viridis(10), xlim = range(trt)), 
     pch = 21, cex = 1.5, axes = FALSE, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
polygon(1.1 * cos(seq(0, 2*pi, length.out = 100)), 
        1.1 * sin(seq(0, 2*pi, length.out = 100)))
plot(xy, asp = 1, bg = quantCol(n2, viridis(10), xlim = range(trt)), 
     pch = 21, cex = 1.5, axes = FALSE, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
polygon(1.1 * cos(seq(0, 2*pi, length.out = 100)), 
        1.1 * sin(seq(0, 2*pi, length.out = 100)))


## trait evol through time

set.seed(1)
x1 <- cumsum(c(0, rnorm(100, 0.2, 2)))
set.seed(1)
x2 <- cumsum(c(0, rnorm(100, -0.4, 2)))
set.seed(2)
x3 <- x1[61] + cumsum(c(0, rnorm(40, -0.4, 2)))

plot(0:100, x1, ylim = range(x1, x2, x3), type = 'n')
segments(x0 = 0:99, x1 = 1:100, y0 = x1[-101], y1 = x1[-1], 
         col = quantCol(x1[-101], pal = viridis(10), xlim = range(x1, x2, x3)))
segments(x0 = 0:99, x1 = 1:100, y0 = x2[-101], y1 = x2[-1], 
         col = quantCol(x2[-101], pal = viridis(10), xlim = range(x1, x2, x3)))
segments(x0 = 60:99, x1 = 61:100, y0 = x3[-41], y1 = x3[-1], 
         col = quantCol(x3[-41], pal = viridis(10), xlim = range(x1, x2, x3)))


## demography through time
set.seed(123)
y1 <- cumsum(c(10, rnorm(100)))
plot(cumsum(c(0, rexp(y1[-1]))), y1, type = 's')


## coalesent
set.seed(12)
coal <- rphylo(10, 1, 0.7)
plot(coal)

foo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
s <- diff(c(foo$xx[-(1:Ntip(coal))], max(foo$xx)))
m <- s * 2:Ntip(coal)
mut <- sample(Ntip(coal) - 1, 12, prob = m, replace = TRUE)
