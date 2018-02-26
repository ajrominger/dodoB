pdf('fig_hyp.pdf', width = (3 - 0.8) * 3, height = (3 - 0.42) * 2)

split.screen(matrix(c(0, 1/3, 0.5, 1, 
                      1/3, 2/3, 0.5, 1, 
                      2/3, 1, 0.5, 1, 
                      0, 2/3, 0, 0.5, 
                      2/3, 1, 0, 0.5), 
                    ncol = 4, byrow = TRUE))

screen(1)
par(mar = c(2, 2, 2, 0.1), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, frame.plot = TRUE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = expression('Immigration'~(nu)), ylab = expression('Speciation'~(lambda)))
segments(1, 2, 2, 1, lwd = 3)

par(xpd = NA)
rect(par('usr')[1], par('usr')[4], par('usr')[2], par('usr')[4] + 0.3, col = 'gray')
mtext('H1', side = 3, line = 0.2)

screen(2)
par(mar = c(2, 2, 2, 0.1), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, frame.plot = TRUE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = 'Enviro. var.', ylab = expression('Speciation'~(lambda)))
segments(1, 1, 2, 2, lwd = 3)

par(xpd = NA)
rect(par('usr')[1], par('usr')[4], par('usr')[2], par('usr')[4] + 0.3, col = 'gray')
mtext('H2', side = 3, line = 0.2)

screen(3)
par(mar = c(2, 2, 2, 0.1), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, frame.plot = TRUE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = expression('Competition'~(alpha[ij])), ylab = 'Diversity')
segments(1, 1, 2, 2, lwd = 3)

par(xpd = NA)
rect(par('usr')[1], par('usr')[4], par('usr')[2], par('usr')[4] + 0.3, col = 'gray')
mtext('H3', side = 3, line = 0.2)


split.screen(1:2, 4, erase = FALSE)

screen(6, new = FALSE)
par(mar = c(2, 2, 2, 0), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, frame.plot = TRUE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = '', ylab = expression('Enviro. filtering'~(rho)))
segments(1, 2, 2, 1, lwd = 3)

screen(7, new = FALSE)
par(mar = c(2, 2, 2, 0.1), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, frame.plot = TRUE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = '', ylab = expression('Competition'~(alpha[ij])))
segments(1, 2, 2, 1, lwd = 3)


screen(4, new = FALSE)
par(mar = c(2, 2, 2, 0.1), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = expression('Immigration'~(nu)), ylab = '')

par(xpd = NA)
rect(par('usr')[1], par('usr')[4], par('usr')[2], par('usr')[4] + 0.3, col = 'gray')
mtext('H4', side = 3, line = 0.2)

screen(5)
par(mar = c(2, 2, 2, 0.1), mgp = c(0.25, 0, 0))
plot(1:2, type = 'n', axes = FALSE, frame.plot = TRUE, xlim = c(0.75, 2.25), ylim = c(0.75, 2.25), 
     xlab = 'Equilib. diversity', ylab = 'Obs. diversity')
abline(0, 1, lty = 2, col = 'gray', lwd = 2)

set.seed(2)
x <- seq(1, 1.45, length.out = 8) + runif(8, -0.05, 0.05)
y <- x - runif(8, 0.1, 0.35)
points(x, y, col = hsv(0.6, 0.6, 0.4), pch = 16)

set.seed(4)
x <- seq(1.55, 2, length.out = 8) + runif(8, -0.05, 0.05)
y <- x + runif(8, 0.1, 0.35)
points(x, y, col = hsv(0.1, 1, 0.8), pch = 16)

text(1.35, 1.4, labels = '1:1 line', srt = 45, col = 'gray50', cex = 0.8, adj = c(0.5, 0))
text(1.25, 0.8, labels = 'Immigration', col = hsv(0.6, 0.6, 0.4), adj = c(0, 0.5), cex = 0.9)
text(1.75, 2.2, labels = 'Competition', col = hsv(0.1, 1, 0.8), adj = c(1, 0.5), cex = 0.9)

par(xpd = NA)
rect(par('usr')[1], par('usr')[4], par('usr')[2], par('usr')[4] + 0.3, col = 'gray')
mtext('H5', side = 3, line = 0.2)

close.screen(all.screens = TRUE)

dev.off()
