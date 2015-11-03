#!/usr/bin/env Rscript
col.color = '#377eb8'
ext.color = '#e41a1c'
lwd=1.5
lab.line=0.5
ylim=c(0,0.7)
cex.annote=0.7
ytop.adj=-0.07
cex.id=0.9

T = seq(-0.75, 0.75, length.out=500)
col = function(x) plogis(-10*x^2)
ext = function(x) plogis(-1.3+3*x^2)
lam = function(x) col(x) - ext(x)
rts = c(uniroot(lam, interval=c(-1,0))$root, uniroot(lam, interval=c(0,1))$root)
label.col = function() mtext(expression(bold(Colonization)), side=2, line=0.5, col=col.color)
label.ext = function() mtext(expression(bold(Extinction)), side=4, line=0.5, col=ext.color)

figure.width = 6.5
figure.height = 3.5
filename = file.path('img', 'figs', 'fig1.pdf')
pdf(width=figure.width, height=figure.height, file=filename)

layout(matrix(c(1,2), nrow=2))
par(mgp=c(1,.5,0), mar=c(0.5,2,0.5,2), oma=c(1.3,0,0,0))
plot(0,0, type='n', ylim=ylim, ylab="", bty='u', xlim=c(-0.75, 0.75), xlab="", xaxt='n', yaxt='n')
polygon(rep(rts, each=2), c(ylim, rev(ylim)), col="#77777744", border=NA)
lines(T, col(T), col=col.color, lwd=lwd)
lines(T, ext(T), col=ext.color, lwd=lwd)
label.col()
label.ext()
text(mean(rts), ylim[2]+ytop.adj, "Present range\n(at equilibrium)", cex=cex.annote)
text(min(T), ylim[2]+ytop.adj/2, "A.", cex=cex.id)

warming = 0.2
plot(0,0, type='n', ylim=ylim, ylab="", bty='u', xlim=range(T), xlab="", xaxt='n', yaxt='n')
ext.x = rep(c(rts[1], rts[1] + warming), each=2)
pres.x = rep(c(rts[1] + warming, rts[2]), each=2)
col.x = rep(c(rts[2], rts[2]+warming), each=2)
polygon(ext.x, c(ylim, rev(ylim)), col=paste0(ext.color, "22"), border=NA)		
polygon(pres.x, c(ylim, rev(ylim)), col="#55555544", border="#55555544")
polygon(col.x, c(ylim, rev(ylim)), col=paste0(col.color, "22"), border=NA)
lines(T, col(T-warming), col=col.color, lwd=lwd)
lines(T, ext(T-warming), col=ext.color, lwd=lwd)
mtext("Latitude", side=1, cex=0.9)
label.col()
label.ext()
text(mean(pres.x), ylim[2]+ytop.adj, "Continued\npresence", cex=cex.annote)
text(mean(ext.x)+0.1, ylim[2]+ytop.adj-0.017, "Future\nextinction", pos=2, cex=cex.annote)
text(mean(col.x)-0.1, ylim[2]+ytop.adj-0.017, "Future\ncolonization", pos=4, cex=cex.annote)
text(min(T)-0.02, mean(ylim)-0.1, "Warming moves the\nniche poleward", pos=4, cex=cex.annote)
arrows(min(T), mean(ylim), rts[1]-0.1, mean(ylim), lwd=2, length=0.1)
text(min(T), ylim[2]+ytop.adj/2, "B.", cex=cex.id)
dev.off()
