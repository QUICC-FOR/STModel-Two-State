#!/usr/bin/env Rscript
col.color = '#377eb8'
ext.color = '#e41a1c'
cex.axis=0.6
lwd=1.5
lab.line=1.5
cex.annote=0.7
ytop.adj=-0.02
cex.id=0.9
fut.xadj = 1.4
fut.yadj = -0.005
arr.yadj = 0.07
arr.xsadj = -0.3
arr.xeadj = -0.1
arrtxt.xadj = -0.4
arrtxt.yadj = -0.04

rc = readRDS("res/resp_curve/28731-ACE-SAC_respCurve.rds")
clim.map = readRDS('dat/climateGrid_unscaled.rds')[,c('lat', 'lon', 'annual_mean_temp')]

# fit a transformation from temperature to latitude. turns out it's pretty linear
library(gam)
latmod = gam(lat ~ s(annual_mean_temp,3), data=clim.map)
warming = 2.5
T = rc$temp
lat = predict(latmod, newdata=list(annual_mean_temp = T))
lat.warm = predict(latmod, newdata=list(annual_mean_temp = T-warming))

C = rc$col.temp
E = rc$ext.temp
xlims = c(35,51)
ylims = c(0,0.25)
rts = predict(latmod, newdata=list(annual_mean_temp=c(2.15, 12.13402)))
label.col = function() mtext(expression(bold(Colonization)), side=2, line=lab.line, col=col.color)
label.ext = function() mtext(expression(bold(Extinction)), side=4, line=lab.line, col=ext.color)

figure.width = 6.5
figure.height = 3.5
filename = file.path('img', 'figs', 'fig1.pdf')
pdf(width=figure.width, height=figure.height, file=filename)

layout(matrix(c(1,2), nrow=2))
par(mgp=c(1,.3,0), tcl=-0.25, mar=c(0.5,2.5,0.5,2.5), oma=c(1.3,0,0,0), cex.axis=cex.axis)
plot(0,0, type='n', ylim=ylims, ylab="", bty='u', xlim=xlims, xlab="", xaxt='n', yaxt='n')
axis(2, cex.axis=cex.axis)
axis(4, cex.axis=cex.axis)
polygon(rep(rts, each=2), c(ylim, rev(ylim)), col="#77777744", border=NA)
lines(lat, C, col=col.color, lwd=lwd)
lines(lat, E, col=ext.color, lwd=lwd)
label.col()
label.ext()
text(mean(rts), ylims[2]+ytop.adj, "Present range\n(at equilibrium)", cex=cex.annote)
text(min(xlims), ylims[2]+ytop.adj/2, "A.", cex=cex.id)

plot(0,0, type='n', ylim=ylims, ylab="", bty='u', xlim=xlims, xlab="", yaxt='n')
axis(2, cex.axis=cex.axis)
axis(4, cex.axis=cex.axis)
ext.x = rep(c(rts[2], rts[2] + warming), each=2)
pres.x = rep(c(rts[2] + warming, rts[1]), each=2)
col.x = rep(c(rts[1], rts[1]+warming), each=2)
polygon(ext.x, c(ylim, rev(ylim)), col=paste0(ext.color, "22"), border=NA)		
polygon(pres.x, c(ylim, rev(ylim)), col="#55555544", border="#55555544")
polygon(col.x, c(ylim, rev(ylim)), col=paste0(col.color, "22"), border=NA)
lines(lat.warm, C, col=col.color, lwd=lwd)
lines(lat.warm, E, col=ext.color, lwd=lwd)
mtext("Latitude", side=1, cex=0.9, line=0.8)
label.col()
label.ext()
text(mean(pres.x), ylims[2]+ytop.adj, "Continued\npresence", cex=cex.annote)
text(mean(ext.x)+fut.xadj, ylims[2]+ytop.adj+fut.yadj, "Future\nextinction", pos=2, cex=cex.annote)
text(mean(col.x)-fut.xadj, ylims[2]+ytop.adj+fut.yadj, "Future\ncolonization", pos=4, cex=cex.annote)
text(min(xlims)+arrtxt.xadj, mean(ylims)+arrtxt.yadj, "Warming moves the\nniche poleward", pos=4, cex=cex.annote)
arrows(min(xlims) + arr.xsadj, mean(ylims)-arr.yadj, rts[2]+arr.xeadj, mean(ylims)-arr.yadj , lwd=2, length=0.1)
text(min(xlims), ylims[2]+ytop.adj/2, "B.", cex=cex.id)
dev.off()
