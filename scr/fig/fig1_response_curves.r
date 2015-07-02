








varNames = readRDS("dat/climVariableNames.rds")
e1.lab = varNames[which(varNames[,1] == env1),2]
e2.lab = varNames[which(varNames[,1] == env2),2]

pdf(w=6.5, h=8, file="img/response_curves.pdf")
par(mfrow=c(4,2), bty='n', mar=c(3,3,0.5, 0.5), mgp=c(2,0.5,0), tck=-0.03)
for(i in 1:4)
{
#par(fig=c(0, 1, 0, 1), bty='o', mar=c(0.5, 0.5, 0.5, 0.5))
#plot(0,0, type='n',xaxt='n', yaxt='n', xlab='', ylab='')
#par(fig=c(0.04, 0.5, 0.04, 0.96), bty='n', new=TRUE, mar=c(4,4, 0.1, 0.1))
plot(x1.us, env1.yc[,1], xlab=e1.lab, ylab="Probability", ylim=c(0,0.25), col='blue', type='l')
polygon(c(x1.us, rev(x1.us)), c(env1.yc[,2], rev(env1.yc[,3])), col="#0000FF22", border=NA)
lines(x1.us, env1.ye[,1], col='red')
polygon(c(x1.us, rev(x1.us)), c(env1.ye[,2], rev(env1.ye[,3])), col="#FF000022", border=NA)

#par(fig=c(0.5, 0.96, 0.04, 0.96), bty='n', new=TRUE)
plot(x2.us, env2.yc[,1], xlab=e2.lab, ylab="", ylim=c(0,0.25), col='blue', type='l')
polygon(c(x2.us, rev(x2.us)), c(env2.yc[,2], rev(env2.yc[,3])), col="#0000FF22", border=NA)
lines(x2.us, env2.ye[,1], col='red')
polygon(c(x2.us, rev(x2.us)), c(env2.ye[,2], rev(env2.ye[,3])), col="#FF000022", border=NA)
}

dev.off()