#!/usr/bin/Rscript
varNames = readRDS("dat/climVariableNames.rds")
spNames = c('18032-ABI-BAL', '28728-ACE-RUB', '28731-ACE-SAC', '19489-BET-PAP',
		'19462-FAG-GRA', '32931-FRA-AME', '32929-FRA-PEN', '183295-PIC-GLA', 
		'183319-PIN-BAN', '183385-PIN-STR', '195773-POP-TRE', '19290-QUE-ALB', 
		'19447-QUE-VEL', '19049-ULM-AME')
spInfoAll = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE, colClasses='character')

pdf(w=6.5, h=8, file="img/response_curves.pdf")
par(mfrow=c(4,2), bty='n', mar=c(3,3,0.5,0.5), oma=c(0,0,0,2), mgp=c(2,0.5,0), tck=-0.03)
for(spName in spNames)
{


	spInfo = spInfoAll[spInfoAll$spName == spName,]
	e1 = spInfo$env1
	e2 = spInfo$env2
	ylim1 = as.numeric(spInfo$rylim1)
	ylim2 = as.numeric(spInfo$rylim2)
	spLab = bquote(italic(.(spInfo$genus)~.(spInfo$species)))

	e1.lab = varNames[which(varNames[,1] == e1),2]
	e2.lab = varNames[which(varNames[,1] == e2),2]

	resp = readRDS(file.path('species', spName, 'res', paste(spName, 'responseCurves.rds', sep='_')))
	

	with(resp,
	{
		plot(env1, col1.mean, xlab=e1.lab, ylab="Probability", ylim=c(0,ylim1), col='blue', type='l')
		polygon(c(env1, rev(env1)), c(col1.lo, rev(col1.up)), col="#0000FF22", border=NA)
		lines(env1, ext1.mean, col='red')
		polygon(c(env1, rev(env1)), c(ext1.lo, rev(ext1.up)), col="#FF000022", border=NA)

		plot(env2, col2.mean, xlab=e2.lab, ylab="", ylim=c(0,ylim2), col='blue', type='l')
		polygon(c(env2, rev(env2)), c(col2.lo, rev(col2.up)), col="#0000FF22", border=NA)
		lines(env2, ext2.mean, col='red')
		polygon(c(env2, rev(env2)), c(ext2.lo, rev(ext2.up)), col="#FF000022", border=NA)
		xpd = par()$xpd
		par(xpd=NA)
		mtext(spLab, side=4, line=0.5, cex=0.8)
		par(xpd=xpd)
	})
}

dev.off()