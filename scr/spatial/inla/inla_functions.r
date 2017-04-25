ext_fit <- function(dat.ext)
{
	coo.ext <- as.matrix(dat.ext[,c('x', 'y')])
	bound.ext <- inla.nonconvex.hull(coo.ext)
	mesh.ext <- inla.mesh.2d(boundary = bound.ext, max.edge=c(75000,150000), cutoff=25000)
	A.ext <- inla.spde.make.A(mesh=mesh.ext, loc=coo.ext)
	
	# extract coordinates and set up spde model
	spde.ext <- inla.spde2.matern(mesh=mesh.ext, alpha=1.5)
	
	# setup inla stack
	stk.ext <- inla.stack(tag='extinction', data=list(y=dat.ext$ext), A=list(A.ext, 1),
												effects=list(s=1:spde.ext$n.spde, list(intercept=rep(1, nrow(dat.ext)), 
																															 t=dat.ext$annual_mean_temp, p=dat.ext$tot_annual_pp)))
	
	# build 2 models:
	#	climate effects only
	#	climate effects + space
	
	formulas <- list(
		f.ext.c = y ~ 0 + intercept + t + p + I(t^2) + I(p^2),
		f.ext.cs = y ~ 0 + intercept + t + p + I(t^2) + I(p^2) + f(s, model=spde.ext))
	
	models.e <- lapply(formulas, function(f) {
		inla(f, data=inla.stack.data(stk.ext), control.predictor=list(A=inla.stack.A(stk.ext)), 
				 family = "binomial", control.family = list(link = "logit"), control.compute=list(config=TRUE))
	})
	
	list(models = models.e, data = stk.ext, formulas = formulas, mesh = mesh.ext)
}


do_inla_sample <- function(mod, nms)
{
	samples <- inla.posterior.sample(n=1000, result=mod, add.names=FALSE)
	xnames <- rownames(samples[[1]]$latent)
	idx <- lapply(nms, function(n) {
		if(n == 's') {
			which(substr(xnames, 1, nchar(n)) == n)
		} else {
			which(xnames == n)
		}
	})
	
	mat.samples <- sapply(samples, function(spl)
	{
		if('s' %in% nms) {
			c(m=spl$latent[idx[[1]]], t=spl$latent[idx[[2]]], p=spl$latent[idx[[3]]], t2=spl$latent[idx[[4]]], 
				p2=spl$latent[idx[[5]]], s=spl$latent[idx[[6]]])
		} else {
			c(m=spl$latent[idx[[1]]], t=spl$latent[idx[[2]]], p=spl$latent[idx[[3]]], t2=spl$latent[idx[[4]]], 
				p2=spl$latent[idx[[5]]])		
		}
	})
	mat.samples
}



get_inla_samples <- function(mod, new.coords)
{
	A.pr <- inla.spde.make.A(mesh=mod$mesh, loc=new.coords)
	
	## draw samples from the model and rearrange into a matrix for prediction
	pr.names <- c("intercept", "t", "p", "I(t^2)", "I(p^2)")
	predict.names <- list(pr.names, c(pr.names, 's'))
	samples <- mapply(do_inla_sample, mod$models, predict.names, SIMPLIFY=FALSE)
	list(samples=samples, A.mat=A.pr)
}


inla_lp <- function(inla_samples, mcmc_samples, a_mat, temp, precip)
{
	
	predGrid <- cbind(m=1, t=temp, p=precip, t2=temp^2, p2=precip^2)
	if(is.null(a_mat))
	{
		spatial_pr <- as.matrix(predGrid[,1:5] %*% inla_samples[[2]][1:5,])
	} else 
	{
		predGrid <- cbind(predGrid, s=a_mat)
		spatial_pr <- as.matrix(predGrid %*% inla_samples[[2]])
	}
	mcmc_pr <- as.matrix(predGrid[,1:5] %*% t(as.matrix(mcmc_samples)))
	nonspatial_pr <- as.matrix(predGrid[,1:5] %*% inla_samples[[1]])
	list(mcmc=mcmc_pr, nonspatial=nonspatial_pr, spatial=spatial_pr)
}



grid_summary <- function(lst, coords)
{
	res <- data.frame(x = coords[,'x'], y = coords[,'y'])
	res <- cbind(res, as.data.frame(sapply(lst, rowMeans)))
	res <- cbind(res, as.data.frame(sapply(lst, function(x) apply(x, 1, quantile, 0.025))))
	res <- cbind(res, as.data.frame(sapply(lst, function(x) apply(x, 1, quantile, 0.025))))
	colnames(res) <- as.vector(outer(names(lst), c('mean', 'lower', 'upper'), paste, sep='.'))
	res
}



inla_random_field <- function(mesh, mod, xlim, ylim, dims=c(200,200))
{
	gproj <- inla.mesh.projector(mesh, xlim=xlim, ylim=ylim, dims=dims)
	meanfield <- inla.mesh.project(gproj, mod$summary.random$s$mean)
	sdfield <- inla.mesh.project(gproj, mod$summary.random$s$sd)
	xx = seq(xlim[1], xlim[2], length.out=dims[1])
	yy = seq(ylim[1], ylim[2], length.out=dims[2])
	list(mean=raster(list(x=xx, y=yy, z=meanfield)), sd=raster(list(x=xx, y=yy, z=sdfield)))
}



choose_range <- function(mi, ma)
{
	if((mi <= 0 & ma <= 0) | (mi >= 0 & ma >= 0))
		return(c(mi, ma))
	
	if(abs(mi) > ma) {
		return(c(mi, abs(mi)))
	} else {
		return(c(-1 * ma, ma))
	}
}


inla_rf_plots <- function(projs, basedata, col.mean, col.sd, ...)
{
	if(missing(col.mean))
		col.mean <- colorRampPalette(c('#b2182b', '#ef8a62', '#fddbc7', '#ffffff', '#d1e5f0', '#67a9cf', '#2166ac'))(100)
	if(missing(col.sd))
		col.sd <- colorRampPalette(c('#f2f0f7', '#cbc9e2', '#9e9ac8', '#756bb1', '#54278f'))(100)
	
	miv <- minValue(projs$mean)
	mav <- maxValue(projs$mean)
	zl <- choose_range(miv, mav)
	
	plot(projs$mean, col=col.mean, xaxt='n', yaxt='n', zlim=zl, ...)
	if(!missing(basedata))
		plot(basedata$ocean, add=TRUE)
	plot(projs$sd, col=col.sd, xaxt='n', yaxt='n')
	if(!missing(basedata))
		plot(basedata$ocean, add=TRUE)
}


inla_rc_lims <- function(spName)
{
	xlims <- list("28728-ACE-RUB" = c(-5,20),
								"28731-ACE-SAC" = c(0,15),
								"19462-FAG-GRA" = c(0,20),
								"32931-FRA-AME" = c(0,20),
								"32945-FRA-NIG" = c(0,20),
								"19287-QUE-MAC" = c(0,15),
								"19408-QUE-RUB" = c(0,20),
								"183397-TSU-CAN" = c(0,15),
								"19481-BET-ALL" = c(-5,15),
								"183375-PIN-RES" = c(-5,15),
								"183385-PIN-STR" = c(-5,15),
								"22463-POP-GRA" = c(-5,15),
								"18032-ABI-BAL" = c(-5, 10),
								"19489-BET-PAP" = c(-5, 10),
								"183412-LAR-LAR" = c(-5, 10),
								"183295-PIC-GLA" = c(-5, 10),
								"183302-PIC-MAR" = c(-5, 10),
								"18034-PIC-RUB" = c(-5, 10),
								"183319-PIN-BAN" = c(-5, 10),
								"195773-POP-TRE" = c(-5, 10),
								"505490-THU-OCC" = c(-5, 10))
	ylims <- list("28728-ACE-RUB" = c(0, 0.4),
								"28731-ACE-SAC" = c(0, 0.5),
								"19462-FAG-GRA" = c(0, 0.4),
								"32931-FRA-AME" = c(0, 0.4),
								"32945-FRA-NIG" = c(0, 0.5),
								"19287-QUE-MAC" = c(0, 0.8),
								"19408-QUE-RUB" = c(0, 0.5),
								"183397-TSU-CAN" = c(0, 0.4),
								"19481-BET-ALL" = c(0, 0.6),
								"183375-PIN-RES" = c(0, 1),
								"183385-PIN-STR" = c(0, 1),
								"22463-POP-GRA" = c(0, 0.8),
								"18032-ABI-BAL" = c(0, 1),
								"19489-BET-PAP" = c(0, 0.4),
								"183412-LAR-LAR" = c(0, 0.6),
								"183295-PIC-GLA" = c(0, 0.6),
								"183302-PIC-MAR" = c(0, 0.4),
								"18034-PIC-RUB" = c(0, 0.6),
								"183319-PIN-BAN" = c(0, 0.3),
								"195773-POP-TRE" = c(0, 0.4),
								"505490-THU-OCC" = c(0, 0.5))
	xlimsP <- list("28728-ACE-RUB" = c(500,2000),
								 "28731-ACE-SAC" = c(500,2000),
								 "19462-FAG-GRA" = c(500,2000),
								 "32931-FRA-AME" = c(500,2000),
								 "32945-FRA-NIG" = c(500,2000),
								 "19287-QUE-MAC" = c(500,2000),
								 "19408-QUE-RUB" = c(500,2000),
								 "183397-TSU-CAN" = c(500,2000),
								 "19481-BET-ALL" = c(500,2000),
								 "183375-PIN-RES" = c(500,2000),
								 "183385-PIN-STR" = c(500,2000),
								 "22463-POP-GRA" = c(500,2000),
								 "18032-ABI-BAL" = c(500,2000),
								 "19489-BET-PAP" = c(500,2000),
								 "183412-LAR-LAR" = c(500,2000),
								 "183295-PIC-GLA" = c(500,2000),
								 "183302-PIC-MAR" = c(500,2000),
								 "18034-PIC-RUB" = c(500,2000),
								 "183319-PIN-BAN" = c(500,2000),
								 "195773-POP-TRE" = c(500,2000),
								 "505490-THU-OCC" = c(500,2000))
	ylimsP <- list("28728-ACE-RUB" = c(0, 0.4),
								"28731-ACE-SAC" = c(0, 0.4),
								"19462-FAG-GRA" = c(0, 0.8),
								"32931-FRA-AME" = c(0, 0.6),
								"32945-FRA-NIG" = c(0, 0.9),
								"19287-QUE-MAC" = c(0, 0.9),
								"19408-QUE-RUB" = c(0, 0.6),
								"183397-TSU-CAN" = c(0, 0.6),
								"19481-BET-ALL" = c(0, 0.4),
								"183375-PIN-RES" = c(0, 0.4),
								"183385-PIN-STR" = c(0, 0.4),
								"22463-POP-GRA" = c(0, 0.6),
								"18032-ABI-BAL" = c(0, 0.4),
								"19489-BET-PAP" = c(0, 0.4),
								"183412-LAR-LAR" = c(0, 0.4),
								"183295-PIC-GLA" = c(0, 0.8),
								"183302-PIC-MAR" = c(0, 0.4),
								"18034-PIC-RUB" = c(0, 0.5),
								"183319-PIN-BAN" = c(0, 0.2),
								"195773-POP-TRE" = c(0, 0.6),
								"505490-THU-OCC" = c(0, 0.3))
	xl <- if(spName %in% names(xlims)) {xlims[[spName]]} else c(-5, 25)
	xlp <- if(spName %in% names(xlimsP)) {xlimsP[[spName]]} else c(500, 2000)
	yl <- if(spName %in% names(ylims)) {ylims[[spName]]} else c(0, 1)
	ylp <- if(spName %in% names(ylimsP)) {ylimsP[[spName]]} else c(0, 1)
	list(x=xl, xp=xlp, y=yl, yp=ylp)
}

# not used, kept for historicity
# pr_to_rate <- function(pr, t, prev, continuous=TRUE)
# {
# 	if(missing(prev)) prev <- 1
# 	if(continuous)
# 	{
# 		return(-log(1 - pr)/(prev*t))
# 	} else {
# 		return((1 - (1 - pr)^(t/5)) / prev)
# 	}
# }


