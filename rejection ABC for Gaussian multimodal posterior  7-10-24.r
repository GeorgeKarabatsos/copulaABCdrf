rm(list = ls())
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
library('ks')	# for kde()

parameterNames		=	c('theta1', 'theta2', 'theta3', 'theta4', 'theta5')
simConditionNames	=	c('simResults:')

postEtheta.RejectionABC.simulations.mean	=	matrix(NA, nrow = 1, ncol = 5)
postEtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 1, ncol = 5)
postMEDtheta.RejectionABC.simulations.mean	=	matrix(NA, nrow = 1, ncol = 5)
postMEDtheta.RejectionABC.simulations.sd	=	matrix(NA, nrow = 1, ncol = 5)
postMODEtheta.RejectionABC.simulations.mean	=	matrix(NA, nrow = 1, ncol = 5)
postMODEtheta.RejectionABC.simulations.sd	=	matrix(NA, nrow = 1, ncol = 5)
MLEtheta.RejectionABC.simulations.mean		=	matrix(NA, nrow = 1, ncol = 5)
MLEtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 1, ncol = 5)
postSDtheta.RejectionABC.simulations.mean	=	matrix(NA, nrow = 1, ncol = 5)
postSDtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 1, ncol = 5)
cover95.RejectionABC.simulations.mean		=	matrix(NA, nrow = 1, ncol = 5)
cover50.RejectionABC.simulations.mean		=	matrix(NA, nrow = 1, ncol = 5)
cover95.exact.simulations.mean				=	matrix(NA, nrow = 1, ncol = 5)
cover50.exact.simulations.mean				=	matrix(NA, nrow = 1, ncol = 5)
KS.Dn.RejectionABC.simulations.mean		=	matrix(NA, nrow = 1, ncol = 5)
KS.Dn.RejectionABC.simulations.sd			=	matrix(NA, nrow = 1, ncol = 5)
KS.stat.RejectionABC.simulations.mean		=	matrix(NA, nrow = 1, ncol = 5)
KS.stat.RejectionABC.simulations.sd		=	matrix(NA, nrow = 1, ncol = 5)
MAE.postEtheta.RejectionABC.simulations 	= 	matrix(NA, nrow = 1, ncol = 5)
MSE.postEtheta.RejectionABC.simulations 	= 	matrix(NA, nrow = 1, ncol = 5)
MAE.postMEDtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 1, ncol = 5)
MSE.postMEDtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 1, ncol = 5)
MAE.postMODEtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 1, ncol = 5) 	
MSE.postMODEtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 1, ncol = 5)	
MAE.MLEtheta.RejectionABC.simulations		= 	matrix(NA, nrow = 1, ncol = 5)	
MSE.MLEtheta.RejectionABC.simulations		= 	matrix(NA, nrow = 1, ncol = 5)	

#KS.simulations.RejectionABC					=	matrix(NA, nrow = 1, ncol = 5)
#KS.stat.simulations.RejectionABC			=	matrix(NA, nrow = 1, ncol = 5)

setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Gaussian Multimodal Posterior')

RDataFiles		=	list.files(path = ".", pattern = "\\.RData$")

postEtheta.RejectionABC.simulations 			=	c()
postMEDtheta.RejectionABC.simulations 			=	c()
postMODEtheta.RejectionABC.simulations 			=	c()
MLEtheta.RejectionABC.simulations 				=	c()
postSDtheta.RejectionABC.simulations 			=	c()
cover95.RejectionABC.simulations				=	c()
cover50.RejectionABC.simulations				=	c()
cover95.exact.simulations						=	c()
cover50.exact.simulations						=	c()
KS.Dn.RejectionABC.simulations					=	c()
KS.stat.RejectionABC.simulations				=	c()
L1error.postEtheta.RejectionABC.simulations		=	c()
L2error.postEtheta.RejectionABC.simulations		=	c()
L1error.postMEDtheta.RejectionABC.simulations	=	c()
L2error.postMEDtheta.RejectionABC.simulations	=	c()
L1error.postMODEtheta.RejectionABC.simulations	= 	c()
L2error.postMODEtheta.RejectionABC.simulations	= 	c()
L1error.MLEtheta.RejectionABC.simulations		= 	c()
L2error.MLEtheta.RejectionABC.simulations		= 	c()

for (fileNum in 1 : length(RDataFiles))	{	# fileNum = 1 # for testing
	load(RDataFiles[fileNum])
	#
	# Do the rejection ABC analysis:
	euclideanDistances			=	as.matrix(pdist(sy.table, Y = sx))
	theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .01),]
	# Above, retained 3% of thetas with summary statistics sy nearest to observed dataset summaries sx.
	
	postEtheta.RejectionABC		=	colMeans(theta.table.RejectionABC)
	postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
	post.kde5d.RejectionABC		=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate
	postMODEtheta.RejectionABC	= 	theta.table.RejectionABC[which.max(post.kde5d.RejectionABC),]
	
	priorPDFs.RejectionABC		=	apply(	cbind(	dunif(theta.table.RejectionABC[,1],min=-3,max=3),dunif(theta.table.RejectionABC[,2],min=-4,max=4),
											dunif(theta.table.RejectionABC[,3],min=-3,max=3),dunif(theta.table.RejectionABC[,4],min=-3,max=3),
											dunif(theta.table.RejectionABC[,5],min=-3,max=3)), 1, prod)
	likelihoods.RejectionABC	=	post.kde5d.RejectionABC / priorPDFs.RejectionABC
	MLEtheta.RejectionABC		=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]
	postSDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, sd)
	#
	postEtheta.RejectionABC.simulations 		= 	rbind(postEtheta.RejectionABC.simulations, 		postEtheta.RejectionABC)
	postMEDtheta.RejectionABC.simulations		= 	rbind(postMEDtheta.RejectionABC.simulations, 	postMEDtheta.RejectionABC)
	postMODEtheta.RejectionABC.simulations	= 	rbind(postMODEtheta.RejectionABC.simulations, 	postMODEtheta.RejectionABC)
	MLEtheta.RejectionABC.simulations			=	rbind(MLEtheta.RejectionABC.simulations, 		MLEtheta.RejectionABC)
	postSDtheta.RejectionABC.simulations		= 	rbind(postSDtheta.RejectionABC.simulations, 	postSDtheta.RejectionABC)

	# Do weighted one-sample (two-tailed) Kolmogorov Smirnov tests (Monahan,2011,p.358)
	KS.Dn.RejectionABC.simulations0		=	c()
	KS.stat.RejectionABC.simulations0	=	c()
	cover95.RejectionABC.simulations0	=	c()
	cover50.RejectionABC.simulations0	=	c()
	cover95.exact.simulations0			=	c()
	cover50.exact.simulations0			=	c()
	#
	for (k in 1:p) { # k = 1;	# For testing
		theta.with.weights.RejectionABC	=	cbind(theta.table.RejectionABC[,k], rep(1/length(theta.table.RejectionABC[,k]),length(theta.table.RejectionABC[,k])))
		# above using rep(1/length(theta.table.RejectionABC[,k]),length(theta.table.RejectionABC[,k])) in place of theta.DRFweights.table[,k]
		sorted.theta.with.weights.RejectionABC	=	theta.with.weights.RejectionABC[order(theta.with.weights.RejectionABC[,1]), ]
		# Univariate local-polynomial (degree =2) likelihood kernel density estimation (from MCMC samples)
		# using plug-in bandwidth of Sheather and Jones (1991) :
		fit 		= 	kde1d(Samples.theta[,k])# Samples.theta are samples from exact posterior distribution.
		F 			= 	pkde1d(sorted.theta.with.weights.RejectionABC[,1], fit)#cdf estimates of true marginal posterior distribution (from MCMC samples).
		Fk			=	cumsum(sorted.theta.with.weights.RejectionABC[,2])
		Fk.minus.1	=	c(0, Fk[1 : (length(Fk) - 1) ] )
		Dn 			= 	max(apply(cbind(abs(Fk - F), abs(F - Fk.minus.1)), 1, max))# Kolmogorov-Smirnov statistic
		S.2			=	sum(sorted.theta.with.weights.RejectionABC[,2] ^ 2)
		KSstat		=	(1 / S.2) * Dn
		# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)
		KS.Dn.RejectionABC.simulations0	=	c(KS.Dn.RejectionABC.simulations0, 	Dn		)
		KS.stat.RejectionABC.simulations0	=	c(KS.stat.RejectionABC.simulations0, 	KSstat	)
		#
		# Check for 95%(50%) coverage of RejectionABC samples:
		Q.025		=	quantile(theta.table.RejectionABC[,k], .025)
		Q.975		=	quantile(theta.table.RejectionABC[,k], .975)
		Q.25		=	quantile(theta.table.RejectionABC[,k], .25)
		Q.75		=	quantile(theta.table.RejectionABC[,k], .75)
		cover95.RejectionABC.simulations0 	= 	c(cover95.RejectionABC.simulations0, as.numeric((Q.025 <= truth[k]) & (Q.975 >= truth[k])) )
		cover50.RejectionABC.simulations0 	= 	c(cover50.RejectionABC.simulations0, as.numeric((Q.25 <= truth[k])  & (Q.75 >= truth[k])) )
		# Check for 95%(50%) coverage from sample from the exact posterior distribution:
		cover95.exact.simulations0 			= 	c(cover95.exact.simulations0, as.numeric((qkde1d(.025, fit) <= truth[k]) & (qkde1d(.975, fit) >= truth[k])))
		cover50.exact.simulations0 			= 	c(cover50.exact.simulations0, as.numeric((qkde1d(.25, fit) <= truth[k]) & (qkde1d(.75, fit) >= truth[k])))
	}
	cover95.RejectionABC.simulations		=	rbind(cover95.RejectionABC.simulations, 	cover95.RejectionABC.simulations0)
	cover50.RejectionABC.simulations		=	rbind(cover50.RejectionABC.simulations, 	cover50.RejectionABC.simulations0)
	cover95.exact.simulations				=	rbind(cover95.exact.simulations, 	cover95.exact.simulations0)
	cover50.exact.simulations				=	rbind(cover50.exact.simulations, 	cover50.exact.simulations0)
	KS.Dn.RejectionABC.simulations			=	rbind(KS.Dn.RejectionABC.simulations, 		KS.Dn.RejectionABC.simulations0)
	KS.stat.RejectionABC.simulations		=	rbind(KS.stat.RejectionABC.simulations, 	KS.stat.RejectionABC.simulations0)
	L1error.postEtheta.RejectionABC.simulations		=	rbind(L1error.postEtheta.RejectionABC.simulations, 	abs(postEtheta.RejectionABC		- truth) 		)
	L2error.postEtheta.RejectionABC.simulations		=	rbind(L2error.postEtheta.RejectionABC.simulations,    	(postEtheta.RejectionABC 	- truth)^2 	)
	L1error.postMEDtheta.RejectionABC.simulations	=	rbind(L1error.postMEDtheta.RejectionABC.simulations, 	abs(postMEDtheta.RejectionABC	- truth) 		)
	L2error.postMEDtheta.RejectionABC.simulations	=	rbind(L2error.postMEDtheta.RejectionABC.simulations, 		(postMEDtheta.RejectionABC 	- truth)^2 	)
	L1error.postMODEtheta.RejectionABC.simulations	=	rbind(L1error.postMODEtheta.RejectionABC.simulations, 	abs(postMODEtheta.RejectionABC	- truth) 		)
	L2error.postMODEtheta.RejectionABC.simulations	=	rbind(L2error.postMODEtheta.RejectionABC.simulations, 		(postMODEtheta.RejectionABC 	- truth)^2 	)
	L1error.MLEtheta.RejectionABC.simulations			=	rbind(L1error.MLEtheta.RejectionABC.simulations, 		abs(MLEtheta.RejectionABC	- truth) 		)
	L2error.MLEtheta.RejectionABC.simulations			=	rbind(L2error.MLEtheta.RejectionABC.simulations, 			(MLEtheta.RejectionABC 	- truth)^2 	)
	if (fileNum == length(RDataFiles))	{
		postEtheta.RejectionABC.simulations.mean		=	apply(postEtheta.RejectionABC.simulations, 2, mean)
		postEtheta.RejectionABC.simulations.sd		=	apply(postEtheta.RejectionABC.simulations, 2, sd)
		postMEDtheta.RejectionABC.simulations.mean	=	apply(postMEDtheta.RejectionABC.simulations, 2, mean)
		postMEDtheta.RejectionABC.simulations.sd		=	apply(postMEDtheta.RejectionABC.simulations, 2, sd)
		postMODEtheta.RejectionABC.simulations.mean	=	apply(postMODEtheta.RejectionABC.simulations, 2, mean)
		postMODEtheta.RejectionABC.simulations.sd		=	apply(postMODEtheta.RejectionABC.simulations, 2, sd)
		MLEtheta.RejectionABC.simulations.mean		=	apply(MLEtheta.RejectionABC.simulations, 2, mean)
		MLEtheta.RejectionABC.simulations.sd			=	apply(MLEtheta.RejectionABC.simulations, 2, sd)
		postSDtheta.RejectionABC.simulations.mean		=	apply(postSDtheta.RejectionABC.simulations, 2, mean)
		postSDtheta.RejectionABC.simulations.sd		=	apply(postSDtheta.RejectionABC.simulations, 2, sd)
		cover95.exact.simulations.mean			=	apply(cover95.exact.simulations, 2, mean)
		cover50.exact.simulations.mean			=	apply(cover50.exact.simulations, 2, mean)
		cover95.RejectionABC.simulations.mean			=	apply(cover95.RejectionABC.simulations, 2, mean)
		cover50.RejectionABC.simulations.mean			=	apply(cover50.RejectionABC.simulations, 2, mean)
		KS.Dn.RejectionABC.simulations.mean			=	apply(KS.Dn.RejectionABC.simulations, 2, mean)
		KS.Dn.RejectionABC.simulations.sd				=	apply(KS.Dn.RejectionABC.simulations, 2, sd)
		KS.stat.RejectionABC.simulations.mean			=	apply(KS.stat.RejectionABC.simulations, 2, mean)
		KS.stat.RejectionABC.simulations.sd			=	apply(KS.stat.RejectionABC.simulations, 2, sd)
		MAE.postEtheta.RejectionABC.simulations		=	apply(L1error.postEtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
		MSE.postEtheta.RejectionABC.simulations		=	apply(L2error.postEtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
		MAE.postMEDtheta.RejectionABC.simulations		=	apply(L1error.postMEDtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
		MSE.postMEDtheta.RejectionABC.simulations		=	apply(L2error.postMEDtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
		MAE.postMODEtheta.RejectionABC.simulations	=	apply(L1error.postMODEtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
		MSE.postMODEtheta.RejectionABC.simulations	=	apply(L2error.postMODEtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
		MAE.MLEtheta.RejectionABC.simulations			=	apply(L1error.MLEtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
		MSE.MLEtheta.RejectionABC.simulations			=	apply(L2error.MLEtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
	}
}
postEtheta.RejectionABC.simulations.mean	=	round(matrix(postEtheta.RejectionABC.simulations.mean, nrow = 1), 2)
postEtheta.RejectionABC.simulations.sd	=	round(matrix(postEtheta.RejectionABC.simulations.sd, nrow = 1), 2)
postMEDtheta.RejectionABC.simulations.mean=	round(matrix(postMEDtheta.RejectionABC.simulations.mean, nrow = 1), 2)
postMEDtheta.RejectionABC.simulations.sd	=	round(matrix(postMEDtheta.RejectionABC.simulations.sd, nrow = 1), 2)
postMODEtheta.RejectionABC.simulations.mean=	round(matrix(postMODEtheta.RejectionABC.simulations.mean, nrow = 1), 2)
postMODEtheta.RejectionABC.simulations.sd	=	round(matrix(postMODEtheta.RejectionABC.simulations.sd, nrow = 1), 2)
MLEtheta.RejectionABC.simulations.mean	=	round(matrix(MLEtheta.RejectionABC.simulations.mean, nrow = 1), 2)
MLEtheta.RejectionABC.simulations.sd		=	round(matrix(MLEtheta.RejectionABC.simulations.sd, nrow = 1), 2)
postSDtheta.RejectionABC.simulations.mean	=	round(matrix(postSDtheta.RejectionABC.simulations.mean, nrow = 1), 2)
postSDtheta.RejectionABC.simulations.sd	=	round(matrix(postSDtheta.RejectionABC.simulations.sd, nrow = 1), 2)
cover95.RejectionABC.simulations.mean		=	round(matrix(cover95.RejectionABC.simulations.mean, nrow = 1), 2)
cover50.RejectionABC.simulations.mean		=	round(matrix(cover50.RejectionABC.simulations.mean, nrow = 1), 2)
cover95.exact.simulations.mean				=	round(matrix(cover95.exact.simulations.mean, nrow = 1), 2)
cover50.exact.simulations.mean				=	round(matrix(cover50.exact.simulations.mean, nrow = 1), 2)
KS.Dn.RejectionABC.simulations.mean			=	round(matrix(KS.Dn.RejectionABC.simulations.mean, nrow = 1), 2)
KS.Dn.RejectionABC.simulations.sd			=	round(matrix(KS.Dn.RejectionABC.simulations.sd, nrow = 1), 2)
KS.stat.RejectionABC.simulations.mean		=	round(matrix(KS.stat.RejectionABC.simulations.mean, nrow = 1), 2)
KS.stat.RejectionABC.simulations.sd			=	round(matrix(KS.stat.RejectionABC.simulations.sd, nrow = 1), 2)
MAE.postEtheta.RejectionABC.simulations		= 	round(matrix(MAE.postEtheta.RejectionABC.simulations, nrow = 1), 2)
MSE.postEtheta.RejectionABC.simulations		= 	round(matrix(MSE.postEtheta.RejectionABC.simulations, nrow = 1), 2)
MAE.postMEDtheta.RejectionABC.simulations	= 	round(matrix(MAE.postMEDtheta.RejectionABC.simulations, nrow = 1), 2)
MSE.postMEDtheta.RejectionABC.simulations	= 	round(matrix(MSE.postMEDtheta.RejectionABC.simulations, nrow = 1), 2)
MAE.postMODEtheta.RejectionABC.simulations	= 	round(matrix(MAE.postMODEtheta.RejectionABC.simulations, nrow = 1), 2)
MSE.postMODEtheta.RejectionABC.simulations	= 	round(matrix(MSE.postMODEtheta.RejectionABC.simulations, nrow = 1), 2)
MAE.MLEtheta.RejectionABC.simulations		= 	round(matrix(MAE.MLEtheta.RejectionABC.simulations, nrow = 1), 2)
MSE.MLEtheta.RejectionABC.simulations		= 	round(matrix(MSE.MLEtheta.RejectionABC.simulations, nrow = 1), 2)


postEtheta.RejectionABC.simulations.mean.sd 	= 	matrix(paste(postEtheta.RejectionABC.simulations.mean, ' (', postEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
postMEDtheta.RejectionABC.simulations.mean.sd 	= 	matrix(paste(postMEDtheta.RejectionABC.simulations.mean, ' (', postMEDtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
postMODEtheta.RejectionABC.simulations.mean.sd 	= 	matrix(paste(postMODEtheta.RejectionABC.simulations.mean, ' (', postMODEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
MLEtheta.RejectionABC.simulations.mean.sd 		= 	matrix(paste(MLEtheta.RejectionABC.simulations.mean, ' (', MLEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
postSDtheta.RejectionABC.simulations.mean.sd 	=	matrix(paste(postSDtheta.RejectionABC.simulations.mean, '(', postSDtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
cover95.50.RejectionABC.simulations.mean 		= 	matrix(paste(cover95.RejectionABC.simulations.mean, ' (', cover50.RejectionABC.simulations.mean, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
cover95.50.exact.simulations.mean 				= 	matrix(paste(cover95.exact.simulations.mean, 		' (', cover50.exact.simulations.mean, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
KS.Dn.RejectionABC.simulations.mean.sd 			= 	matrix(paste(KS.Dn.RejectionABC.simulations.mean, ' (', KS.Dn.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
KS.stat.RejectionABC.simulations.mean.sd 		= 	matrix(paste(KS.stat.RejectionABC.simulations.mean, ' (', KS.stat.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
MAE.MSE.postEtheta.RejectionABC.simulations 	= 	matrix(paste(MAE.postEtheta.RejectionABC.simulations, ' (', MSE.postEtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
MAE.MSE.postMEDtheta.RejectionABC.simulations 	= 	matrix(paste(MAE.postMEDtheta.RejectionABC.simulations, ' (', MSE.postMEDtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
MAE.MSE.postMODEtheta.RejectionABC.simulations 	= 	matrix(paste(MAE.postMODEtheta.RejectionABC.simulations, ' (', MSE.postMODEtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
MAE.MSE.MLEtheta.RejectionABC.simulations 		= 	matrix(paste(MAE.MLEtheta.RejectionABC.simulations, ' (', MSE.MLEtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))

rownames(postEtheta.RejectionABC.simulations.mean)		=	simConditionNames
rownames(postEtheta.RejectionABC.simulations.sd)		=	simConditionNames
rownames(postEtheta.RejectionABC.simulations.mean.sd)	=	simConditionNames
rownames(postMEDtheta.RejectionABC.simulations.mean)	=	simConditionNames
rownames(postMEDtheta.RejectionABC.simulations.sd)		=	simConditionNames
rownames(postMEDtheta.RejectionABC.simulations.mean.sd)=	simConditionNames
rownames(postMODEtheta.RejectionABC.simulations.mean)	=	simConditionNames
rownames(postMODEtheta.RejectionABC.simulations.sd)		=	simConditionNames
rownames(postMODEtheta.RejectionABC.simulations.mean.sd)=	simConditionNames
rownames(MLEtheta.RejectionABC.simulations.mean)		=	simConditionNames
rownames(MLEtheta.RejectionABC.simulations.sd)			=	simConditionNames
rownames(MLEtheta.RejectionABC.simulations.mean.sd)		=	simConditionNames
rownames(postSDtheta.RejectionABC.simulations.mean)		=	simConditionNames
rownames(postSDtheta.RejectionABC.simulations.sd)		=	simConditionNames
rownames(postSDtheta.RejectionABC.simulations.mean.sd) =	simConditionNames
rownames(cover95.RejectionABC.simulations.mean)			=	simConditionNames
rownames(cover50.RejectionABC.simulations.mean)			=	simConditionNames
rownames(cover95.50.RejectionABC.simulations.mean)		=	simConditionNames
rownames(cover95.exact.simulations.mean)				=	simConditionNames
rownames(cover50.exact.simulations.mean)				=	simConditionNames
rownames(cover95.50.exact.simulations.mean)				=	simConditionNames
rownames(KS.Dn.RejectionABC.simulations.mean)			=	simConditionNames
rownames(KS.Dn.RejectionABC.simulations.sd)				=	simConditionNames
rownames(KS.Dn.RejectionABC.simulations.mean.sd)		=	simConditionNames
rownames(KS.stat.RejectionABC.simulations.mean)			=	simConditionNames
rownames(KS.stat.RejectionABC.simulations.sd)			=	simConditionNames
rownames(KS.stat.RejectionABC.simulations.mean.sd) 		=	simConditionNames
rownames(MAE.postEtheta.RejectionABC.simulations)		=	simConditionNames
rownames(MSE.postEtheta.RejectionABC.simulations)		=	simConditionNames
rownames(MAE.MSE.postEtheta.RejectionABC.simulations) 	=	simConditionNames
rownames(MAE.postMEDtheta.RejectionABC.simulations)		=	simConditionNames
rownames(MSE.postMEDtheta.RejectionABC.simulations)		=	simConditionNames
rownames(MAE.MSE.postMEDtheta.RejectionABC.simulations)	=	simConditionNames
rownames(MAE.postMODEtheta.RejectionABC.simulations)	=	simConditionNames
rownames(MSE.postMODEtheta.RejectionABC.simulations)	=	simConditionNames
rownames(MAE.MSE.postMODEtheta.RejectionABC.simulations)=	simConditionNames
rownames(MAE.MLEtheta.RejectionABC.simulations)			=	simConditionNames
rownames(MSE.MLEtheta.RejectionABC.simulations)			=	simConditionNames
rownames(MAE.MSE.MLEtheta.RejectionABC.simulations)		=	simConditionNames

colnames(postEtheta.RejectionABC.simulations.mean)		=	parameterNames
colnames(postEtheta.RejectionABC.simulations.sd)		=	parameterNames
colnames(postEtheta.RejectionABC.simulations.mean.sd)	=	parameterNames
colnames(postMEDtheta.RejectionABC.simulations.mean)	=	parameterNames
colnames(postMEDtheta.RejectionABC.simulations.sd)		=	parameterNames
colnames(postMEDtheta.RejectionABC.simulations.mean.sd)=	parameterNames
colnames(postMODEtheta.RejectionABC.simulations.mean)	=	parameterNames
colnames(postMODEtheta.RejectionABC.simulations.sd)		=	parameterNames
colnames(postMODEtheta.RejectionABC.simulations.mean.sd)=	parameterNames
colnames(MLEtheta.RejectionABC.simulations.mean)		=	parameterNames
colnames(MLEtheta.RejectionABC.simulations.sd)			=	parameterNames
colnames(MLEtheta.RejectionABC.simulations.mean.sd)		=	parameterNames
colnames(postSDtheta.RejectionABC.simulations.mean)		=	parameterNames
colnames(postSDtheta.RejectionABC.simulations.sd)		=	parameterNames
colnames(postSDtheta.RejectionABC.simulations.mean.sd)	=	parameterNames
colnames(cover95.RejectionABC.simulations.mean)			=	parameterNames
colnames(cover50.RejectionABC.simulations.mean)			=	parameterNames
colnames(cover95.50.RejectionABC.simulations.mean)		=	parameterNames
colnames(cover95.exact.simulations.mean)				=	parameterNames
colnames(cover50.exact.simulations.mean)				=	parameterNames
colnames(cover95.50.exact.simulations.mean)				=	parameterNames
colnames(KS.Dn.RejectionABC.simulations.mean)			=	parameterNames
colnames(KS.Dn.RejectionABC.simulations.sd)				=	parameterNames
colnames(KS.Dn.RejectionABC.simulations.mean.sd)		=	parameterNames
colnames(KS.stat.RejectionABC.simulations.mean)			=	parameterNames
colnames(KS.stat.RejectionABC.simulations.sd)			=	parameterNames
colnames(KS.stat.RejectionABC.simulations.mean.sd) 		=	parameterNames
colnames(MAE.postEtheta.RejectionABC.simulations)		=	parameterNames
colnames(MSE.postEtheta.RejectionABC.simulations)		=	parameterNames
colnames(MAE.MSE.postEtheta.RejectionABC.simulations)	=	parameterNames
colnames(MAE.postMEDtheta.RejectionABC.simulations)		=	parameterNames
colnames(MSE.postMEDtheta.RejectionABC.simulations)		=	parameterNames
colnames(MAE.MSE.postMEDtheta.RejectionABC.simulations)=	parameterNames
colnames(MAE.postMODEtheta.RejectionABC.simulations)	=	parameterNames
colnames(MSE.postMODEtheta.RejectionABC.simulations)	=	parameterNames
colnames(MAE.MSE.postMODEtheta.RejectionABC.simulations)=	parameterNames
colnames(MAE.MLEtheta.RejectionABC.simulations)			=	parameterNames
colnames(MSE.MLEtheta.RejectionABC.simulations)			=	parameterNames
colnames(MAE.MSE.MLEtheta.RejectionABC.simulations)		=	parameterNames

#postEtheta.RejectionABC.simulations.mean
#postEtheta.RejectionABC.simulations.sd
noquote(postEtheta.RejectionABC.simulations.mean.sd)
#postMEDtheta.RejectionABC.simulations.mean
#postMEDtheta.RejectionABC.simulations.sd
noquote(postMEDtheta.RejectionABC.simulations.mean.sd)
#postMODEtheta.RejectionABC.simulations.mean
#postMODEtheta.RejectionABC.simulations.sd
noquote(postMODEtheta.RejectionABC.simulations.mean.sd)
#MLEtheta.RejectionABC.simulations.mean
#MLEtheta.RejectionABC.simulations.sd
noquote(MLEtheta.RejectionABC.simulations.mean.sd)
#postSDtheta.RejectionABC.simulations.mean
#postSDtheta.RejectionABC.simulations.sd
noquote(postSDtheta.RejectionABC.simulations.mean.sd)
#cover95.RejectionABC.simulations.mean
#cover50.RejectionABC.simulations.mean
noquote(cover95.50.RejectionABC.simulations.mean)
#cover95.exact.simulations.mean
#cover50.exact.simulations.mean
noquote(cover95.50.exact.simulations.mean)
#KS.Dn.RejectionABC.simulations.mean
#KS.Dn.RejectionABC.simulations.sd
noquote(KS.Dn.RejectionABC.simulations.mean.sd)
#KS.stat.RejectionABC.simulations.mean
#KS.stat.RejectionABC.simulations.sd
noquote(KS.stat.RejectionABC.simulations.mean.sd)
#MAE.postEtheta.RejectionABC.simulations
#MSE.postEtheta.RejectionABC.simulations
noquote(MAE.MSE.postEtheta.RejectionABC.simulations)
#MAE.postMEDtheta.RejectionABC.simulations
#MSE.postMEDtheta.RejectionABC.simulations
noquote(MAE.MSE.postMEDtheta.RejectionABC.simulations)
#MAE.postMODEtheta.RejectionABC.simulations
#MSE.postMODEtheta.RejectionABC.simulations
noquote(MAE.MSE.postMODEtheta.RejectionABC.simulations)
#MAE.MLEtheta.RejectionABC.simulations
#MSE.MLEtheta.RejectionABC.simulations
noquote(MAE.MSE.MLEtheta.RejectionABC.simulations)