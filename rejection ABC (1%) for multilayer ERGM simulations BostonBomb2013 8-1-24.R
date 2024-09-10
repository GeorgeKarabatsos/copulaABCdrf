rm(list = ls())
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
library('ks')	# for kde()
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('kde1d');			library('kde1d')

for (simStudy in 1 : 2)	{ # simStudy = 1 # For testing.
	if (simStudy == 1)	{
		postEtheta.RejectionABC.simulations.mean	= 	matrix(NA, nrow = 2, ncol = 18)
		postEtheta.RejectionABC.simulations.sd		= 	matrix(NA, nrow = 2, ncol = 18)
		postMEDtheta.RejectionABC.simulations.mean	= 	matrix(NA, nrow = 2, ncol = 18)
		postMEDtheta.RejectionABC.simulations.sd	= 	matrix(NA, nrow = 2, ncol = 18)
		postMODEtheta.RejectionABC.simulations.mean	= 	matrix(NA, nrow = 2, ncol = 18)
		postMODEtheta.RejectionABC.simulations.sd	= 	matrix(NA, nrow = 2, ncol = 18)
		MLEtheta.RejectionABC.simulations.mean		= 	matrix(NA, nrow = 2, ncol = 18)
		MLEtheta.RejectionABC.simulations.sd		= 	matrix(NA, nrow = 2, ncol = 18)
		postSDtheta.RejectionABC.simulations.mean	= 	matrix(NA, nrow = 2, ncol = 18)
		postSDtheta.RejectionABC.simulations.sd		= 	matrix(NA, nrow = 2, ncol = 18)
		cover95.RejectionABC.simulations.mean		= 	matrix(NA, nrow = 2, ncol = 18)
		cover50.RejectionABC.simulations.mean		= 	matrix(NA, nrow = 2, ncol = 18)
		MAE.postEtheta.RejectionABC.simulations 	= 	matrix(NA, nrow = 2, ncol = 18)
		MSE.postEtheta.RejectionABC.simulations 	= 	matrix(NA, nrow = 2, ncol = 18)
		MAE.postMEDtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 2, ncol = 18)
		MSE.postMEDtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 2, ncol = 18)
		MAE.postMODEtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 2, ncol = 18)
		MSE.postMODEtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 2, ncol = 18)
		MAE.MLEtheta.RejectionABC.simulations		= 	matrix(NA, nrow = 2, ncol = 18)
		MSE.MLEtheta.RejectionABC.simulations		= 	matrix(NA, nrow = 2, ncol = 18)
		simConditionNames	=	c(	'ERGMml,  n = 100 nodes, n_sim = 10',
									'ERGMml,  n = 100 nodes, n_sim = 20')
		rownames(postEtheta.RejectionABC.simulations.mean)	=	simConditionNames
		rownames(postEtheta.RejectionABC.simulations.sd)	=	simConditionNames
		rownames(postMEDtheta.RejectionABC.simulations.mean)	=	simConditionNames
		rownames(postMEDtheta.RejectionABC.simulations.sd)	=	simConditionNames
		rownames(postMODEtheta.RejectionABC.simulations.mean)=	simConditionNames
		rownames(postMODEtheta.RejectionABC.simulations.sd)	=	simConditionNames
		rownames(MLEtheta.RejectionABC.simulations.mean)	=	simConditionNames
		rownames(MLEtheta.RejectionABC.simulations.sd)		=	simConditionNames
		rownames(postSDtheta.RejectionABC.simulations.mean)	=	simConditionNames
		rownames(postSDtheta.RejectionABC.simulations.sd)	=	simConditionNames
		rownames(cover95.RejectionABC.simulations.mean)		=	simConditionNames
		rownames(cover50.RejectionABC.simulations.mean)		=	simConditionNames
		rownames(MAE.postEtheta.RejectionABC.simulations)	=	simConditionNames
		rownames(MSE.postEtheta.RejectionABC.simulations)	=	simConditionNames
		rownames(MAE.postMEDtheta.RejectionABC.simulations)	=	simConditionNames
		rownames(MSE.postMEDtheta.RejectionABC.simulations)	=	simConditionNames
		rownames(MAE.postMODEtheta.RejectionABC.simulations)	=	simConditionNames
		rownames(MSE.postMODEtheta.RejectionABC.simulations)	=	simConditionNames
		rownames(MAE.MLEtheta.RejectionABC.simulations)		=	simConditionNames
		rownames(MSE.MLEtheta.RejectionABC.simulations)		=	simConditionNames
		colnames(postEtheta.RejectionABC.simulations.mean)	=	paste('theta', 1 : 18, sep = '')
		colnames(postEtheta.RejectionABC.simulations.sd)	=	paste('theta', 1 : 18, sep = '')
		colnames(postMEDtheta.RejectionABC.simulations.mean)	=	paste('theta', 1 : 18, sep = '')
		colnames(postMEDtheta.RejectionABC.simulations.sd)	=	paste('theta', 1 : 18, sep = '')
		colnames(postMODEtheta.RejectionABC.simulations.mean)=	paste('theta', 1 : 18, sep = '')
		colnames(postMODEtheta.RejectionABC.simulations.sd)	=	paste('theta', 1 : 18, sep = '')
		colnames(MLEtheta.RejectionABC.simulations.mean)	=	paste('theta', 1 : 18, sep = '')
		colnames(MLEtheta.RejectionABC.simulations.sd)		=	paste('theta', 1 : 18, sep = '')
		colnames(postSDtheta.RejectionABC.simulations.mean)	=	paste('theta', 1 : 18, sep = '')
		colnames(postSDtheta.RejectionABC.simulations.sd)	=	paste('theta', 1 : 18, sep = '')
		colnames(cover95.RejectionABC.simulations.mean)		=	paste('theta', 1 : 18, sep = '')
		colnames(cover50.RejectionABC.simulations.mean)		=	paste('theta', 1 : 18, sep = '')
		colnames(MAE.postEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
		colnames(MSE.postEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
		colnames(MAE.postMEDtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
		colnames(MSE.postMEDtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
		colnames(MAE.postMODEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
		colnames(MSE.postMODEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
		colnames(MAE.MLEtheta.RejectionABC.simulations)		=	paste('theta', 1 : 18, sep = '')
		colnames(MSE.MLEtheta.RejectionABC.simulations)		=	paste('theta', 1 : 18, sep = '')
	}

	if (simStudy == 1) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\Multilayer-valued ERGM\\n = 100 n_sim = 10')}
	if (simStudy == 2) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\Multilayer-valued ERGM\\n = 100 n_sim = 20')}

	RDataFiles		=	list.files(path = ".", pattern = "\\.RData$")

	postEtheta.RejectionABC.simulations 			=	c()
	postMEDtheta.RejectionABC.simulations 			=	c()
	postMODEtheta.RejectionABC.simulations 			=	c()
	MLEtheta.RejectionABC.simulations 				=	c()
	postSDtheta.RejectionABC.simulations 			=	c()
	cover95.RejectionABC.simulations				=	c()
	cover50.RejectionABC.simulations				=	c()
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

	# Loop through the 10 replicas for the given network model, n nodes, and n_sim.
	for (fileNum in 1 : length(RDataFiles)){# fileNum = 1 # For testing.
		load(RDataFiles[fileNum])
		#
		# Do the rejection ABC analysis:
		euclideanDistances		=	 as.matrix(pdist(sy.table, Y = matrix(X.test, nrow = 1) ))
		theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .01),]
		# Above, retained 1% of thetas with summary statistics sy nearest to observed dataset summaries sx.		
		#
		postEtheta.RejectionABC	=	colMeans(theta.table.RejectionABC)
		postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
		#
		# kde() gives the Kernel density estimate for 1- to 6-dimensional data:
		# See:   https://www.rdocumentation.org/packages/ks/versions/1.10.7/topics/kde
		# Therefore, given that the model parameters is 18 dimensional, will not use kde().
		# post.kde.RejectionABC	=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate}
		#
		# Instead, will use a product kernel multidimensional density estimator (e.g. Wasserman, 2006, All of Nonparametric Statistics):
		for (k in 1:length(truth)) { # k = 1;	# For testing
			if (k == 1) {post.kde.RejectionABC = matrix(NA, nrow = nrow(theta.table.RejectionABC), ncol =  length(truth))}
			fit.RejectionABC			= 	kde1d(theta.table.RejectionABC[,k])
			post.kde.RejectionABC[,k]	=	dkde1d(theta.table.RejectionABC[,k], fit.RejectionABC)	}
		post.kde.RejectionABC			=	apply(post.kde.RejectionABC, 1, prod)
		#
		postMODEtheta.RejectionABC		= 	theta.table.RejectionABC[which.max(post.kde.RejectionABC),]		
		priorPDFs.RejectionABC		=	dmvnorm(theta.table.RejectionABC, mean = muPrior, sigma = SigmaPrior)
		likelihoods.RejectionABC		=	post.kde.RejectionABC / priorPDFs.RejectionABC
		MLEtheta.RejectionABC			=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]
		postSDtheta.RejectionABC		=	apply(theta.table.RejectionABC, 2, sd)	
		#
		postEtheta.RejectionABC.simulations	=	rbind(postEtheta.RejectionABC.simulations, postEtheta.RejectionABC)
		postMEDtheta.RejectionABC.simulations	=	rbind(postMEDtheta.RejectionABC.simulations, postMEDtheta.RejectionABC)
		postMODEtheta.RejectionABC.simulations	=	rbind(postMODEtheta.RejectionABC.simulations, postMODEtheta.RejectionABC)
		MLEtheta.RejectionABC.simulations		=	rbind(MLEtheta.RejectionABC.simulations, MLEtheta.RejectionABC)
		postSDtheta.RejectionABC.simulations	=	rbind(postSDtheta.RejectionABC.simulations, postSDtheta.RejectionABC)
		
		# Compute 95%(50% coverage for each individual parameter.
		postQtheta.RejectionABC				= 	apply(theta.table.RejectionABC, 2, quantile, probs = c(.025, .25, .75, .975))
		cover95.RejectionABC.simulations0		=	c()
		cover50.RejectionABC.simulations0		=	c()
		#
		for (k in 1:length(truth)) { # k = 1;	# For testing
			if (k == 1){postQtheta.RejectionABC= 	apply(theta.table.RejectionABC, 2, quantile, probs = c(.025, .25, .75, .975))}
			# Check for 95%(50%) coverage:
			q.lo	=	postQtheta.RejectionABC[1,k]
			q.hi	=	postQtheta.RejectionABC[4,k]
			cover95.RejectionABC.simulations0 	= 	c(cover95.RejectionABC.simulations0, as.numeric(q.lo <= truth[k] & q.hi >= truth[k]))
			q.lo	=	postQtheta.RejectionABC[2,k];
			q.hi	=	postQtheta.RejectionABC[3,k]
			cover50.RejectionABC.simulations0 	= 	c(cover50.RejectionABC.simulations0, as.numeric(q.lo <= truth[k] & q.hi >= truth[k]))
		}
		cover95.RejectionABC.simulations				=	rbind(cover95.RejectionABC.simulations, 	cover95.RejectionABC.simulations0)
		cover50.RejectionABC.simulations				=	rbind(cover50.RejectionABC.simulations, 	cover50.RejectionABC.simulations0)
		L1error.postEtheta.RejectionABC.simulations		=	rbind(L1error.postEtheta.RejectionABC.simulations, 	abs(postEtheta.RejectionABC		- truth) 		)
		L2error.postEtheta.RejectionABC.simulations		=	rbind(L2error.postEtheta.RejectionABC.simulations,    	(postEtheta.RejectionABC 	- truth)^2 	)
		L1error.postMEDtheta.RejectionABC.simulations	=	rbind(L1error.postMEDtheta.RejectionABC.simulations, 	abs(postMEDtheta.RejectionABC	- truth) 		)
		L2error.postMEDtheta.RejectionABC.simulations	=	rbind(L2error.postMEDtheta.RejectionABC.simulations, 		(postMEDtheta.RejectionABC 	- truth)^2 	)
		L1error.postMODEtheta.RejectionABC.simulations	=	rbind(L1error.postMODEtheta.RejectionABC.simulations, 	abs(postMODEtheta.RejectionABC	- truth) 		)
		L2error.postMODEtheta.RejectionABC.simulations	=	rbind(L2error.postMODEtheta.RejectionABC.simulations, 		(postMODEtheta.RejectionABC 	- truth)^2 	)
		L1error.MLEtheta.RejectionABC.simulations		=	rbind(L1error.MLEtheta.RejectionABC.simulations, 		abs(MLEtheta.RejectionABC	- truth) 		)
		L2error.MLEtheta.RejectionABC.simulations		=	rbind(L2error.MLEtheta.RejectionABC.simulations, 			(MLEtheta.RejectionABC 	- truth)^2 	)
		#
		if (fileNum == length(RDataFiles))	{
			postEtheta.RejectionABC.simulations.mean0	=	matrix(apply(postEtheta.RejectionABC.simulations, 2, mean), nrow = 1)
			postEtheta.RejectionABC.simulations.sd0		=	matrix(apply(postEtheta.RejectionABC.simulations, 2, sd), nrow = 1)
			postMEDtheta.RejectionABC.simulations.mean0	=	matrix(apply(postMEDtheta.RejectionABC.simulations, 2, mean), nrow = 1)
			postMEDtheta.RejectionABC.simulations.sd0	=	matrix(apply(postMEDtheta.RejectionABC.simulations, 2, sd), nrow = 1)
			postMODEtheta.RejectionABC.simulations.mean0=	matrix(apply(postMODEtheta.RejectionABC.simulations, 2, mean), nrow = 1)
			postMODEtheta.RejectionABC.simulations.sd0	=	matrix(apply(postMODEtheta.RejectionABC.simulations, 2, sd), nrow = 1)
			MLEtheta.RejectionABC.simulations.mean0		=	matrix(apply(MLEtheta.RejectionABC.simulations, 2, mean), nrow = 1)
			MLEtheta.RejectionABC.simulations.sd0		=	matrix(apply(MLEtheta.RejectionABC.simulations, 2, sd), nrow = 1)
			postSDtheta.RejectionABC.simulations.mean0	=	matrix(apply(postSDtheta.RejectionABC.simulations, 2, mean), nrow = 1)
			postSDtheta.RejectionABC.simulations.sd0	=	matrix(apply(postSDtheta.RejectionABC.simulations, 2, sd), nrow = 1)
			cover95.RejectionABC.simulations.mean0		=	matrix(apply(cover95.RejectionABC.simulations, 2, mean), nrow = 1)
			cover50.RejectionABC.simulations.mean0		=	matrix(apply(cover50.RejectionABC.simulations, 2, mean), nrow = 1)
			MAE.postEtheta.RejectionABC.simulations0	=	matrix(apply(L1error.postEtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Absolute Error (MAE)
			MSE.postEtheta.RejectionABC.simulations0	=	matrix(apply(L2error.postEtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Squared Error (MSE)
			MAE.postMEDtheta.RejectionABC.simulations0	=	matrix(apply(L1error.postMEDtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Absolute Error (MAE)
			MSE.postMEDtheta.RejectionABC.simulations0	=	matrix(apply(L2error.postMEDtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Squared Error (MSE)
			MAE.postMODEtheta.RejectionABC.simulations0	=	matrix(apply(L1error.postMODEtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Absolute Error (MAE)
			MSE.postMODEtheta.RejectionABC.simulations0	=	matrix(apply(L2error.postMODEtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Squared Error (MSE)
			MAE.MLEtheta.RejectionABC.simulations0		=	matrix(apply(L1error.MLEtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Absolute Error (MAE)
			MSE.MLEtheta.RejectionABC.simulations0		=	matrix(apply(L2error.MLEtheta.RejectionABC.simulations, 2, mean), nrow = 1) # Mean Squared Error (MSE)
		}
	}
	postEtheta.RejectionABC.simulations.mean[simStudy, 1:length(truth)]	=	round(postEtheta.RejectionABC.simulations.mean0, 2)
	postEtheta.RejectionABC.simulations.sd[simStudy, 1:length(truth)]		=	round(postEtheta.RejectionABC.simulations.sd0, 2)
	postMEDtheta.RejectionABC.simulations.mean[simStudy, 1:length(truth)]	=	round(postMEDtheta.RejectionABC.simulations.mean0, 2)
	postMEDtheta.RejectionABC.simulations.sd[simStudy, 1:length(truth)]	=	round(postMEDtheta.RejectionABC.simulations.sd0, 2)
	postMODEtheta.RejectionABC.simulations.mean[simStudy, 1:length(truth)]=	round(postMODEtheta.RejectionABC.simulations.mean0, 2)
	postMODEtheta.RejectionABC.simulations.sd[simStudy, 1:length(truth)]	=	round(postMODEtheta.RejectionABC.simulations.sd0, 2)
	MLEtheta.RejectionABC.simulations.mean[simStudy, 1:length(truth)]		=	round(MLEtheta.RejectionABC.simulations.mean0, 2)
	MLEtheta.RejectionABC.simulations.sd[simStudy, 1:length(truth)]		=	round(MLEtheta.RejectionABC.simulations.sd0, 2)
	postSDtheta.RejectionABC.simulations.mean[simStudy, 1:length(truth)]	=	round(postSDtheta.RejectionABC.simulations.mean0, 2)
	postSDtheta.RejectionABC.simulations.sd[simStudy, 1:length(truth)]	=	round(postSDtheta.RejectionABC.simulations.sd0, 2)
	cover95.RejectionABC.simulations.mean[simStudy, 1:length(truth)]		=	round(cover95.RejectionABC.simulations.mean0, 2)
	cover50.RejectionABC.simulations.mean[simStudy, 1:length(truth)]		=	round(cover50.RejectionABC.simulations.mean0, 2)
	MAE.postEtheta.RejectionABC.simulations[simStudy, 1:length(truth)]	= 	round(MAE.postEtheta.RejectionABC.simulations0, 2)
	MSE.postEtheta.RejectionABC.simulations[simStudy, 1:length(truth)]	= 	round(MSE.postEtheta.RejectionABC.simulations0, 2)
	MAE.postMEDtheta.RejectionABC.simulations[simStudy, 1:length(truth)]	= 	round(MAE.postMEDtheta.RejectionABC.simulations0, 2)
	MSE.postMEDtheta.RejectionABC.simulations[simStudy, 1:length(truth)]	= 	round(MSE.postMEDtheta.RejectionABC.simulations0, 2)
	MAE.postMODEtheta.RejectionABC.simulations[simStudy, 1:length(truth)]	= 	round(MAE.postMODEtheta.RejectionABC.simulations0, 2)
	MSE.postMODEtheta.RejectionABC.simulations[simStudy, 1:length(truth)]	= 	round(MSE.postMODEtheta.RejectionABC.simulations0, 2)
	MAE.MLEtheta.RejectionABC.simulations[simStudy, 1:length(truth)]		= 	round(MAE.MLEtheta.RejectionABC.simulations0, 2)
	MSE.MLEtheta.RejectionABC.simulations[simStudy, 1:length(truth)]		= 	round(MSE.MLEtheta.RejectionABC.simulations0, 2)	
}

#postEtheta.RejectionABC.simulations.mean
#postEtheta.RejectionABC.simulations.sd
postEtheta.RejectionABC.simulations.meansd	=	matrix(paste(postEtheta.RejectionABC.simulations.mean, ' (', postEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(postEtheta.RejectionABC.simulations.meansd)	=	simConditionNames
colnames(postEtheta.RejectionABC.simulations.meansd)	=	paste('theta', 1 : 18, sep = '')
noquote(postEtheta.RejectionABC.simulations.meansd)
#
#postMEDtheta.RejectionABC.simulations.mean
#postMEDtheta.RejectionABC.simulations.sd	
postMEDtheta.RejectionABC.simulations.meansd	=	matrix(paste(postMEDtheta.RejectionABC.simulations.mean, ' (', postMEDtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(postMEDtheta.RejectionABC.simulations.meansd)	=	simConditionNames
colnames(postMEDtheta.RejectionABC.simulations.meansd)	=	paste('theta', 1 : 18, sep = '')
noquote(postMEDtheta.RejectionABC.simulations.meansd)
#
#postMODEtheta.RejectionABC.simulations.mean	
#postMODEtheta.RejectionABC.simulations.sd	
postMODEtheta.RejectionABC.simulations.meansd	=	matrix(paste(postMODEtheta.RejectionABC.simulations.mean, ' (', postMODEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(postMODEtheta.RejectionABC.simulations.meansd)	=	simConditionNames
colnames(postMODEtheta.RejectionABC.simulations.meansd)	=	paste('theta', 1 : 18, sep = '')
noquote(postMODEtheta.RejectionABC.simulations.meansd)
#
#MLEtheta.RejectionABC.simulations.mean	
#MLEtheta.RejectionABC.simulations.sd		
MLEtheta.RejectionABC.simulations.meansd	=	matrix(paste(MLEtheta.RejectionABC.simulations.mean, ' (', MLEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(MLEtheta.RejectionABC.simulations.meansd)	=	simConditionNames
colnames(MLEtheta.RejectionABC.simulations.meansd)	=	paste('theta', 1 : 18, sep = '')
noquote(MLEtheta.RejectionABC.simulations.meansd)
#
#postSDtheta.RejectionABC.simulations.mean	
#postSDtheta.RejectionABC.simulations.sd		
postSDtheta.RejectionABC.simulations.meansd	=	matrix(paste(postSDtheta.RejectionABC.simulations.mean, ' (', postSDtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(postSDtheta.RejectionABC.simulations.meansd)	=	simConditionNames
colnames(postSDtheta.RejectionABC.simulations.meansd)	=	paste('theta', 1 : 18, sep = '')
noquote(postSDtheta.RejectionABC.simulations.meansd)
#
#cover95.RejectionABC.simulations.mean		
#cover50.RejectionABC.simulations.mean		
cover95.50.RejectionABC.simulations.mean	=	matrix(paste(cover95.RejectionABC.simulations.mean	, ' (', cover50.RejectionABC.simulations.mean, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(cover95.50.RejectionABC.simulations.mean)	=	simConditionNames
colnames(cover95.50.RejectionABC.simulations.mean)	=	paste('theta', 1 : 18, sep = '')
noquote(cover95.50.RejectionABC.simulations.mean)
#MAE.postEtheta.RejectionABC.simulations		
#MSE.postEtheta.RejectionABC.simulations	
MAE.MSE.postEtheta.RejectionABC.simulations		=	matrix(paste(MAE.postEtheta.RejectionABC.simulations, ' (', MSE.postEtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(MAE.MSE.postEtheta.RejectionABC.simulations)	=	simConditionNames
colnames(MAE.MSE.postEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
noquote(MAE.MSE.postEtheta.RejectionABC.simulations)
#MAE.postMEDtheta.RejectionABC.simulations	
#MSE.postMEDtheta.RejectionABC.simulations	
MAE.MSE.postMEDtheta.RejectionABC.simulations	=	matrix(paste(MAE.postMEDtheta.RejectionABC.simulations, ' (', MSE.postMEDtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(MAE.MSE.postMEDtheta.RejectionABC.simulations)	=	simConditionNames
colnames(MAE.MSE.postMEDtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
noquote(MAE.MSE.postMEDtheta.RejectionABC.simulations)
#MAE.postMODEtheta.RejectionABC.simulations
#MSE.postMODEtheta.RejectionABC.simulations	
MAE.MSE.postMODEtheta.RejectionABC.simulations	=	matrix(paste(MAE.postMODEtheta.RejectionABC.simulations, ' (', MSE.postMODEtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(MAE.MSE.postMODEtheta.RejectionABC.simulations)	=	simConditionNames
colnames(MAE.MSE.postMODEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
noquote(MAE.MSE.postMODEtheta.RejectionABC.simulations)
#MAE.MLEtheta.RejectionABC.simulations		
#MSE.MLEtheta.RejectionABC.simulations		
MAE.MSE.MLEtheta.RejectionABC.simulations	=	matrix(paste(MAE.MLEtheta.RejectionABC.simulations, ' (', MSE.MLEtheta.RejectionABC.simulations, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean),ncol=ncol(postEtheta.RejectionABC.simulations.mean))
rownames(MAE.MSE.MLEtheta.RejectionABC.simulations)	=	simConditionNames
colnames(MAE.MSE.MLEtheta.RejectionABC.simulations)	=	paste('theta', 1 : 18, sep = '')
noquote(MAE.MSE.MLEtheta.RejectionABC.simulations)

# Present transposed tables of the results: 
noquote(t(postEtheta.RejectionABC.simulations.meansd))
noquote(t(postMEDtheta.RejectionABC.simulations.meansd))
noquote(t(postMODEtheta.RejectionABC.simulations.meansd))
noquote(t(MLEtheta.RejectionABC.simulations.meansd))
noquote(t(postSDtheta.RejectionABC.simulations.meansd))
noquote(t(MAE.MSE.postEtheta.RejectionABC.simulations))
noquote(t(MAE.MSE.postMEDtheta.RejectionABC.simulations))
noquote(t(MAE.MSE.postMODEtheta.RejectionABC.simulations))
noquote(t(MAE.MSE.MLEtheta.RejectionABC.simulations))
