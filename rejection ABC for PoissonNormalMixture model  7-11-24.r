rm(list = ls())
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
library('ks')	# for kde()

for (simStudy in 1 : 8)	{
	#   simStudy = 1 # for testing
	if (simStudy == 1)	{
		parameterNames		=	c('lambda', 'mu')
		simConditionNames	=	c(	'n_sim = 10', 'n_sim = 25', 'n_sim = 33', 'n_sim = 50',
									'n_sim = 66', 'n_sim = 75', 'n_sim = 90', 'n_sim = 100')
		postEtheta.RejectionABC.simulations.mean		=	matrix(NA, nrow = 8, ncol = 2)
		postEtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 8, ncol = 2)
		postMEDtheta.RejectionABC.simulations.mean	=	matrix(NA, nrow = 8, ncol = 2)
		postMEDtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 8, ncol = 2)
		postMODEtheta.RejectionABC.simulations.mean	=	matrix(NA, nrow = 8, ncol = 2)
		postMODEtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 8, ncol = 2)
		MLEtheta.RejectionABC.simulations.mean		=	matrix(NA, nrow = 8, ncol = 2)
		MLEtheta.RejectionABC.simulations.sd			=	matrix(NA, nrow = 8, ncol = 2)
		postSDtheta.RejectionABC.simulations.mean		=	matrix(NA, nrow = 8, ncol = 2)
		postSDtheta.RejectionABC.simulations.sd		=	matrix(NA, nrow = 8, ncol = 2)
		cover95.RejectionABC.simulations.mean			=	matrix(NA, nrow = 8, ncol = 2)
		cover50.RejectionABC.simulations.mean			=	matrix(NA, nrow = 8, ncol = 2)
		cover95.exact.simulations.mean				=	matrix(NA, nrow = 8, ncol = 2)
		cover50.exact.simulations.mean				=	matrix(NA, nrow = 8, ncol = 2)
		KS.Dn.RejectionABC.simulations.mean			=	matrix(NA, nrow = 8, ncol = 2)
		KS.Dn.RejectionABC.simulations.sd				=	matrix(NA, nrow = 8, ncol = 2)
		KS.stat.RejectionABC.simulations.mean			=	matrix(NA, nrow = 8, ncol = 2)
		KS.stat.RejectionABC.simulations.sd			=	matrix(NA, nrow = 8, ncol = 2)
		KS.RejectionABC.simulations						=	matrix(NA, nrow = 8, ncol = 2)
		KS.stat.RejectionABC.simulations				=	matrix(NA, nrow = 8, ncol = 2)
		MAE.postEtheta.RejectionABC.simulations 		= 	matrix(NA, nrow = 8, ncol = 2)
		MSE.postEtheta.RejectionABC.simulations 		= 	matrix(NA, nrow = 8, ncol = 2)
		MAE.postMEDtheta.RejectionABC.simulations		= 	matrix(NA, nrow = 8, ncol = 2)
		MSE.postMEDtheta.RejectionABC.simulations		= 	matrix(NA, nrow = 8, ncol = 2)
		MAE.postMODEtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 8, ncol = 2)
		MSE.postMODEtheta.RejectionABC.simulations	= 	matrix(NA, nrow = 8, ncol = 2)
		MAE.MLEtheta.RejectionABC.simulations			= 	matrix(NA, nrow = 8, ncol = 2)
		MSE.MLEtheta.RejectionABC.simulations			= 	matrix(NA, nrow = 8, ncol = 2)
	}	
	if (simStudy == 1) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 10')}
	if (simStudy == 2) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 25')}
	if (simStudy == 3) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 33')}
	if (simStudy == 4) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 50')}
	if (simStudy == 5) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 66')}
	if (simStudy == 6)	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 75')}
	if (simStudy == 7) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 90')}
	if (simStudy == 8) 	{setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 100')}

	RDataFiles		=	list.files(path = ".", pattern = "\\.RData$")

	postEtheta.RejectionABC.simulations 				=	c()
	postMEDtheta.RejectionABC.simulations 			=	c()
	postMODEtheta.RejectionABC.simulations 			=	c()
	MLEtheta.RejectionABC.simulations 					=	c()
	postSDtheta.RejectionABC.simulations 				=	c()
	cover95.RejectionABC.simulations					=	c()
	cover50.RejectionABC.simulations					=	c()
	cover95.exact.simulations							=	c()
	cover50.exact.simulations							=	c()	
	KS.Dn.RejectionABC.simulations						=	c()
	KS.stat.RejectionABC.simulations					=	c()
	L1error.postEtheta.RejectionABC.simulations		=	c()
	L2error.postEtheta.RejectionABC.simulations		=	c()
	L1error.postMEDtheta.RejectionABC.simulations	=	c()
	L2error.postMEDtheta.RejectionABC.simulations	=	c()
	L1error.postMODEtheta.RejectionABC.simulations	= 	c()
	L2error.postMODEtheta.RejectionABC.simulations	= 	c()
	L1error.MLEtheta.RejectionABC.simulations			= 	c()
	L2error.MLEtheta.RejectionABC.simulations			= 	c()

	# Loop through the 10 replicas for the given n_sim simulation condition 
	for (fileNum in 1 : length(RDataFiles))	{# fileNum = 1 # For testing.
		load(RDataFiles[fileNum])
		#
		# Do the rejection ABC analysis:
		euclideanDistances			=	as.matrix(pdist(sy.table, Y = sx))
		theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .01),]
		# Above, retained 1%,2%, or 3% of thetas with summary statistics sy nearest to observed dataset summaries sx.
		
		postEtheta.RejectionABC		=	colMeans(theta.table.RejectionABC)
		postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
		post.kde2d.RejectionABC		=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate
		postMODEtheta.RejectionABC	= 	theta.table.RejectionABC[which.max(post.kde2d.RejectionABC),]
		priorPDFs.RejectionABC		=	apply(cbind(dgamma(theta.table.RejectionABC[,1], shape = 1/2, rate = 0.1), dunif(theta.table.RejectionABC[,2], min = -10, max = 10)), 1, prod)
		likelihoods.RejectionABC	=	post.kde2d.RejectionABC / priorPDFs.RejectionABC
		MLEtheta.RejectionABC		=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]
		postSDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, sd)	
		#
		postEtheta.RejectionABC.simulations	=	rbind(postEtheta.RejectionABC.simulations, postEtheta.RejectionABC)
		postMEDtheta.RejectionABC.simulations	=	rbind(postMEDtheta.RejectionABC.simulations, postMEDtheta.RejectionABC)
		postMODEtheta.RejectionABC.simulations=	rbind(postMODEtheta.RejectionABC.simulations, postMODEtheta.RejectionABC)
		MLEtheta.RejectionABC.simulations		=	rbind(MLEtheta.RejectionABC.simulations, MLEtheta.RejectionABC)
		postSDtheta.RejectionABC.simulations	=	rbind(postSDtheta.RejectionABC.simulations, postSDtheta.RejectionABC)

		# Do weighted one-sample (two-tailed) Kolmogorov Smirnov tests (Monahan,2011,p.358)
		# and compute 95%(50% coverage for each individual parameter.
		KS.Dn.RejectionABC.simulations0			=	c()
		KS.stat.RejectionABC.simulations0		=	c()
		cover95.RejectionABC.simulations0 		= 	c()
		cover50.RejectionABC.simulations0 		= 	c()
		cover95.exact.simulations0				=	c()
		cover50.exact.simulations0				=	c()
		#
		for (k in 1:p) { # k = 1;	# For testing
			theta.with.weights.RejectionABC	=	cbind(theta.table.RejectionABC[,k], rep(1/length(theta.table.RejectionABC[,k]),length(theta.table.RejectionABC[,k])))
			# above using rep(1/length(theta.table.RejectionABC[,k]),length(theta.table.RejectionABC[,k])) in place of theta.DRFweights.table[,k]
			sorted.theta.with.weights.RejectionABC	=	theta.with.weights.RejectionABC[order(theta.with.weights.RejectionABC[,1]), ]
			if (k == 1){F	=	pgamma(sorted.theta.with.weights.RejectionABC[,1],shape = 0.5 + sum.x.Pois, rate = 0.1 + n)	}
			if (k == 2){F	=	pNormScaleMixPost(sorted.theta.with.weights.RejectionABC[,1])								}
			Fk			=	cumsum(sorted.theta.with.weights.RejectionABC[,2])
			Fk.minus.1	=	c(0, Fk[1 : (length(Fk) - 1) ] )
			Dn 			= 	max(apply(cbind(abs(Fk - F), abs(F - Fk.minus.1)), 1, max))# Kolmogorov-Smirnov statistic
			S.2			=	sum(sorted.theta.with.weights.RejectionABC[,2] ^ 2)
			KSstat		=	(1 / S.2) * Dn
			# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)
			KS.Dn.RejectionABC.simulations0		=	c(KS.Dn.RejectionABC.simulations0, 	Dn		)
			KS.stat.RejectionABC.simulations0		=	c(KS.stat.RejectionABC.simulations0, 	KSstat	)
			#
			# Check for 95%(50%) coverage of RejectionABC samples:
			Q.025		=	quantile(theta.table.RejectionABC[,k], .025)
			Q.975		=	quantile(theta.table.RejectionABC[,k], .975)
			Q.25		=	quantile(theta.table.RejectionABC[,k], .25)
			Q.75		=	quantile(theta.table.RejectionABC[,k], .75)
			cover95.RejectionABC.simulations0 	= 	c(cover95.RejectionABC.simulations0, as.numeric((Q.025 <= truth[k]) & (Q.975 >= truth[k])) )
			cover50.RejectionABC.simulations0 	= 	c(cover50.RejectionABC.simulations0, as.numeric((Q.25 <= truth[k])  & (Q.75 >= truth[k])) )
			#
			# Check for 95%(50%) coverage of exact posterior distribution:
			if (k == 1){
				cover95.exact.simulations0 	= 	c(cover95.exact.simulations0, as.numeric((qgamma(.025, shape = 0.5 + sum.x.Pois, rate = 0.1 + n) <= truth[1]) & (qgamma(.975, shape = 0.5 + sum.x.Pois, rate = 0.1 + n) >= truth[1])))
				cover50.exact.simulations0 	= 	c(cover50.exact.simulations0, as.numeric((qgamma(.25,  shape = 0.5 + sum.x.Pois, rate = 0.1 + n) <= truth[1]) & (qgamma( .75, shape = 0.5 + sum.x.Pois, rate = 0.1 + n) >= truth[1])))
			}
			if (k == 2){
				# Find quantile using linear interpolation (or extrapolation, as appropriate).
				# https://en.wikipedia.org/wiki/Extrapolation
				q		=	.025
				#p1		=	pNormScaleMixPost(theta.table.RejectionABC[,2])
				mu.grid=	seq(from = -10, to = 10, by = .01)# Prior is mu ~ Uniform(-10, 10)
				p1		=	pNormScaleMixPost(mu.grid)
				x1 		=	p1[which.min(abs(p1 - q))]
				x2		=	p1[(p1 != x1)]
				x2 		=	x2[which.min(abs(x2 - q))]
				#y1		=	theta.table.RejectionABC[p1 == x1, 2]
				#y2		=	theta.table.RejectionABC[p1 == x2, 2]
				y1		=	mu.grid[p1 == x1]
				y2		=	mu.grid[p1 == x2]
				q.lo	=	y1 + ((q - x1)/(x2 - x1))*(y2 - y1)
				# Find quantile using linear interpolation (or extrapolation, as appropriate).
				q		=	.975
				#p1		=	pNormScaleMixPost(theta.table.RejectionABC[,2])
				x1 		=	p1[which.min(abs(p1 - q))]
				x2		=	p1[(p1 != x1)]
				x2 		=	x2[which.min(abs(x2 - q))]
				#y1		=	theta.table.RejectionABC[p1 == x1, 2]
				#y2		=	theta.table.RejectionABC[p1 == x2, 2]
				y1		=	mu.grid[p1 == x1]
				y2		=	mu.grid[p1 == x2]
				q.hi	=	y1 + ((q - x1)/(x2 - x1))*(y2 - y1)
				cover95.exact.simulations0 	= 	c(cover95.exact.simulations0, as.numeric(q.lo <= truth[2] & q.hi >= truth[2]))
				# Find quantile using linear interpolation (or extrapolation, as appropriate).
				# https://en.wikipedia.org/wiki/Extrapolation
				q		=	.25
				#p1		=	pNormScaleMixPost(theta.table.RejectionABC[,2]) 
				x1 		=	p1[which.min(abs(p1 - q))]
				x2		=	p1[(p1 != x1)]
				x2 		=	x2[which.min(abs(x2 - q))]
				#y1		=	theta.table.RejectionABC[p1 == x1, 2]
				#y2		=	theta.table.RejectionABC[p1 == x2, 2]
				y1		=	mu.grid[p1 == x1]
				y2		=	mu.grid[p1 == x2]
				q.lo	=	y1 + ((q - x1)/(x2 - x1))*(y2 - y1)
				# Find quantile using linear interpoation (or extrapolation, as appropriate).
				q		=	.75
				#p1		=	pNormScaleMixPost(theta.table.RejectionABC[,2])
				x1 		=	p1[which.min(abs(p1 - q))]
				x2		=	p1[(p1 != x1)]
				x2 		=	x2[which.min(abs(x2 - q))]
				#y1		=	theta.table.RejectionABC[p1 == x1, 2]
				#y2		=	theta.table.RejectionABC[p1 == x2, 2]
				y1		=	mu.grid[p1 == x1]
				y2		=	mu.grid[p1 == x2]
				q.hi	=	y1 + ((q - x1)/(x2 - x1))*(y2 - y1)
				cover50.exact.simulations0 	= 	c(cover50.exact.simulations0, as.numeric(q.lo <= truth[2] & q.hi >= truth[2]))
			}
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
			postEtheta.RejectionABC.simulations.mean0		=	apply(postEtheta.RejectionABC.simulations, 2, mean)
			postEtheta.RejectionABC.simulations.sd0		=	apply(postEtheta.RejectionABC.simulations, 2, sd)
			postMEDtheta.RejectionABC.simulations.mean0	=	apply(postMEDtheta.RejectionABC.simulations, 2, mean)
			postMEDtheta.RejectionABC.simulations.sd0		=	apply(postMEDtheta.RejectionABC.simulations, 2, sd)
			postMODEtheta.RejectionABC.simulations.mean0	=	apply(postMODEtheta.RejectionABC.simulations, 2, mean)
			postMODEtheta.RejectionABC.simulations.sd0	=	apply(postMODEtheta.RejectionABC.simulations, 2, sd)
			MLEtheta.RejectionABC.simulations.mean0		=	apply(MLEtheta.RejectionABC.simulations, 2, mean)
			MLEtheta.RejectionABC.simulations.sd0			=	apply(MLEtheta.RejectionABC.simulations, 2, sd)
			postSDtheta.RejectionABC.simulations.mean0	=	apply(postSDtheta.RejectionABC.simulations, 2, mean)
			postSDtheta.RejectionABC.simulations.sd0		=	apply(postSDtheta.RejectionABC.simulations, 2, sd)
			cover95.exact.simulations.mean0				=	apply(cover95.exact.simulations, 2, mean)
			cover50.exact.simulations.mean0				=	apply(cover50.exact.simulations, 2, mean)						
			cover95.RejectionABC.simulations.mean0		=	apply(cover95.RejectionABC.simulations, 2, mean)
			cover50.RejectionABC.simulations.mean0		=	apply(cover50.RejectionABC.simulations, 2, mean)
			KS.Dn.RejectionABC.simulations.mean0			=	apply(KS.Dn.RejectionABC.simulations, 2, mean)
			KS.Dn.RejectionABC.simulations.sd0				=	apply(KS.Dn.RejectionABC.simulations, 2, sd)
			KS.stat.RejectionABC.simulations.mean0		=	apply(KS.stat.RejectionABC.simulations, 2, mean)
			KS.stat.RejectionABC.simulations.sd0			=	apply(KS.stat.RejectionABC.simulations, 2, sd)
			MAE.postEtheta.RejectionABC.simulations0		=	apply(L1error.postEtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
			MSE.postEtheta.RejectionABC.simulations0		=	apply(L2error.postEtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
			MAE.postMEDtheta.RejectionABC.simulations0	=	apply(L1error.postMEDtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
			MSE.postMEDtheta.RejectionABC.simulations0	=	apply(L2error.postMEDtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
			MAE.postMODEtheta.RejectionABC.simulations0	=	apply(L1error.postMODEtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
			MSE.postMODEtheta.RejectionABC.simulations0	=	apply(L2error.postMODEtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
			MAE.MLEtheta.RejectionABC.simulations0		=	apply(L1error.MLEtheta.RejectionABC.simulations, 2, mean) # Mean Absolute Error (MAE)
			MSE.MLEtheta.RejectionABC.simulations0		=	apply(L2error.MLEtheta.RejectionABC.simulations, 2, mean) # Mean Squared Error (MSE)
		}
	}
	postEtheta.RejectionABC.simulations.mean[simStudy,]		=	round(postEtheta.RejectionABC.simulations.mean0, 2)
	postEtheta.RejectionABC.simulations.sd[simStudy,]		=	round(postEtheta.RejectionABC.simulations.sd0, 2)
	postMEDtheta.RejectionABC.simulations.mean[simStudy,]	=	round(postMEDtheta.RejectionABC.simulations.mean0, 2)
	postMEDtheta.RejectionABC.simulations.sd[simStudy,]		=	round(postMEDtheta.RejectionABC.simulations.sd0, 2)
	postMODEtheta.RejectionABC.simulations.mean[simStudy,]	=	round(postMODEtheta.RejectionABC.simulations.mean0, 2)
	postMODEtheta.RejectionABC.simulations.sd[simStudy,]	=	round(postMODEtheta.RejectionABC.simulations.sd0, 2)
	MLEtheta.RejectionABC.simulations.mean[simStudy,]		=	round(MLEtheta.RejectionABC.simulations.mean0, 2)
	MLEtheta.RejectionABC.simulations.sd[simStudy,]			=	round(MLEtheta.RejectionABC.simulations.sd0, 2)
	postSDtheta.RejectionABC.simulations.mean[simStudy,]	=	round(postSDtheta.RejectionABC.simulations.mean0, 2)
	postSDtheta.RejectionABC.simulations.sd[simStudy,]		=	round(postSDtheta.RejectionABC.simulations.sd0, 2)
	cover95.RejectionABC.simulations.mean[simStudy,]		=	round(cover95.RejectionABC.simulations.mean0, 2)
	cover50.RejectionABC.simulations.mean[simStudy,]		=	round(cover50.RejectionABC.simulations.mean0, 2)
	cover95.exact.simulations.mean[simStudy,]				=	round(cover95.exact.simulations.mean0, 2)
	cover50.exact.simulations.mean[simStudy,]				=	round(cover50.exact.simulations.mean0, 2)	
	KS.Dn.RejectionABC.simulations.mean[simStudy,]			=	round(KS.Dn.RejectionABC.simulations.mean0, 2)
	KS.Dn.RejectionABC.simulations.sd[simStudy,]				=	round(KS.Dn.RejectionABC.simulations.sd0, 2)
	KS.stat.RejectionABC.simulations.mean[simStudy,]		=	round(KS.stat.RejectionABC.simulations.mean0, 2)
	KS.stat.RejectionABC.simulations.sd[simStudy,]			=	round(KS.stat.RejectionABC.simulations.sd0, 2)
	MAE.postEtheta.RejectionABC.simulations[simStudy,]		= 	round(MAE.postEtheta.RejectionABC.simulations0, 2)
	MSE.postEtheta.RejectionABC.simulations[simStudy,]		= 	round(MSE.postEtheta.RejectionABC.simulations0, 2)
	MAE.postMEDtheta.RejectionABC.simulations[simStudy,]	= 	round(MAE.postMEDtheta.RejectionABC.simulations0, 2)
	MSE.postMEDtheta.RejectionABC.simulations[simStudy,]	= 	round(MSE.postMEDtheta.RejectionABC.simulations0, 2)
	MAE.postMODEtheta.RejectionABC.simulations[simStudy,]	= 	round(MAE.postMODEtheta.RejectionABC.simulations0, 2)
	MSE.postMODEtheta.RejectionABC.simulations[simStudy,]	= 	round(MSE.postMODEtheta.RejectionABC.simulations0, 2)
	MAE.MLEtheta.RejectionABC.simulations[simStudy,]		= 	round(MAE.MLEtheta.RejectionABC.simulations0, 2)
	MSE.MLEtheta.RejectionABC.simulations[simStudy,]		= 	round(MSE.MLEtheta.RejectionABC.simulations0, 2)
}


postEtheta.RejectionABC.simulations.mean.sd 	= 	matrix(paste(postEtheta.RejectionABC.simulations.mean, ' (', postEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
postMEDtheta.RejectionABC.simulations.mean.sd 	= 	matrix(paste(postMEDtheta.RejectionABC.simulations.mean, ' (', postMEDtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
postMODEtheta.RejectionABC.simulations.mean.sd 	= 	matrix(paste(postMODEtheta.RejectionABC.simulations.mean, ' (', postMODEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
MLEtheta.RejectionABC.simulations.mean.sd 		= 	matrix(paste(MLEtheta.RejectionABC.simulations.mean, ' (', MLEtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
postSDtheta.RejectionABC.simulations.mean.sd 	=	matrix(paste(postSDtheta.RejectionABC.simulations.mean, '(', postSDtheta.RejectionABC.simulations.sd, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
cover95.50.RejectionABC.simulations.mean 		= 	matrix(paste(cover95.RejectionABC.simulations.mean, ' (', cover50.RejectionABC.simulations.mean, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
cover95.50.exact.simulations.mean 				= 	matrix(paste(cover95.exact.simulations.mean, ' (', cover50.exact.simulations.mean, ')', sep=''),nrow=nrow(postEtheta.RejectionABC.simulations.mean), ncol=ncol(postEtheta.RejectionABC.simulations.mean))
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
rownames(postSDtheta.RejectionABC.simulations.mean.sd) 	=	simConditionNames
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
colnames(postMEDtheta.RejectionABC.simulations.mean.sd)	=	parameterNames
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
colnames(MAE.MSE.postMEDtheta.RejectionABC.simulations)	=	parameterNames
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