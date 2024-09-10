#install.packages('igraph');# Using igraph version before January 30, 2024, for more reliable use of sample_pa()
#library('igraph')		# For network data analysis/management
#install.packages('ergm'); 			library('ergm')		# For network data analysis/management
#install.packages('ergm.count'); 	library('ergm.count')# For network data analysis/management
#install.packages('ergm.rank'); 	library('ergm.rank')	# For network data analysis/management
install.packages('drf');			library('drf')		# Distributional Random Forests.
install.packages('nvmix');			library('nvmix')		# For t copula.
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('ddpcr'); 			library('ddpcr')		# For quiet() 
#install.packages('extraDistr');	library('extraDistr')# 	For rtbinom() truncated binomial sampling.
install.packages('goftest');		library('goftest')	# For pAD() and qAD()
#setwd('C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork')
#setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')
#setwd('C:\\Users\\George Karabatsos\\Desktop\\BA\\code ABCnetwork')
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork')

rm(list = ls())
# source('network statistics 9-14-23.R')
# source('simulateNetwork 2-9-24.R')
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Run simulation over replicate datasets:
# ====================================================================================================================================
replicas		=	10 # 
for (r in 1 : replicas) {	# r = 1
	# ====================================================================================================================================
	# Simulate network dataset from Poisson(lambda=3) model and normal scale mixture model.
	# ====================================================================================================================================
	parameterNames	=	c('lambda', 'mu')
	d					=	length(parameterNames)
	lambda 	= 	3;	# Mean (and variance) parameter of Poisson(lambda) distribution	
	mu			=	0	# Mean parameter of (1/2)*N(mu,var=1) + (1/2)*N(mu,var=0.01) normal mixture.
	truth 		= 	c(lambda, mu); # True data-generating model parameters to be estimated. 
	n				=	100	# Sataset sample size ( to be simulated)
	#
	# Poisson model data simulation:
	x 				= 	rpois(n, lambda = 3)
	sum.x.Pois		=	sum(x)
	mean.x.Pois	=	mean(x)
	# For testing:
	#N				=	10000
	#thetaSamples	=	rgamma(N, shape = 0.5 + sum(x), rate = 0.1 + n)
	#Out.ad.test 	=	ad.test(thetaSamples, "pgamma", shape = 0.5 + sum(x), rate = 0.1 + n, estimated = FALSE)
	#cdf.AD 		=	pAD(q = Out.ad.test$statistic, n = N, lower.tail = TRUE, fast = TRUE)# CDF
	#pvalue.AD		=	1 - cdf.AD		# They equal:  c( pvalue.AD, Out.ad.test$p.value )
	#q.AD 			=	qAD(p = .95, n = N, lower.tail = TRUE, fast = TRUE) # equals 2.492229,critical value for alpha = .05
	# 1 - pAD(q = q.AD, n = N, lower.tail = TRUE, fast = TRUE) # equals 0.05000141
	PostMean.lambda	=	(0.5 + sum(x)) / (0.1 + n) # Posterior mean
	PostMed.lambda	=	qgamma(0.5, shape = 0.5 + sum(x), rate = 0.1 + n)
	PostMode.lambda	=	((0.5 + sum(x)) - 1) / (0.1 + n) # Posterior mode
	MLE.lambda			= 	mean.x.Pois	# Sufficient statistic for Poisson model.
	SD.lambda		=	sqrt( (0.5 + sum(x)) / (( 0.1 + n) ^ 2) )
	#
	# Normal scale mixture model data simulation.
	z 				=	rbinom(n, size = 1, prob = 1 / 2)
	x 				= 	z * rnorm(n, mean = 0, sd = sqrt(1)) + (1 - z) * rnorm(n, mean = 0, sd = sqrt(0.01))
	mean.x.NorMix	=	mean(x)
	# For testing:
	#N				=	10000
	#z.Post 		=	rbinom(n = N, size = 1, prob = 1 / 2)
	#thetaSamples 	=	z.Post * rnorm(N, mean = mean.x, sd = sqrt(1)) + (1 - z.Post) * rnorm(N, mean = mean.x, sd = sqrt(0.01)) 
	#pNormScaleMixPost	=	function(q){# Posterior CDF of mu (common mean parameter of normal mixture).
	#							Out.pNormScaleMixPost = (1/2)*pnorm(q, mean = mean.x.NorMix, sd = sqrt(1)) + (1/2)*pnorm(q, mean = mean.x.NorMix, sd = sqrt(0.01))
	#							return(Out.pNormScaleMixPost)}
	#Out.ad.test 	=	ad.test(thetaSamples, "pNormScaleMixPost", estimated = FALSE)
	#cdf.AD 		=	pAD(q = Out.ad.test$statistic, n = N, lower.tail = TRUE, fast = TRUE)# CDF
	#pvalue.AD		=	1 - cdf.AD		# They equal:  c( pvalue.AD, Out.ad.test$p.value )
	#q.AD 			=	qAD(p = .95, n = N, lower.tail = TRUE, fast = TRUE) # equals 2.492229,critical value for alpha = .05
	# 1 - pAD(q = q.AD, n = N, lower.tail = TRUE, fast = TRUE) # equals 0.05000141
	PostMean.mu	=	mean.x.NorMix # Posterior mean
	PostMed.mu		=	mean.x.NorMix # Posterior median
	PostMode.mu	=	mean.x.NorMix # Posterior mode
	MLE.mu			= 	mean.x.NorMix
	SD.mu			=	sqrt( ((1 / 2) * (1 + (PostMean.mu ^2))) + ((1 / 2) * (0.01 + (PostMean.mu ^2) )) - (PostMean.mu ^ 2) )
	# https://en.wikipedia.org/wiki/Mixture_distribution
	#
	# ====================================================================================================================================
	# Summary statistics from data x 
	# ====================================================================================================================================
	sTerms		=	c('MnPois', 'MnNorScMix')
	sx 			=	c( mean.x.Pois, mean.x.NorMix) 
	names(sx)	=	sTerms
	p 			=	length(sTerms)
	postEtheta.exact	=	c(PostMean.lambda, PostMean.mu)
	postMEDtheta.exact	=	c(PostMed.lambda, PostMed.mu)
	posteriorMode.exact	=	c(PostMode.lambda, PostMode.mu)
	MLE.exact			=	c(MLE.lambda, MLE.mu)
	postSDtheta.exact	=	c(SD.lambda, SD.mu)
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for Price model
	# ====================================================================================================================================
	N				=	10000 #	20		# Sample size of the reference table
	n.table 		= 	100;
	theta.table	=	matrix(NA, nrow = N, ncol = d)
	sy.table		=	matrix(NA, nrow = N, ncol = length(sTerms))
	colnames(theta.table)	<-	parameterNames
	colnames(sy.table)		<-	sTerms
	# 
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j=	c(rgamma(1, shape = 1/2, rate = 0.1), runif(1, min = -10, max = 10))
		#
		# Poisson model:
		y 		= 	rpois(n.table, lambda = theta.j[1])
		mean.y.Pois	=	mean(y)
		#
		# Normal mixture model:
		z 		=	rbinom(n.table, size = 1, prob = 1 / 2)
		y 		= 	z * rnorm(n.table, mean = theta.j[2], sd = sqrt(1)) + (1 - z) * rnorm(n.table, mean = theta.j[2], sd = sqrt(0.01))
		mean.y.NorMix	=	mean(y)
		#
		sy		=	c(mean.y.Pois, mean.y.NorMix)
		names(sy)	<-	sTerms
		sy.j		= 	sy
		#
		# Update Reference Table
		theta.table[j,]	=	theta.j
		sy.table[j,]		=	sy.j
		if ((j / 10) == round(j / 10))	{	flush.console()	}
	}
	# Restrict to subset of theta samples corresponding to acceptable 
	isFinite.sy		=	apply(is.finite(sy.table),1,all)
	notNA.sy			=	apply(! is.na(sy.table),1,all)
	isOK.sy			=	isFinite.sy & notNA.sy
	theta.table 		= 	theta.table[isOK.sy, ]
	sy.table			=	sy.table[isOK.sy, ]
	N 					=	sum(isOK.sy)
	if (d == 1){theta.table = as.matrix(theta.table)}
	if (length(sTerms) == 1){sy.table = as.matrix(sy.table)}
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Construct multivariate meta-t copula ABC posterior distribution based on DRF marginals.
	# ====================================================================================================================================
	# For each parameter (dependent variable), and set of predictors (sy's), 
	# Train a distributional random forest with CART splitting rule.
	postEtheta		=	matrix(NA,	nrow = 1, ncol = d)
	postSDtheta	=	matrix(NA,	nrow = 1, ncol = d)
	postQtheta		=	matrix(NA,	nrow = 5, ncol = d)
	theta.DRFweights.table	=	matrix(NA,	nrow = N, ncol = d)
	breaks			=	matrix(NA,	nrow = N, ncol = d)
	u				=	matrix(NA,	nrow = N, ncol = d)
	pi.thetaRT		=	matrix(NA,	nrow = N, ncol = d)# To store marginal posterior PDFs values of unordered theta.jk's from Reference Table 
	PI.theta		=	matrix(NA,	nrow = N, ncol = d)# To store the d marginal posterior CDFs of ordered theta.jk's
	pi.theta		=	matrix(NA,	nrow = N, ncol = d)# To store the d marginal posterior PDFs of ordered theta.jk's
	colnames(postEtheta)		<-	parameterNames
	colnames(postSDtheta)	<-	parameterNames
	colnames(postQtheta)		<-	parameterNames
	colnames(u)				<-	parameterNames
	colnames(PI.theta	)		<-	parameterNames
	colnames(pi.theta	)		<-	parameterNames
	colnames(pi.thetaRT)		<-	parameterNames
	rownames(postQtheta)		<-	c( '2.5%', '25%', '50%', '75%', '97.5%')
	for (k in 1 : d) {# k = 1
		# Train Distributional Random Forest (DRF) on the theta.jk's: 
		drf.forest 		=	drf(X = sy.table, Y = theta.table[,k])#, compute.variable.importance = TRUE # Costly to compute
		X.test				=	sx
		# Use trained DRF to estimate (predict) posterior mean, standard deviation, and quantiles, conditional on sx:
		postEtheta.k		=	predict(drf.forest, newdata	= X.test, functional = "mean")
		postSDtheta.k		=	predict(drf.forest, newdata	= X.test, functional = "sd")
		postQtheta.k		=	predict(drf.forest, newdata	= X.test, functional = "quantile", quantiles = c(.025,.25,.5,.75,.975))
		postQtheta.k		=	c(postQtheta.k$quantile)
		postEtheta[1,k]	=	c(postEtheta.k$mean)
		postSDtheta[1,k]	=	c(postSDtheta.k$sd)
		postQtheta[,k]	=	c(postQtheta.k)
		# From the trained DRF, extract the theta.jk's and corresponding weights used for predictions:
		pred				=	predict(drf.forest, newdata	= X.test)
		theta.k.vals		=	pred$y	# theta.k.vals - theta.table[,k] equals zeros (as desired)
		theta.k.weights	=	pred$weights[1,]
		theta.DRFweights.table[,k]	=	theta.k.weights	
		# Marginal posterior CDFs values of sampled ordered theta.j's:
		ord					=	order(theta.k.vals)# returns a permutation which rearranges its input into ascending order.
		breaks[,k]			=	theta.k.vals[ord]
		PI.theta[,k]		=	cumsum(theta.k.weights[ord])	# marginal posterior CDF (for ordered theta.jk's)
		# Marginal posterior PDFs values of sampled ordered theta.j's:
		binwidths			=	breaks[2:N, k] - breaks[1:(N-1), k]
		binwidths			=	c(binwidths[1], binwidths)# histogram bins are right-closed (left open) intervals.
		pi.theta[,k]		=	theta.k.weights[ord]  / binwidths	# marginal posterior PDF (for ordered theta.jk's)
		# 'Empirical density is equal to the empirical probability divided by the interval length, or binwidth.'
		# https://ocw.tudelft.nl/course-readings/pre-1-3-histogram-probability-density-function/
		# https://stackoverflow.com/questions/74125399/how-to-calculate-the-density-in-results-of-function-hist-in-r
		# Posterior CDF values (u[,k]) and PDF values for each theta.j,k (in original order of sampled theta.j,k's in Reference Table):
		u[,k]				=	PI.theta[order(ord),k]# breaks[order(ord),k] - theta.k.vals # all zeros, as desired (see link below)
		# https://stackoverflow.com/questions/15464793/restoring-original-order-of-a-vector-matrix-in-r
		# Posterior PDF values for each theta.j,k (in original order of the sampled theta.j,k's in Reference Table):
		pi.thetaRT[,k]	=	pi.theta[order(ord),k] # OK (checked)
	}
	# Remove rows for which u_jk = 0 or 1. (they correspond to theta_j's with zero posterior density anyway).
	keep				=	apply((u > 0) & (u < 1), 1, all)
	u.keep				=	u[keep,]
	theta.table.keep	=	theta.table[keep, ]
	pi.thetaRT.keep 	= 	pi.thetaRT[keep,]
	# Estimate copula parameters (df, correlation matrix) of Meta-t posterior distribution:
	OutCopulaFit 		=	fitStudentcopula(u[keep, ], fit.method = "EM-MLE", df.bounds = c(0.1, 1000), verbose = TRUE)
	post.df			=	OutCopulaFit$df
	post.Correlations	=	OutCopulaFit$scale
	# Meta-t posterior densities of (retained) unordered theta.j's from Reference table
	posteriorPDFs		=	dStudentcopula(u.keep,df=post.df,scale=post.Correlations)*apply(pi.thetaRT.keep,1,prod)
	# Find posterior mode and MLE of theta using the sampled values:
	posteriorMode		=	theta.table.keep[posteriorPDFs == max(posteriorPDFs, na.rm = T),]
	priorPDFs			=	apply(cbind(dgamma(theta.table.keep[,1], shape = 1/2, rate = 0.1), dunif(theta.table.keep[,2], min = -10, max = 10)), 1, prod)
	likelihoods		=	posteriorPDFs / priorPDFs
	MLE					=	theta.table.keep[likelihoods == max(likelihoods, na.rm = T),]
	# ------------------------------------------------------------------------------------------------------------------------------------
	#
	SummaryTable		=	rbind(postEtheta, postSDtheta, postQtheta, posteriorMode, MLE, rep(post.df, d), post.Correlations)
	rownames(SummaryTable)[1:2]		=	c('PostMean', 'PostSD')
	rownames(SummaryTable)[nrow(SummaryTable) - d]		=	'post.df'
	rownames(SummaryTable)[(nrow(SummaryTable)- d + 1) : nrow(SummaryTable)]	=	paste('Cor.', parameterNames, sep='')
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Save and update results:
	# ====================================================================================================================================
	end_time 								= 	Sys.time()
	TotalComputationTimeSimulation			=	end_time - start_time
	TotalComputationTimeAllSimulations		=	TotalComputationTimeAllSimulations + TotalComputationTimeSimulation

	# Initiate collections of simulation results:
	if (r == 1){	# r = 1
		postEtheta.simulations				=	c()
		postEtheta.simulations.exact		=	c()

		postMEDtheta.simulations			=	c()
		postMEDtheta.simulations.exact		=	c()

		posteriorMode.simulations			=	c()
		posteriorMode.simulations.exact	=	c()

		MLE.simulations						=	c()
		MLE.simulations.exact				=	c()

		KS.Dn.simulations						=	c()
		KS.stat.simulations					=	c()

		L1error.postEtheta.simulations			=	c()
		L1error.postEtheta.simulations.exact	=	c()

		L2error.postEtheta.simulations			=	c()
		L2error.postEtheta.simulations.exact	=	c()

		L1error.postMEDtheta.simulations			=	c()
		L1error.postMEDtheta.simulations.exact	=	c()

		L2error.postMEDtheta.simulations			=	c()
		L2error.postMEDtheta.simulations.exact	=	c()

		L1error.posteriorMode.simulations			=	c()
		L1error.posteriorMode.simulations.exact	=	c()

		L2error.posteriorMode.simulations			=	c()
		L2error.posteriorMode.simulations.exact	=	c()

		L1error.MLE.simulations				=	c()
		L1error.MLE.simulations.exact		=	c()

		L2error.MLE.simulations				=	c()
		L2error.MLE.simulations.exact		=	c()

		postSDtheta.simulations				=	c()
		postSDtheta.simulations.exact		=	c()

		post95cover.simulations				=	c()
		post95cover.simulations.exact		=	c()
		postIQRcover.simulations			=	c()
		postIQRcover.simulations.exact		=	c()

		post.df.simulations					=	c()
		post.Correlations.simulations		=	c()

		TotalComputationTimePerSimulation	=	c()
	}

	# Update collections of simulation results:
	postEtheta.simulations					=	rbind(postEtheta.simulations, 				postEtheta					)
	postEtheta.simulations.exact			=	rbind(postEtheta.simulations.exact, 		postEtheta.exact			)

	postMEDtheta.simulations				=	rbind(postMEDtheta.simulations, 			postQtheta[3,]				)
	postMEDtheta.simulations.exact			=	rbind(postMEDtheta.simulations.exact, 		postMEDtheta.exact			)

	posteriorMode.simulations				=	rbind(posteriorMode.simulations, 			posteriorMode				)
	posteriorMode.simulations.exact		=	rbind(posteriorMode.simulations.exact, 		posteriorMode.exact			)

	MLE.simulations							=	rbind(MLE.simulations, 						MLE							)
	MLE.simulations.exact					=	rbind(MLE.simulations.exact, 				MLE.exact					)


	# Do weighted one-sample (two-tailed) Kolmogorov Smirnov test (Monohan,2011,p.358)
	theta.with.weights			=	cbind(theta.table[,1], theta.DRFweights.table[,1])
	sorted.theta.with.weights	=	theta.with.weights[order(theta.with.weights[,1]), ]
	F								=	pgamma(sorted.theta.with.weights[,1],shape = 0.5 + sum.x.Pois, rate = 0.1 + n)
	Fk								=	cumsum(sorted.theta.with.weights[,2])
	Fk.minus.1						=	c(0, Fk[1 : (length(Fk) - 1) ] )
	Dn 								= 	max(apply(cbind(abs(Fk - F), abs(F - Fk.minus.1)), 1, max))# Kolmogorov-Smirnov statistic
	S.2								=	sum(sorted.theta.with.weights[,2] ^ 2)
	Dn1								=	Dn
	KSstat1						=	(1 / S.2) * Dn
	# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)
	#
	# Do weighted one-sample (two-tailed) Kolmogorov Smirnov test (Monohan,2011,p.358)
	theta.with.weights			=	cbind(theta.table[,2], theta.DRFweights.table[,2])
	sorted.theta.with.weights	=	theta.with.weights[order(theta.with.weights[,1]), ]
	pNormScaleMixPost	=	function(q){# Posterior CDF of mu (common mean parameter of normal mixture).
							Out.pNormScaleMixPost = (1/2)*pnorm(q, mean = mean.x.NorMix, sd = sqrt(1)) + (1/2)*pnorm(q, mean = mean.x.NorMix, sd = sqrt(0.01))
	return(Out.pNormScaleMixPost)}
	F								=	pNormScaleMixPost(sorted.theta.with.weights[,1])
	Fk								=	cumsum(sorted.theta.with.weights[,2])
	Fk.minus.1						=	c(0, Fk[1 : (length(Fk) - 1) ] )
	Dn 								= 	max(apply(cbind(abs(Fk - F), abs(F - Fk.minus.1)), 1, max))# Kolmogorov-Smirnov statistic
	S.2								=	sum(sorted.theta.with.weights[,2] ^ 2)
	Dn2								=	Dn
	KSstat2						=	(1 / S.2) * Dn
	# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)
	#
	KS.Dn.simulations				=	rbind(KS.Dn.simulations, 	c(Dn1, Dn2)			)
	KS.stat.simulations			=	rbind(KS.stat.simulations, 	c(KSstat1, KSstat2)	)
	# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)

	L1error.postEtheta.simulations			=	rbind(L1error.postEtheta.simulations, 			abs(postEtheta - truth) 		)
	L1error.postEtheta.simulations.exact	=	rbind(L1error.postEtheta.simulations.exact, 	abs(postEtheta.exact - truth) 		)

	L2error.postEtheta.simulations			=	rbind(L2error.postEtheta.simulations,    		(postEtheta - truth)^2 				)
	L2error.postEtheta.simulations.exact	=	rbind(L2error.postEtheta.simulations.exact, 	(postEtheta.exact - truth)^2 		)

	L1error.postMEDtheta.simulations			=	rbind(L1error.postMEDtheta.simulations,			abs(postQtheta[3,] - truth)		)
	L1error.postMEDtheta.simulations.exact	=	rbind(L1error.postMEDtheta.simulations.exact,	abs(postMEDtheta.exact - truth))

	L2error.postMEDtheta.simulations			=	rbind(L2error.postMEDtheta.simulations,  		(postQtheta[3,] - truth)^2 		)
	L2error.postMEDtheta.simulations.exact	=	rbind(L2error.postMEDtheta.simulations.exact,  	(postMEDtheta.exact - truth)^2	)

	L1error.posteriorMode.simulations			=	rbind(L1error.posteriorMode.simulations, 		abs(posteriorMode - truth) 		)
	L1error.posteriorMode.simulations.exact	=	rbind(L1error.posteriorMode.simulations.exact, 	abs(posteriorMode.exact - truth))

	L2error.posteriorMode.simulations			=	rbind(L2error.posteriorMode.simulations,    	(posteriorMode - truth)^2 		)
	L2error.posteriorMode.simulations.exact	=	rbind(L2error.posteriorMode.simulations.exact,  (posteriorMode.exact - truth)^2 )

	L1error.MLE.simulations					=	rbind(L1error.MLE.simulations, 						abs(MLE - truth) 				)
	L1error.MLE.simulations.exact			=	rbind(L1error.MLE.simulations.exact,					abs(MLE.exact - truth) 		)

	L2error.MLE.simulations					=	rbind(L2error.MLE.simulations,    						(MLE - truth)^2 				)
	L2error.MLE.simulations.exact			=	rbind(L2error.MLE.simulations.exact,					(MLE.exact - truth)^2 		)

	postSDtheta.simulations					=	rbind(postSDtheta.simulations, 						postSDtheta					)
	postSDtheta.simulations.exact			=	rbind(postSDtheta.simulations.exact, 					postSDtheta.exact				)
	
	post95cover.simulations					=	rbind(post95cover.simulations,	ifelse((truth > postQtheta[1,]) & (truth < postQtheta[5,]),1,0))
	postIQRcover.simulations				=	rbind(postIQRcover.simulations,ifelse((truth > postQtheta[2,]) & (truth < postQtheta[4,]),1,0))
	post.df.simulations						=	rbind(post.df.simulations, rep(post.df, d))
	post.Correlations.simulations			=	rbind(post.Correlations.simulations, c(post.Correlations))
	TotalComputationTimePerSimulation		=	rbind(TotalComputationTimePerSimulation, rep(TotalComputationTimeSimulation, d))


	if (r == replicas){
		Mean.postEtheta.simulations			=	apply(postEtheta.simulations, 2, mean);				SD.postEtheta.simulations			=	apply(postEtheta.simulations, 2, sd)
		Mean.postEtheta.simulations.exact	=	apply(postEtheta.simulations.exact, 2, mean);		SD.postEtheta.simulations.exact	=	apply(postEtheta.simulations.exact, 2, sd)

		Mean.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, mean);				SD.postMEDtheta.simulations			=	apply(postMEDtheta.simulations, 2, sd)
		Mean.postMEDtheta.simulations.exact	=	apply(postMEDtheta.simulations.exact, 2, mean);	SD.postMEDtheta.simulations.exact	=	apply(postMEDtheta.simulations.exact, 2, sd)

		Mean.posteriorMode.simulations		=	apply(posteriorMode.simulations, 2, mean);			SD.posteriorMode.simulations		=	apply(posteriorMode.simulations, 2, sd)
		Mean.posteriorMode.simulations.exact=	apply(posteriorMode.simulations.exact, 2, mean);	SD.posteriorMode.simulations.exact	=	apply(posteriorMode.simulations.exact, 2, sd)

		Mean.MLE.simulations					=	apply(MLE.simulations, 2, mean);						SD.MLE.simulations					=	apply(MLE.simulations, 2, sd)
		Mean.MLE.simulations.exact			=	apply(MLE.simulations.exact, 2, mean);				SD.MLE.simulations.exact			=	apply(MLE.simulations.exact, 2, sd)

		Mean.KS.Dn.simulations				=	apply(KS.Dn.simulations, 2, mean);						SD.KS.Dn.simulations					=	apply(KS.Dn.simulations, 2, sd)
		Mean.KS.stat.simulations			=	apply(KS.stat.simulations, 2, mean);					SD.KS.stat.simulations				=	apply(KS.stat.simulations, 2, sd)

		MAE.postEtheta.simulations			=	apply(L1error.postEtheta.simulations, 2, mean) # Mean Absolute Error (MAE)
		MAE.postEtheta.simulations.exact	=	apply(L1error.postEtheta.simulations.exact, 2, mean) # Mean Absolute Error (MAE)

		MSE.postEtheta.simulations			=	apply(L2error.postEtheta.simulations, 2, mean) # Mean Squared Error (MSE)
		MSE.postEtheta.simulations.exact	=	apply(L2error.postEtheta.simulations.exact, 2, mean) # Mean Squared Error (MSE)

		MAE.postMEDtheta.simulations		=	apply(L1error.postMEDtheta.simulations, 2, mean) # Mean Absolute Error (MAE)
		MAE.postMEDtheta.simulations.exact	=	apply(L1error.postMEDtheta.simulations.exact, 2, mean) # Mean Absolute Error (MAE)

		MSE.postMEDtheta.simulations		=	apply(L2error.postMEDtheta.simulations, 2, mean) # Mean Squared Error (MSE)
		MSE.postMEDtheta.simulations.exact	=	apply(L2error.postMEDtheta.simulations.exact, 2, mean) # Mean Squared Error (MSE)

		MAE.posteriorMode.simulations		=	apply(L1error.posteriorMode.simulations, 2, mean) # Mean Absolute Error (MAE)
		MAE.posteriorMode.simulations.exact	=	apply(L1error.posteriorMode.simulations.exact, 2, mean) # Mean Absolute Error (MAE)

		MSE.posteriorMode.simulations		=	apply(L2error.posteriorMode.simulations, 2, mean) # Mean Squared Error (MSE)
		MSE.posteriorMode.simulations.exact	=	apply(L2error.posteriorMode.simulations.exact, 2, mean) # Mean Squared Error (MSE)

		MAE.MLE.simulations					=	apply(L1error.MLE.simulations, 2, mean) # Mean Absolute Error (MAE)
		MAE.MLE.simulations.exact			=	apply(L1error.MLE.simulations.exact, 2, mean) # Mean Absolute Error (MAE)

		MSE.MLE.simulations					=	apply(L2error.MLE.simulations, 2, mean) # Mean Squared Error (MSE)
		MSE.MLE.simulations.exact			=	apply(L2error.MLE.simulations.exact, 2, mean) # Mean Squared Error (MSE)

		Mean.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, mean);		SD.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, sd)
		Mean.postSDtheta.simulations.exact	=	apply(postSDtheta.simulations.exact, 2, mean);	SD.postSDtheta.simulations.exact=	apply(postSDtheta.simulations.exact, 2, sd)

		Prob.post95cover.simulations		=	apply(post95cover.simulations, 2, mean)
		Prob.postIQRcover.simulations		=	apply(postIQRcover.simulations, 2, mean)
		Mean.post.df.simulations			=	apply(post.df.simulations, 2, mean);				SD.post.df.simulations				=	apply(post.df.simulations, 2, sd)
		Mean.post.Correlations.simulations	=	apply(post.Correlations.simulations, 2, mean);		SD.post.Correlations.simulations	=	apply(post.Correlations.simulations, 2, sd)
		Mean.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, mean);	SD.ComputationTimePerSimulation		=	apply(TotalComputationTimePerSimulation, 2, sd)

		SimulationStudySummaryTable			=	rbind(	
		truth,
		Mean.postEtheta.simulations,			SD.postEtheta.simulations,
		Mean.postEtheta.simulations.exact,		SD.postEtheta.simulations.exact,

		Mean.postMEDtheta.simulations,			SD.postMEDtheta.simulations,
		Mean.postMEDtheta.simulations.exact,	SD.postMEDtheta.simulations.exact,

		Mean.posteriorMode.simulations,			SD.posteriorMode.simulations,
		Mean.posteriorMode.simulations.exact,	SD.posteriorMode.simulations.exact,

		Mean.MLE.simulations,					SD.MLE.simulations,
		Mean.MLE.simulations.exact,				SD.MLE.simulations.exact,

		Mean.KS.Dn.simulations,					SD.KS.Dn.simulations,
		Mean.KS.stat.simulations,				SD.KS.stat.simulations,

		MAE.postEtheta.simulations,
		MAE.postEtheta.simulations.exact,

		MSE.postEtheta.simulations,
		MSE.postEtheta.simulations.exact,

		MAE.postMEDtheta.simulations,
		MAE.postMEDtheta.simulations.exact,
		
		MSE.postMEDtheta.simulations,
		MSE.postMEDtheta.simulations.exact,

		MAE.posteriorMode.simulations,
		MAE.posteriorMode.simulations.exact,
		
		MSE.posteriorMode.simulations,
		MSE.posteriorMode.simulations.exact,
		
		MAE.MLE.simulations,
		MAE.MLE.simulations.exact,
		
		MSE.MLE.simulations,
		MSE.MLE.simulations.exact,
		
		Mean.postSDtheta.simulations,			SD.postSDtheta.simulations,
		Mean.postSDtheta.simulations.exact,	SD.postSDtheta.simulations.exact,

		Prob.post95cover.simulations,	
		Prob.postIQRcover.simulations,
		Mean.post.df.simulations,				SD.post.df.simulations,
		Mean.ComputationTimePerSimulation,		SD.ComputationTimePerSimulation,
		sum(TotalComputationTimeAllSimulations)								)
	
		rownames(SimulationStudySummaryTable)=	c(		
		'truth',
		'Mean.postEtheta.simulations',				'SD.postEtheta.simulations',
		'Mean.postEtheta.simulations.exact',		'SD.postEtheta.simulations.exact',

		'Mean.postMEDtheta.simulations',			'SD.postMEDtheta.simulations',
		'Mean.postMEDtheta.simulations.exact',	'SD.postMEDtheta.simulations.exact',

		'Mean.posteriorMode.simulations',			'SD.posteriorMode.simulations',
		'Mean.posteriorMode.simulations.exact',	'SD.posteriorMode.simulations.exact',

		'Mean.MLE.simulations',						'SD.MLE.simulations',
		'Mean.MLE.simulations.exact',				'SD.MLE.simulations.exact',

		'Mean.KS.Dn.simulations',					'SD.KS.Dn.simulations',
		'Mean.KS.stat.simulations',					'SD.KS.stat.simulations',

		'MAE.postEtheta.simulations',
		'MAE.postEtheta.simulations.exact',

		'MSE.postEtheta.simulations',
		'MSE.postEtheta.simulations.exact',

		'MAE.postMEDtheta.simulations',
		'MAE.postMEDtheta.simulations.exact',

		'MSE.postMEDtheta.simulations',
		'MSE.postMEDtheta.simulations.exact',

		'MAE.posteriorMode.simulations',
		'MAE.posteriorMode.simulations.exact',
		
		'MSE.posteriorMode.simulations',
		'MSE.posteriorMode.simulations.exact',

		'MAE.MLE.simulations',
		'MAE.MLE.simulations.exact',

		'MSE.MLE.simulations',
		'MSE.MLE.simulations.exact',

		'Mean.postSDtheta.simulations',				'SD.postSDtheta.simulations',
		'Mean.postSDtheta.simulations.exact',		'SD.postSDtheta.simulations.exact',

		'Prob.post95cover.simulations',	
		'Prob.postIQRcover.simulations',
		'Mean.post.df.simulations',			'SD.post.df.simulations',
		'Mean.ComputationTimePerSimulation','SD.ComputationTimePerSimulation',
		'TotalComputationTimeAllSimulations'									)													
		colnames(SimulationStudySummaryTable)=	c(parameterNames)
		SimulationStudySummaryTableCorrel	=	rbind(	Mean.post.Correlations.simulations,	SD.post.Correlations.simulations)
		rownames(SimulationStudySummaryTableCorrel)=	c('Mean.post.CorrelationMatrix.simulations', 'SD.post.CorrelationsMatrix.simulations')
	}
	outputFileName	=	paste("Poisson normScaleMix simul n = 100 ntable = ", n.table, " ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	print(paste('Replica', r, 'done'))
	flush.console()
	# ====================================================================================================================================
}