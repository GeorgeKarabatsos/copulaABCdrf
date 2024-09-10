install.packages('drf');		library('drf')	# Distributional Random Forests.
install.packages('nvmix');		library('nvmix')	# For t copula.
install.packages('mvtnorm');	library('mvtnorm')# Multivariate Normal and t Distributions
install.packages('MASS');		library('MASS')	# For mvrnorm().
install.packages('ddpcr'); 		library('ddpcr')	# For quiet() 
install.packages('kde1d');		library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork')

rm(list = ls())
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Run simulation over replicate datasets:
# ====================================================================================================================================
replicas		=	1
for (r in 1 : replicas) {	# r = 1
	# ====================================================================================================================================
	# Simulate network dataset from twisted Gaussian model (300-dimensional).
	# ====================================================================================================================================
	parameterNames	=	paste('theta',1:300,sep='')
	d		=	length(parameterNames)
	truth	= 	c(10, rep(0, 299)) # True data-generating model parameters to be estimated.
	theta	=	truth
	# Likelihood x ~ Normal(theta, I_p), with dim(theta) = 300,
	# theta is assigned a twisted normal prior density, proportional to: 
	# exp( -((theta[1]^2) / 200) - ((theta[2] - (b*(theta[1]^2)) + (100*b))^2)/2 - sum(theta[3:d]^2) )	
   	# with b = 0.1, so that (theta[1], theta[2]) are highly correlated a-priori.
	# Likelihood only provides location information about the model parameters,
	# Dependence structure about posterior distribution of theta 
	# comes mostly from the prior distribution.  
	#
	n			=	300	# Dataset sample size (to be simulated)
	#
	# Model data simulation:
	x			=	matrix(truth, ncol = 1)
	# Data summaries:
	sx			=	t(x)
	sTerms		=	paste('x', 1:d, sep='')
	names(sx)	=	sTerms
	p 			=	length(sTerms)
	#
	# Run MCMC to obtain estimate of true posterior distribution of theta.
	S		=	10000	# Number of MCMC samples
	for (s in 1 : S) { # s = 1 # for testing
		if (s == 1)	{# Starting values of theta for MCMC:
			Samples.theta		=	matrix(NA, nrow = S, ncol = d)
			Samples.logprob	=	matrix(NA, nrow = S, ncol = 1)
			theta.0			=	c(truth)
			mu.prior			=	matrix(0, ncol = d)
			diag1				=	diag(1, d)	
			b					=	0.1	# Controls prior correlation between theta[1] and theta[2]
		}
		#
		# Do Gibbs slice sampling update of theta[1]  (using shrinkage method of Neal 2003):
		## Sample the slice variable:
		logprob	=	dmvnorm(t(x), mean = theta.0, sigma = diag1, log = TRUE) + 	
							log(exp( -((theta.0[1]^2) / 200) - ((theta.0[2] - (b*(theta.0[1]^2)) + (100*b))^2)/2 ))
		z			=	logprob - rexp(1, rate = 1)
		#
		Accept			=	0
		Try				=	0
		theta.0.old	=	theta.0[1]
		L				=	-50
		R				=	50
		while (Accept == 0) {
			U			= 	runif(1, min = 0, max = 1)
			Try			= 	Try + 1
			theta.0[1]	= 	L + U*(R-L)
			logprob	= 	dmvnorm(t(x), mean = theta.0, sigma = diag1, log = TRUE) + 	
								log(exp( -((theta.0[1]^2) / 200) - ((theta.0[2] - (b*(theta.0[1]^2)) + (100*b))^2)/2 ))
			Accept		= 	ifelse(logprob > z, 1, 0)
			if (theta.0[1] <  theta.0.old) {L = theta.0[1]}
			if (theta.0[1] >= theta.0.old) {R = theta.0[1]}
		}
		#
		# Do Gibbs slice sampling update of theta[2] (using shrinkage method of Neal 2003):
		## Sample the slice variable:
		logprob	=	dmvnorm(t(x), mean = theta.0, sigma = diag1, log = TRUE) + 	
							log(exp( -((theta.0[1]^2) / 200) - ((theta.0[2] - (b*(theta.0[1]^2)) + (100*b))^2)/2 ))
		z			=	logprob - rexp(1, rate = 1)
		#
		Accept			=	0
		Try				=	0
		theta.0.old	=	theta.0[2]
		L				=	-50
		R				=	50
		while (Accept == 0) {
			U			= 	runif(1, min = 0, max = 1)
			Try			= 	Try + 1
			theta.0[2]	= 	L + U*(R-L)
			logprob	= 	dmvnorm(t(x), mean = theta.0, sigma = diag1, log = TRUE) + 	
								log(exp( -((theta.0[1]^2) / 200) - ((theta.0[2] - (b*(theta.0[1]^2)) + (100*b))^2)/2 ))
			Accept		= 	ifelse(logprob > z, 1, 0)
			if (theta.0[2] <  theta.0.old) {L = theta.0[2]}
			if (theta.0[2] >= theta.0.old) {R = theta.0[2]}
		}
		#
		# Do Gibbs sampling update of (theta3, ..., theta300):
		theta.0[3:d]	=	rnorm(x[3:d], x[3:d], sqrt(1/2))
		# Above full conditional posterior distribution is consistent with 
		# Bernardo & Smith (1994, p.439 on "normal model, known precision", 
		# and p.441, multivariate normal model).
		#
		Samples.theta[s,]	 	=	theta.0
		Samples.logprob[s]	=	dmvnorm(t(x), mean = theta.0, sigma = diag1, log = TRUE) + 	
									log(exp( -((theta.0[1]^2) / 200) - ((theta.0[2] - (b*(theta.0[1]^2)) + (100*b))^2)/2 - sum(theta.0[3:d]^2) ))	
		if (s/10 == round(s/10,0)){cat("MCMC iteration", s, "of", S, " \n"); flush.console()}
	} 
	PostMean.theta	=	colMeans(Samples.theta) 			# marginal Posterior means
	PostMed.theta		=	apply(Samples.theta, 2, median)	# marginal Posterior medians
	PostMode.theta	=	Samples.theta[Samples.logprob == max(Samples.logprob),]# Biau etal. Multivariate mode estimator
	MLE.theta			= 	PostMode.theta # MLE is posterior mode under uniform prior
	SD.theta			=	apply(Samples.theta, 2, sd)	# marginal Posterior standard deviations
	#
	postEtheta.exact		=	PostMean.theta
	postMEDtheta.exact	=	PostMed.theta
	posteriorMode.exact	=	PostMode.theta
	MLE.exact				=	MLE.theta
	postSDtheta.exact		=	SD.theta
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for Price model
	# ====================================================================================================================================
	N				=	10000 #
	n.table 		= 	1;	# Sample size of the reference table
	theta.table	=	matrix(NA, nrow = N, ncol = d)
	sy.table		=	matrix(NA, nrow = N, ncol = length(sTerms))
	colnames(theta.table)	<-	parameterNames
	colnames(sy.table)		<-	sTerms
	# 
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j	=	c(rnorm(1, mean = 0, sd = sqrt(100)), rnorm(d - 1, mean = 0, sd = 1) )
		theta.j[2]	=	theta.j[2] + b*(theta.j[1]^2) - 100*b # a sample from the twisted normal prior (Li et al. 2017, p.81)
		#
		# Simulate from model's prior predictive:
		y			=	rmvnorm(n.table, mean = theta.j, sigma = diag1)
		# Calculate summary statistics from the generated pseudo-data:
		sy.j		=	y
		#
		# Update Reference Table
		theta.table[j,]	=	theta.j
		sy.table[j,]		=	sy.j
		if ((j / 10) == round(j / 10))	{	flush.console()	}
	}
	# Restrict to subset of theta samples corresponding to acceptable 
	isFinite.sy	=	apply(is.finite(sy.table),1,all)
	notNA.sy		=	apply(! is.na(sy.table),1,all)
	isOK.sy		=	isFinite.sy & notNA.sy
	theta.table 	= 	theta.table[isOK.sy, ]
	sy.table		=	sy.table[isOK.sy, ]
	N 				=	sum(isOK.sy)
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
		if (k <= 2)	{IndRF 	= 	1:2}	# Using selected summary statistics.
		if (k  > 2)	{IndRF 	= 	k}	 	# Using selected summary statistics.
		sy.table0			=	as.matrix(sy.table[,IndRF])
		drf.forest 		=	drf(X = sy.table0, Y = theta.table[,k])#, compute.variable.importance = TRUE # Costly to compute
		X.test				=	matrix(sx[IndRF], nrow = 1) # Using selected summary statistics.
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
		cat("Distribution Random Forest", k, "of", d, "trained and predicted.", " \n"); flush.console()
	}
	# Remove rows for which u_jk = 0 or 1. (they correspond to theta_j's with zero posterior density anyway).
	# keep					=	apply((u > 0) & (u < 1), 1, all)
	# u.keep				=	u[keep,]
	# theta.table.keep	=	theta.table[keep, ]
	# pi.thetaRT.keep 	= 	pi.thetaRT[keep,]
	#
	# I found that pi.thetaRT.keep produced many zeros and u produced many zeros and ones
	# leading to all zeros for the posteriorPDFs computations, below.
	# Therefore, I tried using kernel smoothing for marginal densities and inverse cdfs, below,
	# which mitigated this issue.
	for (k in 1:p)	{
		if (k==1){		pi.thetaRT.kernel		=	matrix(NA, nrow = N, ncol = p)
						u.kernel				=	matrix(NA, nrow = N, ncol = p)			}
		fitKernelRT 				= 	kde1d(theta.table[,k], weights = theta.DRFweights.table[,k])
		pi.thetaRT.kernel[,k]	=	dkde1d(theta.table[,k], fitKernelRT)
		u.kernel[,k]				=	pkde1d(theta.table[,k], fitKernelRT)
	}
	#
	#"multiplied by N / (N + 1) to avoid evaluating the [copula] density at the edges of the unit square."
	# (from p. 499 of Genest & Neslehova, 2007, Astin Bulletin)
	# See also Genest et al. (1995) and Okhrin 2012, Ch.17, p.484, "Fitting High-Dimensional Copulae to Data"
	u.kernel			=	u.kernel * ( N / ( N + 1 ) )
	keep				=	apply((u.kernel > 0) & (u.kernel < 1), 1, all) & (apply(pi.thetaRT.kernel,1,prod)>0) # sum(keep)
	u.keep				=	u.kernel[keep,]
	theta.table.keep	=	theta.table[keep, ]
	pi.thetaRT.keep 	= 	pi.thetaRT.kernel[keep,]
	#
	# Estimate copula parameters (df, correlation matrix) of Meta-t posterior distribution:
	OutCopulaFit  	=	fitStudentcopula(u.keep,fit.method = "EM-MLE", df.bounds = c(0.1, 1000), verbose = TRUE)
	post.df			=	OutCopulaFit$df
	post.Correlations	=	OutCopulaFit$scale
	# Meta-t posterior densities of (retained) unordered theta.j's from Reference table
	posteriorPDFs		=	dStudentcopula(u.keep,df=post.df,scale=post.Correlations)*apply(pi.thetaRT.keep,1,prod)
	#
	# Find posterior mode and MLE of theta using the sampled values:
	posteriorMode		=	theta.table.keep[posteriorPDFs == max(posteriorPDFs, na.rm = T),]
	priorPDFs 			= 	exp( -((theta.table.keep[,1]^2) / 200) - ((theta.table.keep[,2] - (b*(theta.table.keep[,1]^2)) + (100*b))^2)/2 - rowSums(theta.table.keep[,3:d]^2) )	
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
	end_time 									= 	Sys.time()
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
	#
	# Do weighted one-sample (two-tailed) Kolmogorov Smirnov tests (Monahan,2011,p.358)
	KS.Dn.simulations0				=	c()
	KS.stat.simulations0				=	c()
	for (k in 1:p) { # k = 1;	# For testing
		theta.with.weights			=	cbind(theta.table[,k], theta.DRFweights.table[,k])
		sorted.theta.with.weights	=	theta.with.weights[order(theta.with.weights[,1]), ]
		# Univariate local-polynomial (degree =2) likelihood kernel density estimation (from MCMC samples)
		# using plug-in bandwidth of Sheather and Jones (1991) :
		fit 		= 	kde1d(Samples.theta[,k])
		F 			= 	pkde1d(sorted.theta.with.weights[,1], fit)#cdf estimates of true marginal posterior distribution (from MCMC samples).
		Fk			=	cumsum(sorted.theta.with.weights[,2])
		Fk.minus.1	=	c(0, Fk[1 : (length(Fk) - 1) ] )
		Dn 			= 	max(apply(cbind(abs(Fk - F), abs(F - Fk.minus.1)), 1, max))# Kolmogorov-Smirnov statistic
		S.2			=	sum(sorted.theta.with.weights[,2] ^ 2)
		KSstat		=	(1 / S.2) * Dn
		# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)
		KS.Dn.simulations0	=	c(KS.Dn.simulations0, Dn		)
		KS.stat.simulations0	=	c(KS.stat.simulations0, KSstat	)
	}
	KS.Dn.simulations			=	rbind(KS.Dn.simulations, 	KS.Dn.simulations0	)
	KS.stat.simulations		=	rbind(KS.stat.simulations, 	KS.stat.simulations0	)
	#
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
	outputFileName	=	paste("Twisted prior model p = 300 ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	outputFileNameSummaryTable	=	paste("TwistedNormal SimulationStudySummaryTable ",gsub("\\:", "_", round(end_time)),".csv",sep="")
	write.csv(SimulationStudySummaryTable, file = outputFileNameSummaryTable)
	outputFileNameSummaryTableCorrel	=	paste("TwistedNormal SimulationStudySummaryTableCorrel ",gsub("\\:", "_", round(end_time)),".csv",sep="")
	write.csv(SimulationStudySummaryTableCorrel, file = outputFileNameSummaryTableCorrel)
	print(paste('Replica', r, 'done'))
	flush.console()
	# ====================================================================================================================================
}