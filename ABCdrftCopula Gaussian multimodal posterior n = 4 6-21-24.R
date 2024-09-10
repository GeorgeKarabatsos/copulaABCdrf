install.packages('drf');		library('drf')		# Distributional Random Forests.
install.packages('nvmix');		library('nvmix')		# For t copula.
install.packages('mvtnorm');	library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('MASS');		library('MASS')		# For mvrnorm().
install.packages('ddpcr'); 		library('ddpcr')		# For quiet() 
install.packages('kde1d');		library('kde1d')# for kernel CDF() used in Kolmogorov Smirnov test
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork')

rm(list = ls())
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Run simulation over replicate datasets:
# ====================================================================================================================================
replicas		=	10
for (r in 1 : replicas) {	# r = 1
	# ====================================================================================================================================
	# Simulate network dataset from bivariate Gaussian model having a 4-mode posterior distribution.
	# ====================================================================================================================================
	parameterNames=	c('theta1', 'theta2', 'theta3', 'theta4', 'theta5')
	d		=	length(parameterNames)
	truth	= 	c(-0.7, -2.9, -1.0, -0.9, 0.6); # True data-generating model parameters to be estimated.
	theta	=	truth
	# theta is assigned a Uniform prior on [-3,3] x [-4,4] x [-3,3] x [-3,3] x [-3,3]
	n		=	4	# Dataset sample size ( to be simulated)
	#
	# Model data simulation:
	mu		=	theta[1:2];  s1	=	theta[3]^2;	s2	=	theta[4]^2;	rho	=	tanh(theta[5])
	Sigma	=	matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), nrow=2, ncol=2)
	# Simulate data from the model:
	x		=	rmvnorm(n, mean = mu, sigma = Sigma)
	# Data summaries:
	xbar		=	matrix(colMeans(x), nrow = 1)
	CovMatrix	=	cov.wt(x, method = "ML")$cov[c(1,4,3)]# variances are first two entries of cov.wt row vector output.
	#
	# Run MCMC to obtain estimate of true posterior distribution of theta.
	S		=	10000	# Number of MCMC samples
	for (s in 1 : S) {
		if (s == 1)	{
			# Starting values of theta for MCMC:
			mu.0	=	mu;		s1.0 =	s1;	 s2.0 = s2;	rho.0 = rho
			Sigma.0		=	matrix(c(s1.0^2, rho.0 * s1.0* s2.0, rho.0 * s1.0 * s2.0, s2.0^2), nrow=2, ncol=2)
			theta.0		=	c(mu, sqrt(s1.0), sqrt(s2.0), atanh(rho.0))
			Samples.theta		=	matrix(NA, nrow = S, ncol = 5)
			Samples.logprob	=	matrix(NA, nrow = S, ncol = 1)
		}
		# Do Gibbs sampling update of mu.0 with rejection sampling.
		Accept		=	0
		while (Accept == 0) {
			mu.0	=	rmvnorm(1, mean = xbar, sigma = Sigma.0 / n)
			Accept	=	ifelse((mu.0[1] >= -3) & (mu.0[1] <= 3) & (mu.0[2] >= -4) & (mu.0[2] <= 4), 1, 0)
		}
		theta.0[1:2]	=	mu.0
		#
		# Do Gibbs sampling update of (s1, s2, rho) with rejection slice sampling:
		## Sample the slice variable:
		logprob	=	sum(dmvnorm(x, mean = mu.0, sigma = Sigma.0, log = TRUE)) + sum(dunif(theta.0[3:5], min = -3, max = 3, log = TRUE))
		z			=	logprob - rexp(1, rate = 1)
		#
		Accept			=	0
		try				=	0
		theta.0.old	=	theta.0
		while (Accept == 0) {
			try			=	try + 1
			theta.0[3]	=	runif(1, min = -3, max = 3);	s1.0	=	theta.0[3]^2;
			theta.0[4]	=	runif(1, min = -3, max = 3);	s2.0	=	theta.0[4]^2;
			theta.0[5]	=	runif(1, min = -3, max = 3);	rho.0	=	tanh(theta.0[5])
			Sigma.0	=	matrix(c(s1.0^2, rho.0 * s1.0 * s2.0, rho.0 * s1.0 * s2.0, s2.0^2), nrow=2, ncol=2)
			logprob	=	sum(dmvnorm(x, mean = mu.0, sigma = Sigma.0, log = TRUE)) + sum(dunif(theta.0[3:5], min = -3, max = 3, log = TRUE))
			Accept		=	ifelse(logprob > z, 1, 0)
			if ((Accept == 0) & (try == 300)) {theta.0 = theta.0.old}
			if(try == 300) { break }
		}
		Samples.theta[s,]	 	=	theta.0
		Samples.logprob[s]	=	logprob
		if (s/10 == round(s/10,0)){cat("MCMC iteration", s, "of", S, " \n"); flush.console()}
	} 
	PostMean.theta	=	colMeans(Samples.theta) 			# marginal Posterior means
	PostMed.theta		=	apply(Samples.theta, 2, median)	# marginal Posterior medians
	PostMode.theta	=	Samples.theta[Samples.logprob == max(Samples.logprob),]# Biau etal. Multivariate mode estimator
	MLE.theta			= 	PostMode.theta # MLE is posterior mode under uniform prior
	SD.theta			=	apply(Samples.theta, 2, sd)	# marginal Posterior standard deviations
	#
	# ====================================================================================================================================
	# Summary statistics from data x 
	# ====================================================================================================================================
	sTerms		=	c('xbar1', 'xbar2', 'var1', 'var2', 'Cov12')
	sx 			=	c( xbar, 	CovMatrix  	)	# variances are first two entries of CovMatrix row vector.
	names(sx)	=	sTerms
	p 			=	length(sTerms)
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
	n.table 		= 	n;	# Sample size of the reference table
	theta.table	=	matrix(NA, nrow = N, ncol = d)
	sy.table		=	matrix(NA, nrow = N, ncol = length(sTerms))
	colnames(theta.table)	<-	parameterNames
	colnames(sy.table)		<-	sTerms
	# 
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j	=	c(runif(1, min = -3, max = 3),runif(1, min = -4, max = 4), runif(3, min = -3, max = 3))
		#
		# Simulate from model's prior predictive:
		mu.j		=	theta.j[1:2];  s1.j	=	theta.j[3]^2;	  s2.j	=	theta.j[4]^2;	  rho.j  =  tanh(theta.j[5])
		Sigma.j	=	matrix(c(s1.j^2, rho.j * s1.j * s2.j, rho.j * s1.j * s2.j, s2.j^2), nrow = 2, ncol = 2)
		# Simulate data from the model:
		y			=	rmvnorm(n.table, mean = mu.j, sigma = Sigma.j)
		# Calculate summary statistics from the generated pseudo-data:
		ybar.j			=	matrix(colMeans(y), nrow = 1)
		CovMatrix.j	=	cov.wt(y, method = "ML")$cov[c(1,4,3)]# variances are first two entries of cov.wt row vector output.
		sy			=	c(ybar.j, CovMatrix.j)
		names(sy)	<-	sTerms
		sy.j		= 	sy
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
	priorPDFs			=	apply(	cbind(	dunif(theta.table.keep[,1],min=-3,max=3),dunif(theta.table.keep[,2],min=-4,max=4),
											dunif(theta.table.keep[,3],min=-3,max=3),dunif(theta.table.keep[,4],min=-3,max=3),
											dunif(theta.table.keep[,5],min=-3,max=3)), 1, prod)
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
		# theta is assigned a Uniform prior on [-3,3] x [-4,4] x [-3,3] x [-3,3] x [-3,3]
		# Univariate local-polynomial (degree =2) likelihood kernel density estimation (from MCMC samples)
		# using plug-in bandwidth of Sheather and Jones (1991) :
		fit 		= 	kde1d(Samples.theta[,k], xmin = ifelse(k != 2, -3, -4), xmax = ifelse(k != 2,  3,  4))
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
	outputFileName	=	paste("Multimodal Gaussian simul n = 4 ntable = 4 ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	print(paste('Replica', r, 'done'))
	flush.console()
	# ====================================================================================================================================
}