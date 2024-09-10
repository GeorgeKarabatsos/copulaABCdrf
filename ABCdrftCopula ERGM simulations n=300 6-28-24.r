install.packages('igraph');			library('igraph')		# For network data analysis/management
install.packages('ergm'); 			library('ergm')		# For network data analysis/management
install.packages('ergm.count'); 	library('ergm.count')# For network data analysis/management
install.packages('ergm.rank'); 	library('ergm.rank')	# For network data analysis/management
install.packages('drf');			library('drf')		# Distributional Random Forests.
install.packages('nvmix');			library('nvmix')		# For t copula.
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('ddpcr'); 			library('ddpcr')		# For quiet() 
#setwd('C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork')
setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')

rm(list = ls())
source('network statistics 9-14-23.R')
source('simulateNetwork 9-14-23.R')
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Run simulation over replicate datasets:
# ====================================================================================================================================
replicas		=	10
for (r in 1 : replicas) {  # r = 1
	# ====================================================================================================================================
	# Simulate network dataset from ERGM (network dataset to be analyzed using ABC)
	# ====================================================================================================================================
	# Network data 
	Nnodes 			= 	300
								sigma = -0.2	; 	tau = 0.5; 	thetadegree	= .8;
	truth 				= 	c(	sigma			, 	tau		,	thetadegree			)
	parameterNames	=	c('kstar(2)'	, 	'triangles',	'degree1.5'	)
	d					=	length(parameterNames)
	init.G 			= 	network(Nnodes, directed = FALSE, density = 0.1)
	#init.G 			= 	matrix(rbinom(Nnodes ^ 2, 1, 0.005), Nnodes, Nnodes); diag(init.G) = 0;
	#init.G 			= 	as.network.matrix(init.G, directed = FALSE)
	G 					= 	simulate(init.G ~ kstar(2) + triangles + degree1.5, nsim = 1, coef = truth)# in 'network' format
	#
	# Monte Carlo Maximum likelihood estimates and standard error:
	MCMLEout			=	ergm(formula = G ~ kstar(2) + triangles + degree1.5, estimate = "MLE", verbose = FALSE)
	MCMLE				=	coef(MCMLEout)
	seMCMLE			=	sqrt(diag(vcov(MCMLEout)))
	#seMPLEgodambe		=	sqrt(diag(vcov(	ergm(	formula = G ~ kstar(2) + triangles, estimate = "MPLE", 
	#												control = control.ergm(MPLE.covariance.method = "Godambe")))))
	# ==========================================================================================
	# Set up the g prior for ERG model:  
	# ==========================================================================================
	# Now set up the g-prior for the ERG model, given by 
	# theta 	~ 	Normal(0, g*invJ)
	# The g prior, such as the one above, was originally proposed for generalized linear models
	# (Wang & George, 2007; see also Sabanés Bové & Held, 2011, p.392).
	# In the prior above, invJ is the inverse matrix (the inverse of matrix J) representing the 
	# (estimated) sampling covariance matrix of the Monte Carlo Maximum Likelihood estimate
	# (i.e., the MCMLE) of the model parameter theta from the network dataset
	# (Hummel, Hunter, & Handcock, 2012, and references therein).
	# J is the observed Fisher information matrix (i.e., negative Hessian matrix) 
	# evaluated at the observed data point (network) and at the MCMLE of theta.
	# For the fixed choice g = 1 (i.e., sample size is one network observation)
	# the g-prior above coincides with the unit-information prior (Kass & Wasserman, 1995).
	# References
	# Hummel, R., Hunter, D., & Handcock, M. (2012). Improving simulation-based 
	#   algorithms for fitting ERGMs. Journal of Computational and Graphical 
	#   Statistics, 21, 920-939.
	# Kass, R., & Wasserman, L. (1995). A reference Bayesian test for nested 
	#   hypotheses and its relationship to the Schwarz criterion. Journal of the
	#   American Statistical Association, 90, 928-934.
	# Sabanés Bové, D., & Held, L. Hyper-g priors for generalized linear models. 
	#   Bayesian Analysis, 6, 387-410. 
	# Wang, X., & George, E. (2007). Adaptive Bayesian criteria in variable selection
	#   for generalized linear models. Statistica Sinica, 17, 667-690.
	# ------------------------------------------------------------------------------------------
	# g-prior setup:
	# ------------------------------------------------------------------------------------------
	invJ 				=	vcov(	ergm(	formula = G ~ kstar(2) + triangles + degree1.5, estimate = "MPLE"))
	# Command line above extracts invJ(thetahat) (i.e.,the inverse of matrix J(thetahat)), 
	# where J(thetahat) is the negative Hessian matrix of the log-likelihood evaluated at the MLE thetahat. 
	muPrior 			= 	matrix(rep(0, d), nrow = 1) # mean of g prior
	g 					=	100000;
	ginvJ 				= 	g * invJ 		# covariance matrix of g-prior.
	SigmaPrior			=	ginvJ			# covariance matrix of g-prior.
	# SigmaPrior		=	diag(10,d, d) # Alternative prior covariance matrix.
	# ====================================================================================================================================
	# Calculate MPLE summary statistics from network G (using ERGM sufficient statistics) 
	# ====================================================================================================================================
	outMPLEx			=	MPLEergm(G, c('kstar(2)', 'triangles', 'degree1.5'))
	MPLEx 				= 	outMPLEx$MPLE[1,]
	seMPLE				=	sqrt(outMPLEx$MPLE[2,])
	MPLEterms			=	outMPLEx$ergmTerms
	p 					=	length(MPLEterms)
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for all network summaries from ergm package.
	# ------------------------------------------------------------------------------------------------------------------------------------
	# MPLE_ergm = c();		Terms_ergm = c()
	# out	=	MPLEergm(G);  MPLE_ergm[[1]] = 	out$MPLE;	Terms_ergm[[1]]	=	out$Terms
	# out2	=	MPLEergm(G, Terms_ergm[[1]]); # Using all (selected) terms.
	# rbind(out$MPLE[1,], out2$MPLE[1,]); # 	MPLEs 					# Results match.
	# rbind(out$MPLE[2,], out2$MPLE[2,])	#	variances(MPLEs)		# Results match.
	# c(ncol(out$MPLE), ncol(out2$MPLE), all(out$MPLE == out2$MPLE), all(out2$ergmTerms == out$ergmTerms))
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for selected network summaries from ergm package:
	# ------------------------------------------------------------------------------------------------------------------------------------
	# MPLE_ergm = c();		Terms_ergm = c()
	# MPLEergm(G, c('edges', 'triangles')); # Equals:#   MPLEergm(G, c('triangles', 'edges'))
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for all network summaries from igraph package:
	# ------------------------------------------------------------------------------------------------------------------------------------
	# MPLEigraph(G) # Computation takes a while when Nnodes = 20
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for selected network summaries from igraph package:
	# ------------------------------------------------------------------------------------------------------------------------------------
	# MPLEigraph(G, c('assortativity_degree', 'diameter', 'var_node_triangles' ))
	# MPLEigraph(G, c('var_node_triangles', 'diameter')); # MPLEigraph(G, 'diameter')
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for ERGM
	# ====================================================================================================================================
	N				=	10000			# Sample size of the reference table
	NnodesTable 	= 	floor(Nnodes * 1);
	theta.table	=	matrix(NA, nrow = N, ncol = d)
	MPLE.table		=	matrix(NA, nrow = N, ncol = length(MPLEterms))
	colnames(theta.table)	<-	parameterNames
	colnames(MPLE.table)		<-	MPLEterms
	# 
	# Generate Reference table:
	for (j in 1 : N) { # j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j		=	rmvnorm(1, mean = muPrior, sigma = SigmaPrior)
		# Simulate network dataset: 
		init.G 		= 	network(NnodesTable, directed = FALSE, density = 0.1)
		#init.G 		= 	matrix(rbinom(NnodesTable ^ 2, 1, 0.005), NnodesTable, NnodesTable); diag(init.G) = 0;
		#init.G 		= 	as.network.matrix(init.G, directed = FALSE)
		Gtable			= 	simulate(init.G ~ kstar(2) + triangles + degree1.5, nsim = 1, coef = theta.j)# in 'network' format
		# Compute MPLE summary statistics of the simulated data:
		outMPLEy			=	MPLEergm(Gtable, MPLEterms)
		MPLEy.j			= 	outMPLEy$MPLE[1,]
		# Update Reference Table
		theta.table[j,]	=	theta.j
		MPLE.table[j,]	=	MPLEy.j
		if ((j / 10) == round(j / 10))	{	flush.console()	}
	}
	# Restrict to subset of theta samples corresponding to nondegenerate netorks (finte MPLEs) 
	isFiniteMPLE		=		apply(is.finite(MPLE.table),1,all)
	theta.table 		= 		theta.table[isFiniteMPLE,]
	MPLE.table			=		MPLE.table[isFiniteMPLE,]
	N 						=		sum(isFiniteMPLE)
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Construct multivariate meta-t copula ABC posterior distribution based on DRF marginals.
	# ====================================================================================================================================
	# For each parameter (dependent variable), and set of predictors (MPLEs), 
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
		drf.forest 		=	drf(X = MPLE.table, Y = theta.table[,k])#, compute.variable.importance = TRUE # Costly to compute
		X.test				=	MPLEx
		# Use trained DRF to estimate (predict) posterior mean, standard deviation, and quantiles, conditional on MPLEx:
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
	priorPDFs			=	dmvnorm(theta.table.keep, mean = muPrior, sigma = SigmaPrior)
	likelihoods		=	posteriorPDFs / priorPDFs
	MLE					=	theta.table.keep[likelihoods == max(likelihoods, na.rm = T),]
	# ------------------------------------------------------------------------------------------------------------------------------------
	# It is tempting to try to optimize (using random meta-t values of theta, or using optim())
	# to get posteriorMode and MLE, using their values obtained above as starting values.
	# However doing so would still use the marginal (discrete) posterior PDF (histogram)
	# and CDF values of the original theta.jk's samples from the Reference Table.
	# That is, the marginal probabilities of any new trial values of theta (random or searched) 
	# would not be considered and not calculated (or difficult to calculate). 
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
		postMEDtheta.simulations			=	c()
		posteriorMode.simulations			=	c()
		MLE.simulations						=	c()
		MCMLE.simulations					=	c()
		MPLE.simulations						=	c()
		L1error.postEtheta.simulations		=	c()
		L2error.postEtheta.simulations		=	c()
		L1error.postMEDtheta.simulations	=	c()
		L2error.postMEDtheta.simulations	=	c()
		L1error.posteriorMode.simulations	=	c()
		L2error.posteriorMode.simulations	=	c()
		L1error.MLE.simulations				=	c()
		L2error.MLE.simulations				=	c()
		L1error.MCMLE.simulations			=	c()
		L2error.MCMLE.simulations			=	c()
		L1error.MPLE.simulations			=	c()
		L2error.MPLE.simulations			=	c()
		postSDtheta.simulations				=	c()
		seMCMLE.simulations					=	c()
		seMPLE.simulations					=	c()
		post95cover.simulations				=	c()
		postIQRcover.simulations			=	c()
		MCMLE95cover.simulations			=	c()
		MPLE95cover.simulations				=	c()
		post.df.simulations					=	c()
		post.Correlations.simulations		=	c()
		TotalComputationTimePerSimulation	=	c()
	}

	# Update collections of simulation results:
	postEtheta.simulations					=	rbind(postEtheta.simulations, 				postEtheta							)
	postMEDtheta.simulations				=	rbind(postMEDtheta.simulations, 			postQtheta[3,]					)
	posteriorMode.simulations				=	rbind(posteriorMode.simulations, 			posteriorMode						)
	MLE.simulations							=	rbind(MLE.simulations, 						MLE									)
	MCMLE.simulations							=	rbind(MCMLE.simulations, 					MCMLE								)
	MPLE.simulations							=	rbind(MPLE.simulations, 					MPLEx								)
	L1error.postEtheta.simulations			=	rbind(L1error.postEtheta.simulations, 	abs(postEtheta - truth) 		)
	L2error.postEtheta.simulations			=	rbind(L2error.postEtheta.simulations,    (postEtheta - truth)^2 			)
	L1error.postMEDtheta.simulations		=	rbind(L1error.postMEDtheta.simulations,	abs(postQtheta[3,] - truth)		)
	L2error.postMEDtheta.simulations		=	rbind(L2error.postMEDtheta.simulations,  (postQtheta[3,] - truth)^2 		)
	L1error.posteriorMode.simulations		=	rbind(L1error.posteriorMode.simulations, abs(posteriorMode - truth) 		)
	L2error.posteriorMode.simulations		=	rbind(L2error.posteriorMode.simulations,    (posteriorMode - truth)^2 	)
	L1error.MLE.simulations					=	rbind(L1error.MLE.simulations, 			abs(MLE - truth) 					)
	L2error.MLE.simulations					=	rbind(L2error.MLE.simulations,    				(MLE - truth)^2 				)
	L1error.MCMLE.simulations				=	rbind(L1error.MCMLE.simulations, 			abs(MCMLE - truth) 				)
	L2error.MCMLE.simulations				=	rbind(L2error.MCMLE.simulations,    			(MCMLE - truth)^2 			)
	L1error.MPLE.simulations				=	rbind(L1error.MPLE.simulations, 			abs(MPLEx - truth) 				)
	L2error.MPLE.simulations				=	rbind(L2error.MPLE.simulations,    			(MPLEx - truth)^2 				)
	postSDtheta.simulations					=	rbind(postSDtheta.simulations, 			postSDtheta						)
	seMCMLE.simulations						=	rbind(seMCMLE.simulations, seMCMLE)
	seMPLE.simulations						=	rbind(seMPLE.simulations, seMPLE)
	post95cover.simulations					=	rbind(post95cover.simulations,	ifelse((truth > postQtheta[1,]) & (truth < postQtheta[5,]),1,0))
	postIQRcover.simulations				=	rbind(postIQRcover.simulations,ifelse((truth > postQtheta[2,]) & (truth < postQtheta[4,]),1,0))
	MCMLE95cover.simulations				=	rbind(MCMLE95cover.simulations,(truth > MCMLE - (1.96 * seMCMLE)) & (truth < MCMLE + (1.96 * seMCMLE)))
	MPLE95cover.simulations					=	rbind(MPLE95cover.simulations,	(truth > MPLEx -  (1.96 * seMPLE )) & (truth < MPLEx   + (1.96 * seMPLE)))
	post.df.simulations						=	rbind(post.df.simulations, rep(post.df, d))
	post.Correlations.simulations			=	rbind(post.Correlations.simulations, c(post.Correlations))
	TotalComputationTimePerSimulation		=	rbind(TotalComputationTimePerSimulation, rep(TotalComputationTimeSimulation, d))


	if (r == replicas){
		Mean.postEtheta.simulations			=	apply(postEtheta.simulations, 2, mean, na.rm = T);			SD.postEtheta.simulations		=	apply(postEtheta.simulations, 2, sd, na.rm = T)
		Mean.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, mean, na.rm = T);			SD.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, sd, na.rm = T)
		Mean.posteriorMode.simulations		=	apply(posteriorMode.simulations, 2, mean, na.rm = T);		SD.posteriorMode.simulations	=	apply(posteriorMode.simulations, 2, sd, na.rm = T)
		Mean.MLE.simulations					=	apply(MLE.simulations, 2, mean, na.rm = T);					SD.MLE.simulations				=	apply(MLE.simulations, 2, sd, na.rm = T)
		Mean.MCMLE.simulations				=	apply(MCMLE.simulations, 2, mean, na.rm = T);					SD.MCMLE.simulations				=	apply(MCMLE.simulations, 2, sd, na.rm = T)
		Mean.MPLE.simulations				=	apply(MPLE.simulations, 2, mean, na.rm = T);					SD.MPLE.simulations				=	apply(MPLE.simulations, 2, sd, na.rm = T)
		MAE.postEtheta.simulations			=	apply(L1error.postEtheta.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.postEtheta.simulations			=	apply(L2error.postEtheta.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.postMEDtheta.simulations		=	apply(L1error.postMEDtheta.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.postMEDtheta.simulations		=	apply(L2error.postMEDtheta.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.posteriorMode.simulations		=	apply(L1error.posteriorMode.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.posteriorMode.simulations		=	apply(L2error.posteriorMode.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.MLE.simulations					=	apply(L1error.MLE.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.MLE.simulations					=	apply(L2error.MLE.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.MCMLE.simulations				=	apply(L1error.MCMLE.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.MCMLE.simulations				=	apply(L2error.MCMLE.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.MPLE.simulations					=	apply(L1error.MPLE.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.MPLE.simulations					=	apply(L2error.MPLE.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		Mean.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, mean, na.rm = T);			SD.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, sd, na.rm = T)
		Mean.seMCMLE.simulations			=	apply(seMCMLE.simulations, 2, mean, na.rm = T);				SD.seMCMLE.simulations		=	apply(seMCMLE.simulations, 2, sd, na.rm = T)
		Mean.seMPLE.simulations				=	apply(seMPLE.simulations, 2, mean, na.rm = T);				SD.seMPLE.simulations		=	apply(seMPLE.simulations, 2, sd, na.rm = T)
		Prob.post95cover.simulations		=	apply(post95cover.simulations, 2, mean, na.rm = T)
		Prob.postIQRcover.simulations		=	apply(postIQRcover.simulations, 2, mean, na.rm = T)
		Prob.MCMLE95cover.simulations		=	apply(MCMLE95cover.simulations, 2, mean, na.rm = T)
		Prob.MPLE95cover.simulations		=	apply(MPLE95cover.simulations, 2, mean, na.rm = T)
		Mean.post.df.simulations			=	apply(post.df.simulations, 2, mean, na.rm = T);				SD.post.df.simulations					=	apply(post.df.simulations, 2, sd, na.rm = T)
		Mean.post.Correlations.simulations	=	apply(post.Correlations.simulations, 2, mean, na.rm = T);		SD.post.Correlations.simulations	=	apply(post.Correlations.simulations, 2, sd, na.rm = T)
		Mean.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, mean, na.rm = T);	SD.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, sd, na.rm = T)

		SimulationStudySummaryTable			=	rbind(	
		truth,
		Mean.postEtheta.simulations,		SD.postEtheta.simulations,
		Mean.postMEDtheta.simulations,		SD.postMEDtheta.simulations,
		Mean.posteriorMode.simulations,	SD.posteriorMode.simulations,
		Mean.MLE.simulations,				SD.MLE.simulations,
		Mean.MCMLE.simulations,				SD.MCMLE.simulations,
		Mean.MPLE.simulations,				SD.MPLE.simulations,
		MAE.postEtheta.simulations,
		MSE.postEtheta.simulations,
		MAE.postMEDtheta.simulations,
		MSE.postMEDtheta.simulations,
		MAE.posteriorMode.simulations,
		MSE.posteriorMode.simulations,
		MAE.MLE.simulations,
		MSE.MLE.simulations,
		MAE.MCMLE.simulations,
		MSE.MCMLE.simulations,
		MAE.MPLE.simulations,
		MSE.MPLE.simulations,
		Mean.postSDtheta.simulations,		SD.postSDtheta.simulations,
		Mean.seMCMLE.simulations,			SD.seMCMLE.simulations	,
		Mean.seMPLE.simulations,			SD.seMPLE.simulations	,
		Prob.post95cover.simulations,	
		Prob.postIQRcover.simulations,
		Prob.MCMLE95cover.simulations,
		Prob.MPLE95cover.simulations,
		Mean.post.df.simulations,			SD.post.df.simulations,
		Mean.ComputationTimePerSimulation,	SD.ComputationTimePerSimulation,
		sum(TotalComputationTimeAllSimulations)								)
	
		rownames(SimulationStudySummaryTable)=	c(		
		'truth',
		'Mean.postEtheta.simulations',		'SD.postEtheta.simulations',
		'Mean.postMEDtheta.simulations',	'SD.postMEDtheta.simulations',
		'Mean.posteriorMode.simulations',	'SD.posteriorMode.simulations',
		'Mean.MLE.simulations',				'SD.MLE.simulations',
		'Mean.MCMLE.simulations',			'SD.MCMLE.simulations',
		'Mean.MPLE.simulations',			'SD.MPLE.simulations',
		'MAE.postEtheta.simulations',
		'MSE.postEtheta.simulations',
		'MAE.postMEDtheta.simulations',
		'MSE.postMEDtheta.simulations',
		'MAE.posteriorMode.simulations',
		'MSE.posteriorMode.simulations',
		'MAE.MLE.simulations',
		'MSE.MLE.simulations',
		'MAE.MCMLE.simulations',
		'MSE.MCMLE.simulations',
		'MAE.MPLE.simulations',
		'MSE.MPLE.simulations',
		'Mean.postSDtheta.simulations',	'SD.postSDtheta.simulations',
		'Mean.seMCMLE.simulations',			'SD.seMCMLE.simulations'	,
		'Mean.seMPLE.simulations',			'SD.seMPLE.simulations'	,
		'Prob.post95cover.simulations',	
		'Prob.postIQRcover.simulations',
		'Prob.MCMLE95cover.simulations',
		'Prob.MPLE95cover.simulations',
		'Mean.post.df.simulations',			'SD.post.df.simulations',
		'Mean.ComputationTimePerSimulation','SD.ComputationTimePerSimulation',
		'TotalComputationTimeAllSimulations'									)													
		colnames(SimulationStudySummaryTable)=	c(parameterNames)
		SimulationStudySummaryTableCorrel	=	rbind(	Mean.post.Correlations.simulations,	SD.post.Correlations.simulations)
		rownames(SimulationStudySummaryTableCorrel)=	c('Mean.post.CorrelationMatrix.simulations', 'SD.post.CorrelationsMatrix.simulations')
	}
	outputFileName	=	paste("ERGM simul n = 300 NnodesTable = ", NnodesTable, " ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	print(paste('Replica', r, 'done'))
	flush.console()
	# ====================================================================================================================================
}