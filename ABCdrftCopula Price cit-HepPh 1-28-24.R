install.packages('igraph');			library('igraph')		# For network data analysis/management
install.packages('ergm'); 			library('ergm')		# For network data analysis/management
install.packages('ergm.count'); 	library('ergm.count')# For network data analysis/management
install.packages('ergm.rank'); 	library('ergm.rank')	# For network data analysis/management
install.packages('drf');			library('drf')		# Distributional Random Forests.
install.packages('nvmix');			library('nvmix')		# For t copula.
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('ddpcr'); 			library('ddpcr')		# For quiet() 
setwd('C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork')
#setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')

rm(list = ls())
source('network statistics 9-14-23.R')
source('simulateNetwork 9-14-23.R')
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Analyze real dataset using Price model:
# ====================================================================================================================================
	# Network data 
	# cit-HepPh network data obtained from:
	# https://networkrepository.com/cit.php
	# Geometric median of subsamping estimates of MPLE:
	#> MPLEx$p	
	#$p
	#gwidegree(cutoff = 28093)           gwidegree.decay                 triangles 
	#             -0.8418461                 0.2365466                 1.6904300
	# ====================================================================================================================================
	# Calculate MPLE summary statistics from network G (using ERGM sufficient statistics) 
	# ====================================================================================================================================
	MPLEx			=	c(-0.8418461, 0.2365466, 1.6904300)
	MPLEterms		=	c('gwidegree', 'gwidegree.decay', 'triangle')
	names(MPLEx)	=	MPLEterms
	Nnodes 			= 	28093 # Number of nodes in network
	#
	# Price model:
	parameterNames	=	c(	'k_0', 	'p')
	d					=	length(parameterNames)
	k_0 = 1; # Prior runif(1, min = 0.9, max = 1.1) (Raynal etal.2022) will be applied later below
	p =	 .02;# Prior runif(1, min = 0, max = .20) applied later below.
	truth 	= c(k_0, p); # True data-generating model parameters to be estimated. 
	alpha 	= 	1; 	# Fixed power parameter (alpha)
	maxOutdegree = Nnodes - 1  # Fixed maximum Outdegree:  The number of nodes a new node attaches to
	#
	# is generated from a binomial distribution B(maxOutdegree, p) (Raynal etal., 2022, p.184)
	# In a directed network, the outdegree of any given node is the number of edges starting from it.
	# To analyze a real directed network dataset G0 (put in ergm package's network format), use:
	# maxOutdegree = max(degree(graph_from_edgelist(as.edgelist(G0), directed = TRUE), mode = "out"))
	# as did Raynal etal.(2022) who considered maxOutdegree = 610 for the Price model and a large network dataset.
	#	
	# For the Price model, Raynal etal.(2022) assigned a uniform prior for p centered on a matching 
	# of the first moment of the empirical network’s out-degree distribution (prior used later).
	# matchedMeanOutdegree = 	mean(degree(G), mode = "out", normalized = T)
	#
	# Raynal et al. (2022,BA) use the following statistics for the Price model.
	#outMPLEx1		=	MPLEigraph(G, c('mean_node_deg_in','var_node_deg_in'))
	#outMPLEx2		=	MPLEergm(G, c('triangles'))
	#MPLEx				=	c(outMPLEx1$MPLE[1,], outMPLEx2$MPLE[1,] )
	#seMPLE			=	sqrt(c(outMPLEx1$MPLE[2,],outMPLEx2$MPLE[2,]))
	#MPLEterms		=	c('mean_node_deg_in','var_node_deg_in','triangles')
	#names(MPLEx)	<-		MPLEterms
	#names(seMPLE)	<-		MPLEterms
	#p 					=	length(MPLEterms)
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for all network summaries from ergm package.
	# ------------------------------------------------------------------------------------------------------------------------------------
	# out	=	MPLEergm(G);  MPLE_ergm[[1]] = 	out$MPLE;	Terms_ergm[[1]]	=	out$Terms
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for selected network summaries from ergm package:
	# ------------------------------------------------------------------------------------------------------------------------------------
	# MPLE_ergm = c();		Terms_ergm = c()
	# MPLEergm(G, c('edges', 'triangles')); # Equals:#   MPLEergm(G, c('triangles', 'edges'))
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Compute MPLEs of ERGM for all network summaries from igraph package:
	# ------------------------------------------------------------------------------------------------------------------------------------
	# MPLEigraph(G) # Computation takes a while when Nnodes = 20
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for Price model
	# ====================================================================================================================================
	N				=	10000	 #	20	# Sample size of the reference table
	NnodesTable 	= 	floor(sqrt(Nnodes))	#  167
	theta.table		=	matrix(NA, nrow = N, ncol = d)
	MPLE.table		=	matrix(NA, nrow = N, ncol = length(MPLEterms))
	colnames(theta.table)	<-	parameterNames
	colnames(MPLE.table)		<-	MPLEterms
	maxOutdegreeTable = NnodesTable - 1  # Fixed maximum Outdegree:  The number of nodes a new node attaches to
	# 
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j		=	c(runif(1, min = 0.9, max = 1.1), runif(1, min = 0, max = .20))
		#
		# Simulate network dataset: 
		Gtable 		= 	simulateNetwork(Nnodes = NnodesTable, model = "Price", parameters = c(theta.j, alpha, maxOutdegreeTable))
		#
		# Compute MPLE summary statistics of the simulated data:
		outMPLEy		=	MPLEergm(Gtable, c('gwidegree(cutoff=Nnodes)','triangles'))
		MPLEy			=	outMPLEy$MPLE[1,]
		names(MPLEy)	<-	MPLEterms
		MPLEy.j		= 	MPLEy
		#
		# Raynal et al. (2022,BA) use the following statistics for the Price model.
		#outMPLEy1		=	MPLEigraph(Gtable, c('mean_node_deg_in','var_node_deg_in'))
		#outMPLEy2		=	MPLEergm(Gtable, c('triangles'))
		#MPLEy				=	c(outMPLEy1$MPLE[1,], outMPLEy2$MPLE[1,] )
		#names(MPLEy)	<-	MPLEterms
		#MPLEy.j			= 	MPLEy
		#
		# Update Reference Table
		theta.table[j,]	=	theta.j
		MPLE.table[j,]		=	MPLEy.j
		if ((j / 10) == round(j / 10))	{	flush.console()	}
	}
	# Restrict to subset of theta samples corresponding to nondegenerate netorks (finte MPLEs) 
	isFiniteMPLE		=		apply(is.finite(MPLE.table),1,all)
	theta.table 		= 		theta.table[isFiniteMPLE,]
	MPLE.table			=		MPLE.table[isFiniteMPLE,]
	N 					=		sum(isFiniteMPLE)
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Construct multivariate meta-t copula ABC posterior distribution based on DRF marginals.
	# ====================================================================================================================================
	# For each parameter (dependent variable), and set of predictors (MPLEs), 
	# Train a distributional random forest with CART splitting rule.
	postEtheta		=	matrix(NA,	nrow = 1, ncol = d)
	postSDtheta		=	matrix(NA,	nrow = 1, ncol = d)
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
	priorPDFs			=	apply(cbind(dunif(theta.table.keep[,1], min = 0.9, max = 1.1), dunif(theta.table.keep[,2], min = 0, max = .20)),1,prod)
	likelihoods			=	posteriorPDFs / priorPDFs
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
	# Save results:
	# ====================================================================================================================================
	end_time 								= 	Sys.time()
	TotalComputationTimeSimulation			=	end_time - start_time
	TotalComputationTimeAllSimulations		=	TotalComputationTimeAllSimulations + TotalComputationTimeSimulation
	#
	outputFileName	=	paste("Price cit-HepPh ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	flush.console()
	# ====================================================================================================================================