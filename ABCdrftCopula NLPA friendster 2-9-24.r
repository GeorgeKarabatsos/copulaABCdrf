#install.packages('igraph');# Using igraph version before January 30, 2024, for more reliable use of sample_pa()
library('igraph')		# For network data analysis/management
install.packages('ergm'); 			library('ergm')		# 	For network data analysis/management
install.packages('ergm.count'); 	library('ergm.count')# 	For network data analysis/management
install.packages('ergm.rank'); 	library('ergm.rank')	# 	For network data analysis/management
install.packages('drf');			library('drf')		# 	Distributional Random Forests.
install.packages('nvmix');			library('nvmix')		# 	For t copula.
install.packages('mvtnorm');		library('mvtnorm')	# 	Multivariate Normal and t Distributions
install.packages('MASS');			library('MASS')		# 	For mvrnorm().
install.packages('ddpcr'); 			library('ddpcr')		# 	For quiet() 
install.packages('extraDistr');	library('extraDistr')# 	For rtbinom() truncated binomial sampling.

setwd('C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork')
#setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')

rm(list = ls())
source('network statistics 9-14-23.R')
source('simulateNetwork 2-9-24.R')
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# Results obtained for friendster dataset: (using k_0 = 0)
# > SummaryTable
	# ====================================================================================================================================
	# Analyze real friendster dataset using Barabasi-Albert model:
	# ====================================================================================================================================
	# FRIENDSTER network dataset:  According to 
	# https://snap.stanford.edu/data/com-Friendster.html
	# https://networkrepository.com/soc-friendster.php
	# FRIENDSTER network dataset has the following (usable)
	# UNC statistics summaries (selected):
	# Dataset statistics
	# Nodes 	65608366
	# Edges 	1806067135
	# Maximum degree	5.2K
	# Minimum degree	1
	# Average degree	55
	# Average clustering coefficient 	0.1623
	# Diameter (longest shortest path) 	32
	# Density	8.39161e-07	# equals 1806067135 / bignchoosek(65608366, 2) approximately (using bignchoosek() from below)
	# https://stackoverflow.com/questions/40527010/r-how-can-i-calculate-large-numbers-in-n-choose-k
	# The answer to the actual question is that R cannot show numbers it cannot represent, and some of the terms in your equation are too big to represent. So it fails. However there are approximations to factorial that can be used - they work with logarithms which get big a lot slower.
	# The most famous one, Sterling's approximation, was not accurate enough, but the Ramanujan's approximation came to the rescue :)
	# ramanujan 		<- 	function(n){n*log(n) - n + log(n*(1 + 4*n*(1+2*n)))/6 + log(pi)/2}
	# nchoosek 		<- 	function(n,k){factorial(n)/(factorial(k)*factorial(n-k))} 
	# bignchoosek 	<- 	function(n,k){exp(ramanujan(n) - ramanujan(k) - ramanujan(n-k))}
	#length(triangles(sim.G)) / 3
	#transitivity(sim.G, type = "global")
	#(3*(length(triangles(sim.G)) / 3)) / nchoosek(Nnodes,3) # Does not equal to above code line.
	#count_triangles(sim.G)
	#transitivity(sim.G, type = "global") # equals transitivity(sim.G, type = "globalundirected")
	#
	# Will use the following 3 summary statistics to analyze the friendster dataset
	# using the Barabasi-Albert model (with free parameters m (count) and alpha (power) ).
	# Density	=	edge_density(sim.G, loops = FALSE) #  equals ecount(sim.G) / nchoosek(Nnodes,2) equals gsize(sim.G) / nchoosek(Nnodes,2)
	# averageClusteringCoeff = transitivity(sim.G, type = "localaverage")# equals: type = "average" and "localaverageundirected" and equals mean(transitivity(sim.G, type = "local"))
	# Diameter 	= 	diameter(sim.G, directed = FALSE)
	# ====================================================================================================================================
	# Calculate summary statistics from network G 
	# ====================================================================================================================================
	Nnodes				=	65608366
	sx					=	c(8.39161e-07, 0.1623, 32)
	#sx 				=	c( Density, averageClusteringCoeff, Diameter) 
	sTerms				=	c('Density', 'AveClusteringCoef', 'Diameter')
	names(sx)			=	sTerms
	#
	# BA model:
	parameterNames	=	c('alpha', 'p')# k_0 set to 0. alpha set to 1 (for sampling stability).
	d					=	length(parameterNames)
	# Prior (alpha, p) ~ runif(1, min = 0, max = 3) * runif(1, min = 0, max = .20) applied later below.
	maxDegree 			= 	5200  	# Fixed maximum Outdegree:  The number of nodes a new node attaches to
	meanDegree			=	55	  	# Fixed mean Outdegree:  The number of nodes a new node attaches to
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for Price model
	# ====================================================================================================================================
	N				=	10000	 #	20	# Sample size of the reference table
	NnodesTable 	= 	1000
	theta.table	=	matrix(NA, nrow = N, ncol = d)
	sy.table		=	matrix(NA, nrow = N, ncol = length(sTerms))
	colnames(theta.table)	<-	parameterNames
	colnames(sy.table)		<-	sTerms
	DegreeTable 	= 	NnodesTable - 1   # Fixed mean Degree:  The number of nodes a new node attaches to
	#
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j		=	c(runif(1, min = 0, max = 3), runif(1, min = 0, max = .20))
		# Above gives a prior sample of alpha and p.
		#
		# Simulate network dataset: 
		Gtable 		= 	simulateNetwork(	Nnodes = NnodesTable, model = "Barabasi-Albert",
												 parameters = c(0, theta.j, DegreeTable))
		# Compute summary statistics of the simulated data:
		Density		=	edge_density(Gtable, loops = FALSE) #  equals ecount(Gtable) / nchoosek(Nnodes,2) equals gsize(Gtable) / nchoosek(Nnodes,2)
		averageClusteringCoeff = transitivity(Gtable, type = "localaverage")# equals: type = "average" and "localaverageundirected" and equals mean(transitivity(Gtable, type = "local"))
		Diameter 		= 	diameter(Gtable, directed = FALSE)
		sy 				=	c( Density, averageClusteringCoeff, Diameter) 
		names(sy)		<-	sTerms
		sy.j			= 	sy
		#
		# Update Reference Table
		theta.table[j,]	=	theta.j
		sy.table[j,]		=	sy.j
		if ((j / 10) == round(j / 10))	{	flush.console()	}
	}
	#
	colnames(theta.table)	=	parameterNames
	colnames(sy.table)		=	sTerms
	#
	# Restrict to subset of theta samples corresponding to nondegenerate networks (finte sy's) 
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
	if (d == 1) {
		u.keep 			= 	as.matrix(u.keep)
		theta.table.keep	=	as.matrix(theta.table.keep)
		pi.thetaRT.keep 	= 	as.matrix(pi.thetaRT.keep)
	}
	# Estimate copula parameters (df, correlation matrix) of Meta-t posterior distribution:
	OutCopulaFit 		=	fitStudentcopula(u[keep, ], fit.method = "EM-MLE", df.bounds = c(0.1, 1000), verbose = TRUE)
	post.df			=	OutCopulaFit$df
	post.Correlations	=	OutCopulaFit$scale
	# Meta-t posterior densities of (retained) unordered theta.j's from Reference table
	posteriorPDFs		=	dStudentcopula(u.keep,df=post.df,scale=post.Correlations)*apply(pi.thetaRT.keep,1,prod)
	# Find posterior mode and MLE of theta using the sampled values:
	posteriorMode		=	theta.table.keep[posteriorPDFs == max(posteriorPDFs, na.rm = T),]
	priorPDFs			=	apply(cbind(dunif(theta.table.keep[,1], min = 0, max = 3), dunif(theta.table.keep[,2], min = 0, max = .20)),1,prod)
	#	alpha 	= 	runif(1, min = 0,   max = 3		)
	#	p 		=	runif(1, min = 0  , max = .20 	)
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
	# Save results:
	# ====================================================================================================================================
	end_time 									= 	Sys.time()
	TotalComputationTimeSimulation			=	end_time - start_time
	TotalComputationTimeAllSimulations		=	TotalComputationTimeAllSimulations + TotalComputationTimeSimulation
	#
	outputFileName	=	paste("Barabasi-Albert friendster ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	flush.console()
	# ====================================================================================================================================