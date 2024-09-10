install.packages('igraph');# Using igraph version before January 30, 2024, for more reliable use of sample_pa()
library('igraph')		# For network data analysis/management
install.packages('ergm'); 			library('ergm')		# For network data analysis/management
install.packages('ergm.count'); 	library('ergm.count')# For network data analysis/management
install.packages('ergm.rank'); 	library('ergm.rank')	# For network data analysis/management
install.packages('drf');			library('drf')		# Distributional Random Forests.
install.packages('nvmix');			library('nvmix')		# For t copula.
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('ddpcr'); 			library('ddpcr')		# For quiet() 
install.packages('extraDistr');	library('extraDistr')# 	For rtbinom() truncated binomial sampling.
#setwd('C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork')
#setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork')

rm(list = ls())
source('network statistics 9-14-23.R')
source('simulateNetwork 2-9-24.R')
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Run simulation over replicate datasets:
# ====================================================================================================================================
replicas		=	10 # 
for (r in 1 : replicas) {	# r = 1
	# ====================================================================================================================================
	# Simulate network dataset from Price model (network dataset to be analyzed using ABC)
	# ====================================================================================================================================
	# Network data 
	# Price model:
	Nnodes 			= 50 # Number of nodes in network (to simulate)
	parameterNames	=	c('alpha', 'p')
	d				=	length(parameterNames)
	k_0 		= 	0; 	
	p			=	.02;
	alpha 		= 	1.2; 	# Fixed power parameter (alpha)
	truth 		= 	c(alpha, p); # True data-generating model parameters to be estimated. 
	maxDegree 	= 	Nnodes - 1  # Fixed maximum Outdegree:  The number of nodes a new node attaches to
	# is generated from a binomial distribution B(maxOutdegree, p) (Raynal etal., 2022, p.184)
	# In a directed network, the outdegree of any given node is the number of edges starting from it.
	# To analyze a real directed network dataset G0 (put in ergm package's network format), use:
	# maxOutdegree = max(degree(graph_from_edgelist(as.edgelist(G0), directed = TRUE), mode = "out"))
	# as did Raynal etal.(2022) who considered maxOutdegree = 610 for the Price model and a large network dataset.
	#	
	G 	= 	simulateNetwork(Nnodes = Nnodes, model = "Barabasi-Albert", parameters = c(0, truth, maxDegree))
	#
	# ====================================================================================================================================
	# Calculate summary statistics from network G (using ERGM sufficient statistics) 
	# ====================================================================================================================================
	sTerms		=	c('Density', 'AveClusteringCoef', 'Diameter')
	Density		=	edge_density(G, loops = FALSE) #  equals ecount(Gtable) / nchoosek(Nnodes,2) equals gsize(Gtable) / nchoosek(Nnodes,2)
	averageClusteringCoeff = transitivity(G, type = "localaverage")# equals: type = "average" and "localaverageundirected" and equals mean(transitivity(Gtable, type = "local"))
	Diameter 	= 	diameter(G, directed = FALSE)
	sx 			=	c( Density, averageClusteringCoeff, Diameter) 
	names(sx)	=	sTerms
	p 			=	length(sTerms)
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for Price model
	# ====================================================================================================================================
	N				=	10000 #	20		# Sample size of the reference table
	NnodesTable 	= 	floor(Nnodes * 3/4);
	theta.table		=	matrix(NA, nrow = N, ncol = d)
	sy.table		=	matrix(NA, nrow = N, ncol = length(sTerms))
	colnames(theta.table)	<-	parameterNames
	colnames(sy.table)		<-	sTerms
	DegreeTable = NnodesTable - 1  # Fixed maximum Outdegree:  The number of nodes a new node attaches to
	# 
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j		=	c(runif(1, min = 0, max = 3), runif(1, min = 0, max = .20))
		#
		# Simulate network dataset: 
		Gtable 		= 	simulateNetwork(Nnodes = NnodesTable, model = "Barabasi-Albert", parameters = c(0, theta.j, DegreeTable))
		#
		# Compute summary statistics of the simulated data:
		Density.y	=	edge_density(Gtable, loops = FALSE) #  equals ecount(Gtable) / nchoosek(Nnodes,2) equals gsize(Gtable) / nchoosek(Nnodes,2)
		averageClusteringCoeff.y = transitivity(Gtable, type = "localaverage")# equals: type = "average" and "localaverageundirected" and equals mean(transitivity(Gtable, type = "local"))
		Diameter.y 	= 	diameter(Gtable, directed = FALSE)
		sy 			=	c( Density.y, averageClusteringCoeff.y, Diameter.y) 
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
	# Estimate copula parameters (df, correlation matrix) of Meta-t posterior distribution:
	OutCopulaFit 		=	fitStudentcopula(u[keep, ], fit.method = "EM-MLE", df.bounds = c(0.1, 1000), verbose = TRUE)
	post.df			=	OutCopulaFit$df
	post.Correlations	=	OutCopulaFit$scale
	# Meta-t posterior densities of (retained) unordered theta.j's from Reference table
	posteriorPDFs		=	dStudentcopula(u.keep,df=post.df,scale=post.Correlations)*apply(pi.thetaRT.keep,1,prod)
	# Find posterior mode and MLE of theta using the sampled values:
	posteriorMode		=	theta.table.keep[posteriorPDFs == max(posteriorPDFs, na.rm = T),]
	priorPDFs			=	apply(cbind(dunif(theta.table.keep[,1], min = 0, max = 3), dunif(theta.table.keep[,2], min = 0, max = .20)),1,prod)
	likelihoods			=	posteriorPDFs / priorPDFs
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
		postMEDtheta.simulations			=	c()
		posteriorMode.simulations			=	c()
		MLE.simulations						=	c()
		L1error.postEtheta.simulations		=	c()
		L2error.postEtheta.simulations		=	c()
		L1error.postMEDtheta.simulations	=	c()
		L2error.postMEDtheta.simulations	=	c()
		L1error.posteriorMode.simulations	=	c()
		L2error.posteriorMode.simulations	=	c()
		L1error.MLE.simulations				=	c()
		L2error.MLE.simulations				=	c()
		postSDtheta.simulations				=	c()
		post95cover.simulations				=	c()
		postIQRcover.simulations			=	c()
		post.df.simulations					=	c()
		post.Correlations.simulations		=	c()
		TotalComputationTimePerSimulation	=	c()
	}

	# Update collections of simulation results:
	postEtheta.simulations					=	rbind(postEtheta.simulations, 				postEtheta							)
	postMEDtheta.simulations				=	rbind(postMEDtheta.simulations, 			postQtheta[3,]					)
	posteriorMode.simulations				=	rbind(posteriorMode.simulations, 			posteriorMode						)
	MLE.simulations							=	rbind(MLE.simulations, 						MLE									)
	L1error.postEtheta.simulations		=	rbind(L1error.postEtheta.simulations, 	abs(postEtheta - truth) 		)
	L2error.postEtheta.simulations		=	rbind(L2error.postEtheta.simulations,    (postEtheta - truth)^2 			)
	L1error.postMEDtheta.simulations		=	rbind(L1error.postMEDtheta.simulations,	abs(postQtheta[3,] - truth)		)
	L2error.postMEDtheta.simulations		=	rbind(L2error.postMEDtheta.simulations,  (postQtheta[3,] - truth)^2 		)
	L1error.posteriorMode.simulations	=	rbind(L1error.posteriorMode.simulations, abs(posteriorMode - truth) 		)
	L2error.posteriorMode.simulations	=	rbind(L2error.posteriorMode.simulations,    (posteriorMode - truth)^2 	)
	L1error.MLE.simulations					=	rbind(L1error.MLE.simulations, 			abs(MLE - truth) 					)
	L2error.MLE.simulations					=	rbind(L2error.MLE.simulations,    				(MLE - truth)^2 				)
	postSDtheta.simulations					=	rbind(postSDtheta.simulations, 			postSDtheta						)
	post95cover.simulations					=	rbind(post95cover.simulations,	ifelse((truth > postQtheta[1,]) & (truth < postQtheta[5,]),1,0))
	postIQRcover.simulations				=	rbind(postIQRcover.simulations,ifelse((truth > postQtheta[2,]) & (truth < postQtheta[4,]),1,0))
	post.df.simulations						=	rbind(post.df.simulations, rep(post.df, d))
	post.Correlations.simulations			=	rbind(post.Correlations.simulations, c(post.Correlations))
	TotalComputationTimePerSimulation		=	rbind(TotalComputationTimePerSimulation, rep(TotalComputationTimeSimulation, d))


	if (r == replicas){
		Mean.postEtheta.simulations			=	apply(postEtheta.simulations, 2, mean, na.rm = T);			SD.postEtheta.simulations		=	apply(postEtheta.simulations, 2, sd, na.rm = T)
		Mean.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, mean, na.rm = T);			SD.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, sd, na.rm = T)
		Mean.posteriorMode.simulations		=	apply(posteriorMode.simulations, 2, mean, na.rm = T);		SD.posteriorMode.simulations	=	apply(posteriorMode.simulations, 2, sd, na.rm = T)
		Mean.MLE.simulations					=	apply(MLE.simulations, 2, mean, na.rm = T);					SD.MLE.simulations				=	apply(MLE.simulations, 2, sd, na.rm = T)
		MAE.postEtheta.simulations			=	apply(L1error.postEtheta.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.postEtheta.simulations			=	apply(L2error.postEtheta.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.postMEDtheta.simulations		=	apply(L1error.postMEDtheta.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.postMEDtheta.simulations		=	apply(L2error.postMEDtheta.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.posteriorMode.simulations		=	apply(L1error.posteriorMode.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.posteriorMode.simulations		=	apply(L2error.posteriorMode.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		MAE.MLE.simulations					=	apply(L1error.MLE.simulations, 2, mean, na.rm = T) # Mean Absolute Error (MAE)
		MSE.MLE.simulations					=	apply(L2error.MLE.simulations, 2, mean, na.rm = T) # Mean Squared Error (MSE)
		Mean.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, mean, na.rm = T);			SD.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, sd, na.rm = T)
		Prob.post95cover.simulations		=	apply(post95cover.simulations, 2, mean, na.rm = T)
		Prob.postIQRcover.simulations		=	apply(postIQRcover.simulations, 2, mean, na.rm = T)
		Mean.post.df.simulations			=	apply(post.df.simulations, 2, mean, na.rm = T);				SD.post.df.simulations					=	apply(post.df.simulations, 2, sd, na.rm = T)
		Mean.post.Correlations.simulations	=	apply(post.Correlations.simulations, 2, mean, na.rm = T);		SD.post.Correlations.simulations	=	apply(post.Correlations.simulations, 2, sd, na.rm = T)
		Mean.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, mean, na.rm = T);	SD.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, sd, na.rm = T)

		SimulationStudySummaryTable			=	rbind(	
		truth,
		Mean.postEtheta.simulations,		SD.postEtheta.simulations,
		Mean.postMEDtheta.simulations,		SD.postMEDtheta.simulations,
		Mean.posteriorMode.simulations,	SD.posteriorMode.simulations,
		Mean.MLE.simulations,				SD.MLE.simulations,
		MAE.postEtheta.simulations,
		MSE.postEtheta.simulations,
		MAE.postMEDtheta.simulations,
		MSE.postMEDtheta.simulations,
		MAE.posteriorMode.simulations,
		MSE.posteriorMode.simulations,
		MAE.MLE.simulations,
		MSE.MLE.simulations,
		Mean.postSDtheta.simulations,		SD.postSDtheta.simulations,
		Prob.post95cover.simulations,	
		Prob.postIQRcover.simulations,
		Mean.post.df.simulations,			SD.post.df.simulations,
		Mean.ComputationTimePerSimulation,	SD.ComputationTimePerSimulation,
		sum(TotalComputationTimeAllSimulations)								)
	
		rownames(SimulationStudySummaryTable)=	c(		
		'truth',
		'Mean.postEtheta.simulations',		'SD.postEtheta.simulations',
		'Mean.postMEDtheta.simulations',	'SD.postMEDtheta.simulations',
		'Mean.posteriorMode.simulations',	'SD.posteriorMode.simulations',
		'Mean.MLE.simulations',				'SD.MLE.simulations',
		'MAE.postEtheta.simulations',
		'MSE.postEtheta.simulations',
		'MAE.postMEDtheta.simulations',
		'MSE.postMEDtheta.simulations',
		'MAE.posteriorMode.simulations',
		'MSE.posteriorMode.simulations',
		'MAE.MLE.simulations',
		'MSE.MLE.simulations',
		'Mean.postSDtheta.simulations',	'SD.postSDtheta.simulations',
		'Prob.post95cover.simulations',	
		'Prob.postIQRcover.simulations',
		'Mean.post.df.simulations',			'SD.post.df.simulations',
		'Mean.ComputationTimePerSimulation','SD.ComputationTimePerSimulation',
		'TotalComputationTimeAllSimulations'									)													
		colnames(SimulationStudySummaryTable)=	c(parameterNames)
		SimulationStudySummaryTableCorrel	=	rbind(	Mean.post.Correlations.simulations,	SD.post.Correlations.simulations)
		rownames(SimulationStudySummaryTableCorrel)=	c('Mean.post.CorrelationMatrix.simulations', 'SD.post.CorrelationsMatrix.simulations')
	}
	outputFileName	=	paste("NLPA simul n = 50 NnodesTable = ", NnodesTable, " ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	print(paste('Replica', r, 'done'))
	flush.console()
	# ====================================================================================================================================
}