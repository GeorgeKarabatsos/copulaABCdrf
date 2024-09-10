install.packages('igraph');			library('igraph')		# For network data analysis/management
install.packages('ergm'); 			library('ergm')		# For network data analysis/management
install.packages('ergm.count'); 	library('ergm.count')# For network data analysis/management
install.packages('ergm.rank'); 	library('ergm.rank')	# For network data analysis/management
install.packages('drf');			library('drf')		# Distributional Random Forests.
install.packages('nvmix');			library('nvmix')		# For t copula.
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('ddpcr'); 			library('ddpcr')		# For quiet() 
install.packages('kde1d');		library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
#setwd('C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork')
#setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Real network data\\BostonBomb2013')

rm(list = ls())
#source('network statistics 9-14-23.R')
#source('simulateNetwork 9-14-23.R')
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Analyze real dataset using ERGM model:
# ====================================================================================================================================
	# ====================================================================================================================================
	# Calculate summary statistics from network G (using ERGM sufficient statistics) 
	# ====================================================================================================================================
	# from file: import Multilayer Network BostonBomb2013 data then calculate summaries 1-27-24
	#
	# format(c( meanIndegreeG1, varIndegreeG1, meanOutdegreeG1, varOutdegreeG1, 
	#		wClusteringCoefG1, assortativity_degreeG1, reciprocityG1), 		scientific = F)
	#  		"   1.35253692785" "8474.74052782986" "   1.35253692785" "   7.61861428910" 
	#  		"   0.00008622621"  "  -0.06229008973" "   0.00379811641"
	#
	# format(c( meanIndegreeG2, varIndegreeG2, meanOutdegreeG2, varOutdegreeG2, 
	#		wClusteringCoefG2, assortativity_degreeG2, reciprocityG2), 		scientific = F)
	# 		"	0.6141674647" "1804.0675198631" "   0.6141674647" "   4.5875615155" 
	# 		"  0.0004401573" "  -0.0233535797" "   0.0379874337"
	#
	# format(c( meanIndegreeG3, varIndegreeG3, meanOutdegreeG3, varOutdegreeG3, 
	#		wClusteringCoefG3, assortativity_degreeG3, reciprocityG3), 		scientific = F)
	#		" 0.1590769317" "31.5665679809" " 0.1590769317" " 0.8667280838"
	#		" 0.0001867074" "-0.0165706996" " 0.0625188215"
	#
	SUMMARYterms	=	
		c(	'meanIndegreeG1', 'varIndegreeG1', 'meanOutdegreeG1', 'varOutdegreeG1', 'wClusteringCoefG1', 'assortativity_degreeG1', 'reciprocityG1',
			'meanIndegreeG2', 'varIndegreeG2', 'meanOutdegreeG2', 'varOutdegreeG2', 'wClusteringCoefG2', 'assortativity_degreeG2', 'reciprocityG2',
			'meanIndegreeG3', 'varIndegreeG3', 'meanOutdegreeG3', 'varOutdegreeG3', 'wClusteringCoefG3', 'assortativity_degreeG3', 'reciprocityG3')
	sx				=	
		c(	1.35253692785	, 8474.74052782986,  1.35253692785	,  7.61861428910  ,		0.00008622621	, -0.06229008973		,  	0.00379811641,
			0.6141674647	, 1804.0675198631 ,  0.6141674647 	,  4.5875615155   ,		0.0004401573	, -0.0233535797		 	,  	0.0379874337,
			0.1590769317	,	31.5665679809 ,  0.1590769317	,  0.8667280838   ,     0.0001867074	, -0.0165706996 		, 	0.0625188215)
	names(sx)	=	SUMMARYterms
	Nnodes 			= 	4377184 # Number of nodes in network
	# ==========================================================================================
	# Set up the prior for ERG model:  
	# ==========================================================================================
	d 					=	18 
	parameterNames	=	c(	'equalto.1.pm.0', 'greaterthan.1', 'mutual.min', 
								'transitiveweights.min.sum.min', 'transitiveweights.min.max.min','CMP')
	parameterNames	=	c(	paste('G1', parameterNames,sep = '_'),
								paste('G2', parameterNames,sep = '_'),
								paste('G3', parameterNames,sep = '_'))
	muPrior 			= 	matrix(rep(0, d), nrow = 1) # mean of g prior
	SigmaPrior			=	diag(10,d, d)	# covariance matrix of g-prior.
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Simulate Reference Table for ERGM
	# ====================================================================================================================================
	N				=	10000;	# Sample size of the reference table
	NnodesTable 	= 	500
	theta.table		=	matrix(NA, nrow = N, ncol = d)
	sy.table		=	matrix(NA, nrow = N, ncol = length(SUMMARYterms))
	colnames(theta.table)	<-	parameterNames
	colnames(sy.table)		<-	SUMMARYterms
	# 
	# Generate Reference table:
	for (j in 1 : N) {	# j = 1		# For testing
		if ((j / 10) == round(j / 10))	{	print(j)	}
		# Simulate model parameters from the prior:
		theta.j		=	rmvnorm(1, mean = muPrior, sigma = SigmaPrior)
		sy.j 			=	c()	
		#
		# Simulate Layer 1 network dataset: 
		#init.G 		= 	network(NnodesTable, directed = TRUE, density = 0.1)
		init.G 		= 	matrix(rbinom(NnodesTable ^ 2, 5, 0.005), NnodesTable, NnodesTable); diag(init.G) = 0;
		init.G 		= 	as.network.matrix(init.G, directed = TRUE,matrix.type="a",ignore.eval=FALSE, names.eval="tweets") # Important!
		#Gtable			= 	simulate(init.G ~ gwidegree(cutoff=Nnodes) + triangles, nsim = 1, coef = theta.j)# in 'network' format
		# 
		# Simulation given input parameters:
		G1sim	=	simulate(	init.G  ~ 
			  equalto(value = 1)
			+ greaterthan(threshold = 1)
			+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
			+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
			+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
			+ CMP, 					# Conway–Maxwell–Poisson
			response = "tweets",  reference = ~ Poisson, output = "edgelist", coef = theta.j[1 : 6])
		G2sim	=	simulate(	init.G  ~ 
			  equalto(value = 1)
			+ greaterthan(threshold = 1)
			+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
			+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
			+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
			+ CMP, 					# Conway–Maxwell–Poisson
			response = "tweets", reference = ~ Poisson, output = "edgelist", coef = theta.j[7 : 12])		
		G3sim	=	simulate(	init.G  ~ 
			  equalto(value = 1)
			+ greaterthan(threshold = 1)
			+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
			+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
			+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
			+ CMP, 					# Conway–Maxwell–Poisson
			response = "tweets", reference = ~ Poisson, output = "edgelist", coef = theta.j[13 : 18])	
		# Convert graphs from 'edgelist' format to 'igraph' format
		#
		G1sim	=	as.matrix(G1sim)
		G1		=	graph_from_edgelist(G1sim[, 1:2], directed = TRUE)
		G1		=	add_vertices(G1, nv = NnodesTable - vcount(G1))
		E(G1)$weight = G1sim[, 3] # is_weighted(G1)  # [1] TRUE
		#
		G2sim	=	as.matrix(G2sim)
		G2		=	graph_from_edgelist(G2sim[, 1:2], directed = TRUE)
		G2		=	add_vertices(G2, nv = NnodesTable - vcount(G2))
		E(G2)$weight = G2sim[, 3] # is_weighted(G2)  # [1] TRUE
		#
		G3sim	=	as.matrix(G3sim)
		G3		=	graph_from_edgelist(G3sim[, 1:2], directed = TRUE)
		G3		=	add_vertices(G3, nv = NnodesTable - vcount(G3))
		E(G3)$weight = G3sim[, 3] # is_weighted(G3)  # [1] TRUE
		#
		# Calculate summary statistics for G1
		indegreeDistributionG1	=	degree_distribution(G1, mode = "in")
		vals					=	0 : (length(indegreeDistributionG1) - 1)
		meanIndegreeG1			=	sum(vals * indegreeDistributionG1)
		varIndegreeG1			=	sum( ((vals - meanIndegreeG1)^2) * indegreeDistributionG1)
		outdegreeDistributionG1	=	degree_distribution(G1, mode = "out")
		vals					=	0 : (length(outdegreeDistributionG1) - 1)
		meanOutdegreeG1			=	sum(vals * outdegreeDistributionG1)
		varOutdegreeG1			=	sum( ((vals - meanOutdegreeG1)^2) * outdegreeDistributionG1)
		wClusteringCoefG1		= 	transitivity(G1, type = "global")
		assortativity_degreeG1	=	assortativity_degree(G1, directed = TRUE)
		reciprocityG1			=	reciprocity(G1, mode = "default")
		sy.j	=	c(sy.j, meanIndegreeG1, 	varIndegreeG1, meanOutdegreeG1, varOutdegreeG1, 
							wClusteringCoefG1, 	assortativity_degreeG1, reciprocityG1)
		# Calculate summary statistics for G2
		indegreeDistributionG2	=	degree_distribution(G2, mode = "in")
		vals					=	0 : (length(indegreeDistributionG2) - 1)
		meanIndegreeG2			=	sum(vals * indegreeDistributionG2)
		varIndegreeG2			=	sum( ((vals - meanIndegreeG2)^2) * indegreeDistributionG2)
		outdegreeDistributionG2	=	degree_distribution(G2, mode = "out")
		vals					=	0 : (length(outdegreeDistributionG2) - 1)
		meanOutdegreeG2			=	sum(vals * outdegreeDistributionG2)
		varOutdegreeG2			=	sum( ((vals - meanOutdegreeG2)^2) * outdegreeDistributionG2)
		wClusteringCoefG2		= 	transitivity(G2, type = "global")
		assortativity_degreeG2	=	assortativity_degree(G2, directed = TRUE)
		reciprocityG2			=	reciprocity(G2, mode = "default")
		sy.j	=	c(sy.j, meanIndegreeG2, 	varIndegreeG2, meanOutdegreeG2, varOutdegreeG2, 
							wClusteringCoefG2, 	assortativity_degreeG2, reciprocityG2)		
		# Calculate summary statistics for G3
		indegreeDistributionG3	=	degree_distribution(G3, mode = "in")
		vals					=	0 : (length(indegreeDistributionG3) - 1)
		meanIndegreeG3			=	sum(vals * indegreeDistributionG3)
		varIndegreeG3			=	sum( ((vals - meanIndegreeG3)^2) * indegreeDistributionG3)
		outdegreeDistributionG3	=	degree_distribution(G3, mode = "out")
		vals					=	0 : (length(outdegreeDistributionG3) - 1)
		meanOutdegreeG3			=	sum(vals * outdegreeDistributionG3)
		varOutdegreeG3			=	sum( ((vals - meanOutdegreeG3)^2) * outdegreeDistributionG3)
		wClusteringCoefG3		= 	transitivity(G3, type = "global")
		assortativity_degreeG3	=	assortativity_degree(G3, directed = TRUE)
		reciprocityG3			=	reciprocity(G3, mode = "default")
		sy.j	=	c(sy.j, meanIndegreeG3, 	varIndegreeG3, meanOutdegreeG3, varOutdegreeG3, 
							wClusteringCoefG3, 	assortativity_degreeG3, reciprocityG3)	
		# Update Reference Table
		theta.table[j,]	=	theta.j
		sy.table[j,]	=	sy.j
		if ((j / 10) == round(j / 10))	{	flush.console()	}
	}
	# Restrict to subset of theta samples corresponding to nondegenerate netorks (finte sys) 
	isFinitesy		=		apply(is.finite(sy.table),1,all)
	theta.table 		= 		theta.table[isFinitesy,]
	sy.table			=		sy.table[isFinitesy,]
	N 					=		sum(isFinitesy)
	# ====================================================================================================================================

	# ====================================================================================================================================
	# Construct multivariate meta-t copula ABC posterior distribution based on DRF marginals.
	# ====================================================================================================================================
	# For each parameter (dependent variable), and set of predictors (sys), 
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
	# keep					=	apply((u > 0) & (u < 1), 1, all)
	# u.keep				=	u[keep,]
	# theta.table.keep	=	theta.table[keep, ]
	# pi.thetaRT.keep 	= 	pi.thetaRT[keep,]
	#
	# I found that pi.thetaRT.keep produced many zeros and u produced many zeros and ones
	# leading to all zeros for the posteriorPDFs computations, below.
	# Therefore, I tried using kernel smoothing for marginal densities and inverse cdfs, below,
	# which mitigated this issue.
	for (k in 1:d)	{
		if (k==1){		pi.thetaRT.kernel		=	matrix(NA, nrow = N, ncol = d)
						u.kernel				=	matrix(NA, nrow = N, ncol = d)			}
		fitKernelRT 				= 	kde1d(theta.table[,k], weights = theta.DRFweights.table[,k])
		pi.thetaRT.kernel[,k]	=	dkde1d(theta.table[,k], fitKernelRT)
		u.kernel[,k]				=	pkde1d(theta.table[,k], fitKernelRT)
	}
	#
	#"multiplied by N / (N + 1) to avoid evaluating the [copula] density at the edges of the unit square."
	# (from p. 499 of Genest & Neslehova, 2007, Astin Bulletin)
	# See also Genest et al. (1995) and Okhrin 2012, Ch.17, p.484, "Fitting High-Dimensional Copulae to Data"
	#u.kernel			=	u.kernel * ( N / ( N + 1 ) ) # Doesn't help according to the simulation studies.
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
	priorPDFs			=	dmvnorm(theta.table.keep, mean = muPrior, sigma = SigmaPrior)
	likelihoods		=	posteriorPDFs / priorPDFs
	MLE					=	theta.table.keep[likelihoods == max(likelihoods, na.rm = T),]
	# ------------------------------------------------------------------------------------------------------------------------------------
	# Results of Mode and MLE (using univariate kernel density estimations of the marginal posteriors:
	# ------------------------------------------------------------------------------------------------------------------------------------
	#> outModeMLE	=	round(cbind(matrix(posteriorMode,ncol=1), matrix(t(MLE),ncol=1)), 2)
	#  rownames(outModeMLE) = names(MLE);	colnames(outModeMLE) = c('Mode', 'MLE')
	#  outModeMLE
	#  > outModeMLE
	#                                    Mode   MLE
	#  G1_equalto.1.pm.0                -4.07 -4.69
	#  G1_greaterthan.1                 -0.67  1.33
	#  G1_mutual.min                    -2.43  0.12
	#  G1_transitiveweights.min.sum.min -2.04  1.24
	#  G1_transitiveweights.min.max.min -0.33 -1.30
	#  G1_CMP                            3.29  5.45
	#  G2_equalto.1.pm.0                -3.81 -5.26
	#  G2_greaterthan.1                  1.68  3.49
	#  G2_mutual.min                     0.57  0.44
	#  G2_transitiveweights.min.sum.min -2.83  3.39
	#  G2_transitiveweights.min.max.min  1.65  1.59
	#  G2_CMP                           -2.82  4.43
	#  G3_equalto.1.pm.0                -4.63 -4.28
	#  G3_greaterthan.1                  0.82 -0.24
	#  G3_mutual.min                     4.05  0.78
	#  G3_transitiveweights.min.sum.min -0.15 -6.25
	#  G3_transitiveweights.min.max.min -3.66 -0.16
	#  G3_CMP                           -1.99  2.29
	#  
	# ------------------------------------------------------------------------------------------------------------------------------------
	# It is tempting to try to optimize (using random meta-t values of theta, or using optim())
	# to get posteriorMode and MLE, using their values obtained above as starting values.
	# However doing so would still use the marginal (discrete) posterior PDF (histogram)
	# and CDF values of the original theta.jk's samples from the Reference Table.
	# BUT HOW ABOUT USING KERNEL DENSITY ESTIAMTION?
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
	outputFileName	=	paste("ERGM BostonBomb2013 ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	# ====================================================================================================================================