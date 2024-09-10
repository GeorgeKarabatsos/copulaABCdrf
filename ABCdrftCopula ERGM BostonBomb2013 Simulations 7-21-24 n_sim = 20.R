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
#setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Real network data\\BostonBomb2013')
#setwd('C:\\Users\\George Karabatsos\\Desktop')
setwd('C:\\Users\\user\\Desktop\\code ABCnetwork')

#rm(list = ls())
start_time 		= 	Sys.time()
TotalComputationTimeAllSimulations		=	0
set.seed(123)

# ====================================================================================================================================
# Run simulation over replicate datasets:
# ====================================================================================================================================
replicas		=	10
for (r in 6 : replicas) {	# r = 1
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
	# Simulate network dataset from multlayer ERGM model:
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
	# Summary statistics of original bostonBomb2013 dataset:
	SUMMARYterms	=	
		c(	'meanIndegreeG1', 'varIndegreeG1', 'meanOutdegreeG1', 'varOutdegreeG1', 'wClusteringCoefG1', 'assortativity_degreeG1', 'reciprocityG1',
			'meanIndegreeG2', 'varIndegreeG2', 'meanOutdegreeG2', 'varOutdegreeG2', 'wClusteringCoefG2', 'assortativity_degreeG2', 'reciprocityG2',
			'meanIndegreeG3', 'varIndegreeG3', 'meanOutdegreeG3', 'varOutdegreeG3', 'wClusteringCoefG3', 'assortativity_degreeG3', 'reciprocityG3')
	#sx				=	
	#	c(	1.35253692785	, 8474.74052782986,  1.35253692785	,  7.61861428910  ,		0.00008622621	, -0.06229008973		,  	0.00379811641,
	#		0.6141674647	, 1804.0675198631 ,  0.6141674647 	,  4.5875615155   ,		0.0004401573	, -0.0233535797		 	,  	0.0379874337,
	#		0.1590769317	,	31.5665679809 ,  0.1590769317	,  0.8667280838   ,     0.0001867074	, -0.0165706996 		, 	0.0625188215)
	# names(sx)	=	SUMMARYterms
	# Nnodes 			= 	4377184 # Number of nodes in network
	#
	# The posterior mean estimates of the multilayer ERGM, obtained from the BostonBomb 2013 dataset,
	# using the copulaABCdrf method, are shown below. 
	# They will be used as the true data-generating parameters for the simulation study. 
	#
	# Layer 1 						Posterior Mean
	# 	equalto.1.pm.0 				-5.65 
	# greaterthan.1  					-2.85 
	# mutual.min 						-3.85
	# transitiveweights.min.sum.min     0.05
	# transitiveweights.min.max.min 	-1.88
	# CMP 							-2.15
	# 
	# Layer 2 						Posterior Mean
	# equalto.1.pm.0 					-6.39
	# greaterthan.1  					-2.76
	# mutual.min  					 1.70
	# transitiveweights.min.sum.min 	 0.23
	# transitiveweights.min.max.min 	-2.48
	# CMP							-2.04
	# 
	# Layer 3 						Posterior Mean
	# equalto.1.pm.0 					-5.67
	# greaterthan.1  					-3.22
	# mutual.min 						 2.62
	# transitiveweights.min.sum.min 	-1.14
	# transitiveweights.min.max.min   	-2.15
	# CMP 							-1.34
	#
	# True data-generating parameters:
	truth 	= c(	-5.65, -2.85, -3.85,  0.05, -1.88, -2.15, 
				-6.39, -2.76,  1.70,  0.23, -2.48, -2.04,
				-5.67, -3.22,  2.62, -1.14, -2.15, -1.34)
	Nnodes	=	100
	#
	# Simulate Layer 1 network dataset: 
	#init.G 		= 	network(Nnodes, directed = TRUE, density = 0.1)
	init.G 	= 	matrix(rbinom(Nnodes ^ 2, 5, 0.005), Nnodes, Nnodes); diag(init.G) = 0;
	init.G 	= 	as.network.matrix(init.G, directed = TRUE,matrix.type="a",ignore.eval=FALSE, names.eval="tweets") # Important!
	#Gtable			= 	simulate(init.G ~ gwidegree(cutoff=Nnodes) + triangles, nsim = 1, coef = theta.j)# in 'network' format
	sx 		=	c()	# Set up vector of summary statistics of the network to be simualted.
	# 
	# Simulation given input parameters:
	G1sim	=	simulate(	init.G  ~ 
			  equalto(value = 1)
			+ greaterthan(threshold = 1)
			+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
			+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
			+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
			+ CMP, 					# Conway–Maxwell–Poisson
			response = "tweets",  reference = ~ Poisson, output = "edgelist", coef = truth[1 : 6])
	G2sim	=	simulate(	init.G  ~ 
			  equalto(value = 1)
			+ greaterthan(threshold = 1)
			+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
			+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
			+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
			+ CMP, 					# Conway–Maxwell–Poisson
			response = "tweets", reference = ~ Poisson, output = "edgelist", coef = truth[7 : 12])		
	G3sim	=	simulate(	init.G  ~ 
			  equalto(value = 1)
			+ greaterthan(threshold = 1)
			+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
			+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
			+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
			+ CMP, 					# Conway–Maxwell–Poisson
			response = "tweets", reference = ~ Poisson, output = "edgelist", coef = truth[13 : 18])	
	#
	# Find MCMLE from each of the valued networks G1sim, G2sim, and G3sim.
	MCMLE1out	=	ergm(	G1sim  ~ 
			 			 equalto(value = 1)
						+ greaterthan(threshold = 1)
						+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
						+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
						+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
						+ CMP, 					# Conway–Maxwell–Poisson
						response = "tweets",  reference = ~ Poisson)
	MCMLE2out	=	ergm(	G2sim  ~ 
						 equalto(value = 1)
						+ greaterthan(threshold = 1)
						+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
						+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
						+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
						+ CMP, 					# Conway–Maxwell–Poisson
						response = "tweets", reference = ~ Poisson)		
	MCMLE3out	=	ergm(formula = G3sim ~
						  equalto(value = 1)
						+ greaterthan(threshold = 1)
						+ mutual(form = "min") 	# "min" easiest to interpret.Most similar to binary mutuality.
						+ transitiveweights(twopath="min",combine= "sum",affect="min")# "sum" analogous to triangles
						+ transitiveweights(twopath="min",combine= "max",affect="min")# "max" analogous to transitiveties
						+ CMP, 					# Conway–Maxwell–Poisson
						response = "tweets", reference = ~ Poisson)	
	MCMLE		=	c( coef(MCMLE1out), coef(MCMLE2out), coef(MCMLE3out))
	seMCMLE		=	c( sqrt(diag(vcov(MCMLE1out))), sqrt(diag(vcov(MCMLE2out))), sqrt(diag(vcov(MCMLE3out))) )	
	#
	# Convert graphs from 'edgelist' format to 'igraph' format
	#
	G1sim	=	as.matrix(G1sim)
	if( is.matrix(G1sim[, 1:2]) == TRUE){G1	=	graph_from_edgelist(G1sim[, 1:2], 					directed = TRUE)}
	if(!is.matrix(G1sim[, 1:2]) == TRUE){G1	=	graph_from_edgelist(matrix(G1sim[, 1:2], nrow = 1), 	directed = TRUE)}
	G1		=	add_vertices(G1, nv = Nnodes - vcount(G1))
	E(G1)$weight = G1sim[, 3] # is_weighted(G1)  # [1] TRUE
	#
	G2sim	=	as.matrix(G2sim)
	if( is.matrix(G2sim[, 1:2]) == TRUE){G2	=	graph_from_edgelist(G2sim[, 1:2], 					directed = TRUE)}
	if(!is.matrix(G2sim[, 1:2]) == TRUE){G2	=	graph_from_edgelist(matrix(G2sim[, 1:2], nrow = 1), 	directed = TRUE)}
	G2		=	graph_from_edgelist(G2sim[, 1:2], directed = TRUE)
	G2		=	add_vertices(G2, nv = Nnodes - vcount(G2))
	E(G2)$weight = G2sim[, 3] # is_weighted(G2)  # [1] TRUE
	#
	G3sim	=	as.matrix(G3sim)
	if( is.matrix(G3sim[, 1:2]) == TRUE){G3	=	graph_from_edgelist(G3sim[, 1:2], 					directed = TRUE)}
	if(!is.matrix(G3sim[, 1:2]) == TRUE){G3	=	graph_from_edgelist(matrix(G3sim[, 1:2], nrow = 1), 	directed = TRUE)}
	G3		=	graph_from_edgelist(G3sim[, 1:2], directed = TRUE)
	G3		=	add_vertices(G3, nv = Nnodes - vcount(G3))
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
	sx	=	c(sx, meanIndegreeG1, 	varIndegreeG1, meanOutdegreeG1, varOutdegreeG1, 
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
	sx	=	c(sx, meanIndegreeG2, 	varIndegreeG2, meanOutdegreeG2, varOutdegreeG2, 
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
	sx	=	c(sx, meanIndegreeG3, 	varIndegreeG3, meanOutdegreeG3, varOutdegreeG3, 
					wClusteringCoefG3, 	assortativity_degreeG3, reciprocityG3)	
	# sx is the vector of summary statistics,
	# representing the simulated network dataset.
	# ====================================================================================================================================
	# Simulate Reference Table for ERGM
	# ====================================================================================================================================
	N			=	10000;	# Sample size of the reference table
	NnodesTable 	= 	20
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
		if( is.matrix(G1sim[, 1:2]) == TRUE){G1	=	graph_from_edgelist(G1sim[, 1:2], 					directed = TRUE)}
		if(!is.matrix(G1sim[, 1:2]) == TRUE){G1	=	graph_from_edgelist(matrix(G1sim[, 1:2], nrow = 1), 	directed = TRUE)}
		G1		=	add_vertices(G1, nv = NnodesTable - vcount(G1))
		E(G1)$weight = G1sim[, 3] # is_weighted(G1)  # [1] TRUE
		#
		G2sim	=	as.matrix(G2sim)
		if( is.matrix(G2sim[, 1:2]) == TRUE){G2	=	graph_from_edgelist(G2sim[, 1:2], 					directed = TRUE)}
		if(!is.matrix(G2sim[, 1:2]) == TRUE){G2	=	graph_from_edgelist(matrix(G2sim[, 1:2], nrow = 1), 	directed = TRUE)}
		G2		=	add_vertices(G2, nv = NnodesTable - vcount(G2))
		E(G2)$weight = G2sim[, 3] # is_weighted(G2)  # [1] TRUE
		#
		G3sim	=	as.matrix(G3sim)
		if( is.matrix(G3sim[, 1:2]) == TRUE){G3	=	graph_from_edgelist(G3sim[, 1:2], 					directed = TRUE)}
		if(!is.matrix(G3sim[, 1:2]) == TRUE){G3	=	graph_from_edgelist(matrix(G3sim[, 1:2], nrow = 1),	directed = TRUE)}
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
	# Compute new method for computing posterior mode and MLE (based on kernel univariate marginal posteriors):
	for (k in 1 : length(truth)){
		if (k==1){	pi.thetaRT.kernel	=	matrix(NA, nrow = N, ncol = length(truth))
					u.kernel			=	matrix(NA, nrow = N, ncol = length(truth))	}
		# For ERGM:
		XMIN = NaN; XMAX = NaN
		fitKernelRT 			= 	kde1d(theta.table[,k], xmin = XMIN, xmax = XMAX, weights = theta.DRFweights.table[,k])
		pi.thetaRT.kernel[,k]	=	dkde1d(theta.table[,k], fitKernelRT)
		u.kernel[,k]			=	pkde1d(theta.table[,k], fitKernelRT)
	}
	#"multiplied by N / (N + 1) to avoid evaluating the [copula] density at the edges of the unit square."
	# (from p. 499 of Genest & Neslehova, 2007, Astin Bulletin)
	# See also Genest et al. (1995) and Okhrin 2012, Ch.17, p.484, "Fitting High-Dimensional Copulae to Data"
	u.kernel		=	u.kernel * ( N / ( N + 1 ) )
	keep			=	apply((u.kernel > 0) & (u.kernel < 1), 1, all) & (apply(pi.thetaRT.kernel,1,prod)>0) # sum(keep)
	u.keep			=	u.kernel[keep,]
	theta.table.keep=	theta.table[keep, ]
	pi.thetaRT.keep = 	pi.thetaRT.kernel[keep,]
	# Estimate copula parameters (df, scale matrix) of Meta-t posterior distribution:
	OutCopulaFit  	=	fitStudentcopula(u.keep,fit.method = "EM-MLE", df.bounds = c(0.1, 1000), verbose = TRUE)
	post.df			=	OutCopulaFit$df
	post.scale		=	OutCopulaFit$scale

	# Take random samples from fitted student copula:
	r.u.kernel		=	rStudentcopula(n = 100000, df = post.df, scale = post.scale)
	for (k in 1 : length(truth)){
		if (k==1){	pi.thetaRT.kernel	=	matrix(NA, nrow = 100000, ncol = length(truth))
					r.theta.kernel	=	matrix(NA, nrow = 100000, ncol = length(truth))	}
		# For ERGM:
		XMIN = NaN; XMAX = NaN
		fitKernelRT 			= 	kde1d(theta.table[,k], xmin = XMIN, xmax = XMAX, weights = theta.DRFweights.table[,k])
		r.theta.kernel[,k]		=	qkde1d(r.u.kernel[,k], fitKernelRT)
		pi.thetaRT.kernel[,k]	=	dkde1d(r.theta.kernel[,k], fitKernelRT)
	}
	posteriorPDFs		=	dStudentcopula(r.u.kernel, df = post.df, scale = post.scale) * apply(pi.thetaRT.kernel, 1, prod)
	# Find posterior mode and MLE of theta using the sampled values:
	posteriorMode		=	r.theta.kernel[which.max(posteriorPDFs), 	]
	#
	# For ERGM:
	priorPDFs		=	dmvnorm(r.theta.kernel, mean = muPrior, sigma = SigmaPrior)
	likelihoods		=	posteriorPDFs / priorPDFs
	MLE				=	r.theta.kernel[which.max(likelihoods),	]
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
		postSDtheta.simulations				=	c()
		seMCMLE.simulations					=	c()
		post95cover.simulations				=	c()
		postIQRcover.simulations			=	c()
		MCMLE95cover.simulations			=	c()
		post.df.simulations					=	c()
		post.scale.simulations				=	c()
		TotalComputationTimePerSimulation	=	c()
	}

	# Update collections of simulation results:
	postEtheta.simulations					=	rbind(postEtheta.simulations, 				postEtheta							)
	postMEDtheta.simulations				=	rbind(postMEDtheta.simulations, 			postQtheta[3,]					)
	posteriorMode.simulations				=	rbind(posteriorMode.simulations, 			posteriorMode						)
	MLE.simulations							=	rbind(MLE.simulations, 						MLE									)
	MCMLE.simulations						=	rbind(MCMLE.simulations, 					MCMLE								)
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
	postSDtheta.simulations					=	rbind(postSDtheta.simulations, 			postSDtheta						)
	seMCMLE.simulations						=	rbind(seMCMLE.simulations, seMCMLE)
	post95cover.simulations					=	rbind(post95cover.simulations,	ifelse((truth > postQtheta[1,]) & (truth < postQtheta[5,]),1,0))
	postIQRcover.simulations				=	rbind(postIQRcover.simulations,ifelse((truth > postQtheta[2,]) & (truth < postQtheta[4,]),1,0))
	MCMLE95cover.simulations				=	rbind(MCMLE95cover.simulations,(truth > MCMLE - (1.96 * seMCMLE)) & (truth < MCMLE + (1.96 * seMCMLE)))
	post.df.simulations						=	rbind(post.df.simulations, rep(post.df, d))
	post.scale.simulations					=	rbind(post.scale.simulations, c(post.scale))
	TotalComputationTimePerSimulation		=	rbind(TotalComputationTimePerSimulation, rep(TotalComputationTimeSimulation, d))


	if (r == replicas){
		Mean.postEtheta.simulations			=	apply(postEtheta.simulations, 2, mean, na.rm = T);			SD.postEtheta.simulations		=	apply(postEtheta.simulations, 2, sd, na.rm = T)
		Mean.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, mean, na.rm = T);			SD.postMEDtheta.simulations		=	apply(postMEDtheta.simulations, 2, sd, na.rm = T)
		Mean.posteriorMode.simulations		=	apply(posteriorMode.simulations, 2, mean, na.rm = T);		SD.posteriorMode.simulations	=	apply(posteriorMode.simulations, 2, sd, na.rm = T)
		Mean.MLE.simulations				=	apply(MLE.simulations, 2, mean, na.rm = T);					SD.MLE.simulations				=	apply(MLE.simulations, 2, sd, na.rm = T)
		Mean.MCMLE.simulations				=	apply(MCMLE.simulations, 2, mean, na.rm = T);					SD.MCMLE.simulations				=	apply(MCMLE.simulations, 2, sd, na.rm = T)
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
		Mean.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, mean, na.rm = T);			SD.postSDtheta.simulations		=	apply(postSDtheta.simulations, 2, sd, na.rm = T)
		Mean.seMCMLE.simulations			=	apply(seMCMLE.simulations, 2, mean, na.rm = T);				SD.seMCMLE.simulations		=	apply(seMCMLE.simulations, 2, sd, na.rm = T)
		Prob.post95cover.simulations		=	apply(post95cover.simulations, 2, mean, na.rm = T)
		Prob.postIQRcover.simulations		=	apply(postIQRcover.simulations, 2, mean, na.rm = T)
		Prob.MCMLE95cover.simulations		=	apply(MCMLE95cover.simulations, 2, mean, na.rm = T)
		Mean.post.df.simulations			=	apply(post.df.simulations, 2, mean, na.rm = T);				SD.post.df.simulations					=	apply(post.df.simulations, 2, sd, na.rm = T)
		Mean.post.scale.simulations			=	apply(post.scale.simulations, 2, mean, na.rm = T);		SD.post.scale.simulations	=	apply(post.scale.simulations, 2, sd, na.rm = T)
		Mean.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, mean, na.rm = T);	SD.ComputationTimePerSimulation	=	apply(TotalComputationTimePerSimulation, 2, sd, na.rm = T)

		SimulationStudySummaryTable			=	rbind(	
		truth,
		Mean.postEtheta.simulations,		SD.postEtheta.simulations,
		Mean.postMEDtheta.simulations,		SD.postMEDtheta.simulations,
		Mean.posteriorMode.simulations,	SD.posteriorMode.simulations,
		Mean.MLE.simulations,				SD.MLE.simulations,
		Mean.MCMLE.simulations,				SD.MCMLE.simulations,
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
		Mean.postSDtheta.simulations,		SD.postSDtheta.simulations,
		Mean.seMCMLE.simulations,			SD.seMCMLE.simulations	,
		Prob.post95cover.simulations,	
		Prob.postIQRcover.simulations,
		Prob.MCMLE95cover.simulations,
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
		'Mean.postSDtheta.simulations',	'SD.postSDtheta.simulations',
		'Mean.seMCMLE.simulations',			'SD.seMCMLE.simulations'	,
		'Prob.post95cover.simulations',	
		'Prob.postIQRcover.simulations',
		'Prob.MCMLE95cover.simulations',
		'Mean.post.df.simulations',			'SD.post.df.simulations',
		'Mean.ComputationTimePerSimulation','SD.ComputationTimePerSimulation',
		'TotalComputationTimeAllSimulations'									)													
		colnames(SimulationStudySummaryTable)=	c(parameterNames)
		SimulationStudySummaryTableScale	=	rbind(	Mean.post.scale.simulations,	SD.post.scale.simulations)
		rownames(SimulationStudySummaryTableScale)=	c('Mean.post.ScaleMatrix.simulations', 'SD.post.ScaleMatrix.simulations')
	}
	# ====================================================================================================================================
	# Save results:
	# ====================================================================================================================================
	end_time 								= 	Sys.time()
	TotalComputationTimeSimulation			=	end_time - start_time
	TotalComputationTimeAllSimulations		=	TotalComputationTimeAllSimulations + TotalComputationTimeSimulation
	#
	outputFileName	=	paste("ERGM simulations BostonBomb2013 n = 300 n_sim = ", NnodesTable, " ", gsub("\\:", "_", round(end_time)), ".RData",sep='')
	save.image(file = outputFileName)
	print(paste('Replica', r, 'done'))
	flush.console()
	# ====================================================================================================================================
}