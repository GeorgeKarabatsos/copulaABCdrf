simulateNetwork 	= 	function(	Nnodes = NULL, 	model, 	parameters, 
									ERGMformula = NULL, isDirected = FALSE, 
									Nnodes1 = NULL, Nnodes2 = NULL, Time = NULL		)

# 	This code file (function) simulates a network from a network model,
# 	given inputs of specified model parameters, and other network parameters.
#   On 2-8-24, I updated the Barabasi-Albert sampling to be parameterized by p instead of m.
#
#	INPUTS:
#	Nnodes			The number of nodes (vertices) used to simulate the network.
#
#	model			The network model from which to simulate a network.
#					Currently, the available choices for input model are:
#					model = 'ERGM' or 'Price' or 'Barabasi-Albert' or 'DMC' or 'DMR' or 'Watts-Strogatz'
#
#	parameters		The parameters of the network model from which to simulate a network.
#
#	ERGMformula		The formula (with network terms) used to specify the ERGM network model.
#					This input is only required for input:	model = 'ERGM' 
#					For any other model (input), the ERGMformula input is ignored.					
#			 
#	isDirected		isDirected = FALSE (the default) if simulating an undirected network.
#					isDirected = TRUE if simulating an undirected network.
#					The isDirected input is only required for input:	model = 'ERGM'
#					For any other model (input), the isDirected input is ignored, 
#					because the other model choices:
#					model = 'Price'  			automatically outputs a simulated directed network graph;
#					model = 'Barabasi-Albert'  	automatically outputs a simulated undirected network graph;
#					model = 'DMC'  				automatically outputs a simulated undirected network graph;
#					model = 'DMR'  				automatically outputs a simulated undirected network graph;
# 					model = 'Watts-Strogatz'	automatically outputs a simulated undirected network graph;
# 					model = 'Kretzschmar-Morris'automatically outputs a simulated undirected bipartite network graph.
#
#	Nnodes1,Nnodes2 The number of nodes for the two groups in a bipartite undirected network graph, respectively.
#					The Nnodes1 and Nnodes2 inputs are only required for input:	model = 'Kretzschmar-Morris'
#					For any other model (input), the Nnodes1 and Nnodes2 inputs are ignored, 
# 					since all other model choices are for non-bipartite networks.
#					However, in the future, this code function file can be modified to simulate 
#					bipartite undrected networks from the ERGM model.
#	
#	Time			simulation time of the network. 
#					The Time input is only required for input:	model = 'Kretzschmar-Morris'
#					For any other model (input), the Time input is ignored.
#
#	OUTPUT:
#	sim.G 			The simulated network graph, either in the 'ergm' or 'igraph' R package format,
#					depending on the choice of model input.



{
	if (model == "ERGM")	{
		init.G = 	matrix(rbinom(Nnodes ^ 2, 1, 0.005), Nnodes, Nnodes); diag(init.G) = 0;
		init.G = 	as.network.matrix(init.G, directed = isDirected)
		# Alternatively:
		#	init.G 		= 	network(Nnodes, directed = isDirected, density = 0.1)
		#	# density 	= ratio of observed edges to the number of possible edges for the given network. 
		ERGMformula_init 	= 	unlist(strsplit(as.character(ERGMformula),split = " "))
		ERGMformula_init 	= 	ERGMformula_init[ERGMformula_init != "~"]
		ERGMformula_init 	= 	ERGMformula_init[ERGMformula_init != "G"]
		ERGMformula_init 	= 	ERGMformula_init[ERGMformula_init != "offset(edges)"]
		ERGMformula_init 	= 	ERGMformula_init[ERGMformula_init != "+"]
		ERGMformula_init 	=	as.formula(paste('init.G ~', paste(ERGMformula_init, collapse = " + ")))
		sim.G 	=	simulate(	ERGMformula_init, nsim = 1,
								coef = c( log(1/Nnodes) + parameters[1], parameters[2:length(parameters)]))
		sim.G 	= 	set.network.attribute(sim.G, 'n', Nnodes)# list.network.attributes(sim.G)
		# Above: Simulating based on offset log(1/Nnodes), using Supplementary Code from Schmid & Hunter (2023)
		# If you are using reciprocity statistic for a ERGM for a directed graph, 
		# you need to make a slightly different offset adjustment (see Krivitsky et al 2015 & Stewart etal. 2019).
		# Alternatively, set nsim > 1 as in Synthetic likelihood method; say, nsim = 100. Then:
		#sim.Gt 	= 	set.network.attribute(sim.G[[100]], 'n', Nnodes)# sim.Gt is the 100th network simulated from ERGM
	}



	if ((model == "Price") | (model == "Barabasi-Albert"))	{
		# Simulate from a Preferential Attachment mechanistic network model. From igraph package of R.
		# Preferential attachment is a family of simple stochastic algorithms for building a graph. 
		# Variants include the Barabási-Albert model (directed = FALSE) and the Price model (directed = TRUE).
		# https://en.wikipedia.org/wiki/Price%27s_model
		# https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model 
		# This is a simple stochastic algorithm to generate a graph. 
		# (For e.g., consider a network of papers (i.e., the vertices or nodes) where a directed edge 
		# between any two papers means that one paper cited another paper)
		# It is a discrete time step model and in each time step a single vertex (e.g., paper) is added.
		# We start with a single vertex (paper) and no edges (citations) in the first time step. 
		# Then we add one vertex (paper) in each time step and the new vertex (new paper) 
		# initiates some edges (cites) to old vertices (old papers).
		# The probability that an old vertex (old paper) is chosen (cited) is given by 
		# P[i] ~ (k_i)^alpha + a, 
		# where k_i is the indegree of vertex i (more precisely the number of adjacent edges of i 
		# which were not initiated by i itself) 
		# and 'alpha' and 'a' (often denoted k_0) are parameters given by the power and zero.appeal arguments. 
		# The number of edges initiated in a time step is given by the m, out.dist and out.seq arguments:
		# m 			Numeric constant, the number of edges to add in each time step.
		# 					This argument is only used if both out.dist and out.seq are omitted or NULL.
		# out.dist		Numeric vector, the distribution of the number of edges to add in each time step.
		#					This argument is only used if the out.seq argument is omitted or NULL.
		# out.seq 		Numeric vector giving the number of edges to add in each time step.
		#					(* used by Raynal etal 2022, Section 5, for Price's (1965) model) 
		#					Its first element is ignored as no edges are added in the first time step.
		# out.pref 	Logical, if true the total degree is used for calculating the citation probability,
		#					otherwise the in-degree is used.[*The Price model uses in-degree (out.pref = FALSE)].]
		# From Raynal et al. (2022, BA, p.184):"The number of nodes a new node attaches to is generated by a 
		# binomial distribution B(610,p), where the upper bound of 610 is motivated by the maximum out-degree
		# in the empirical network [Scientific citation network dataset from the American Physical Society, 2019]. 
		# The prior distribution of (k_0, p) is uniform over the reactangle [0.9,1.1] x [0.019, 0.021]
		# [where k_0 denotes the constant parameter of the Price model, here, k_0 is the 'a' parameter]"
		# Raynal considered the Price model defined by alpha = 1 and parameters (k_0 (= a), p).
		# If out.seq is given and not NULL then it gives the number of edges to add in a vector, the first
		# element is ignored, the second is the number of edges to add in the second time step and so on.
		# If out.seq is not given or null and out.dist is given and not NULL then it is used as a discrete
		# distribution to generate the number of edges in each time step. Its first element is the probability 
		# that no edges will be added, the second is the probability that one edge is added, etc. (out.dist 
		# does not need to sum up to one, it normalized automatically.) out.dist should contain non-negative
		# numbers and at east one element should be positive.
		# Simulate a network from the Price (1965) preferential attachment model:
		# Price model parameters (names) are 'k_0' and 'p'.
		# parameters = c(k_0, p, alpha, maxOutdegree)
		isDirected	=	model == "Price" # Alternatively, Barabasi-Albert model is for an undirected network.
		if (model == "Price")	{
			k_0			=	parameters[1]	# Constant parameter (k_0, i.e., 'a'). Price (1965) proposed k_0 = 1.
			p 			=	parameters[2]	# p =.02 # consistent with Raynal etal.'s (2022, BA, p.184) prior for p (see above).
			alpha 		=	parameters[3]	# Power parameter (fixed); Price's (1965) model assumes alpha = 1.
			maxOutdegree	=	parameters[4]	# the maximum node outdegree for Price model (directed network);
			# Renormalize to account for simulating network size Nnodes:
			maxOutdegree	=	floor((maxOutdegree/(Nnodes - 1)) * (Nnodes-1)) 
			# Above, N0 is the number of nodes of original graph analyzed.
			# The number of nodes a new node attaches to is generated from a binomial distribution 
			# B(maxOutdegree, p) (as in Raynal etal., 2022, p.184):
			# "The number of nodes a new node attaches to is generated from a binomial distribution B(610, p),
			# where the upper bound of 610 is motivated by the maximum out-degree in the empirical 
			# network [being the empirical citation network from the American Physical Society (2019), 
			# containing n_o = 597,819 articles starting 1893 with over 7 million citations 
			# (corresponding to directed edges).]."
			sim.G 		= 	sample_pa(	Nnodes, power = alpha, zero.appeal = k_0, directed = TRUE,
										out.seq 	= 	rbinom(Nnodes, size = maxOutdegree, p),  
										out.dist 	= 	NULL, out.pref = FALSE, start.graph = NULL, m = NULL, algorithm = "psumtree")
		}
		#if (model == "Barabasi-Albert"){
		#	k_0			=	parameters[1]# Constant parameter (k_0, i.e., 'a'). (unknown parameter) 
			# "k0 is allowed to be negative—it can fall anywhere in the range −m < k_0 < Inf 
			# and the probability of attachment will be positive."(p.219 of Newman,2003,SIAM Review,45,No.2,pp.167–256)
		#	alpha 		=	parameters[2]# Power parameter 
		#	m 			=	parameters[3]# The number of edges to add in each time step. (fixed positive integer). 
		#	sim.G 		= 	sample_pa(	Nnodes, power = alpha, zero.appeal = k_0, directed = FALSE, m = m, 
		#							out.dist = NULL, out.seq = NULL,out.pref = FALSE, start.graph = NULL, algorithm = "psumtree")
		#}
		if (model == "Barabasi-Albert"){
			k_0			=	parameters[1]# Constant parameter (k_0, i.e., 'a'). (unknown parameter) 
			# "k0 is allowed to be negative—it can fall anywhere in the range −m < k_0 < Inf 
			# and the probability of attachment will be positive."(p.219 of Newman,2003,SIAM Review,45,No.2,pp.167–256)
			alpha 		=	parameters[2]# Power parameter 
			p 			=	parameters[3]# The number of edges to add in each time step. (fixed positive integer).
			Degree		=	parameters[4]# the maximum node degree for BA model (undirected network);
			# The BA model is a specific case of the more general non-linear preferential attachment (NLPA) model
          	# (Krapivsky 2000). The NLPA algorithm is identical to the Barabasi-Albert model with the 
			# attachment probability based on a positive exponent alpha, and alpha = 1 
			# defines the Barabasi-Albert model: https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
			if ((k_0 != 0)){	OUT.seq 	= 	rbinom (Nnodes, size = Degree, p)								}
			if ((k_0 == 0)){	OUT.seq 	= 	rtbinom(Nnodes, size = Degree, prob = p, a = 0, b = Inf)	}
			# Above line truncates out event of 0. Needed when k_0 = 0 according to email communication 
			# with Szabolcs Horvát on February 9, 2024.
			sim.G 		= 	sample_pa(	Nnodes, power = alpha, zero.appeal = k_0, directed = FALSE, m = NULL, 
										out.seq 	= 	OUT.seq, 
										out.dist 	= 	NULL, out.pref = FALSE, start.graph = NULL, algorithm = "psumtree")
		}
	}



	if (model == "DMC"){
		# Simulate an undirected network graph from Duplication-Mutation-Complementation (DMC) model 
		# with probability parameters = (q_m, q_c) (DMC model of Vazquez et al. 2003; Vazquez 2003).
		# From Chen & Onnela (2019, SciRep, pp.4-5):
		# Duplication-divergence models. Duplication-divergence models are a popular class of models 
		# used for protein-protein interaction networks. Some examples include the 
		# duplication-mutation-complementation (DMC) (Vázquez et al. 2003, Complexus) and 
		# duplication-mutation-random mutation models (DMR) (Solé et al. 2002, Adv. Complex Syst.; 
		# Pastor-Satorras etal. 2003, J. Theor. Biology). 
		# Given a seed network, both DMC and DMR models grow the network according to their
		# respective generative mechanisms until the requisite number of nodes, n, is reached.
		# In both the DMC and DMR models, a new node is first added at the beginning of each
		# step in network generation. An existing node is chosen uniformly at random for duplication,
		# and an edge is then added between the new node and each neighbor of the chosen node.
		# [The above process represents Steps 1 and 2(a) mentioned below].
		# After this, the two models diverge. 
		# For DMC, for each neighbor of the chosen node, one of the edge between the chosen node
		# and the neighbor or the edge between the new node and the neighbor is randomly chosen
		# and then removed with probability q_mod (or q_m). The step is concluded by adding 
     	# an edge between the chosen node and the new node with probability q_con (or q_c).
		# For DMR, each edge connected to the new node is removed independently with probability q_del.
		# The step concludes by adding an edge between the new node and any existing node at 
		# the start of step t with probability qnew/n(t), where n(t) is the number of nodes 
		# in the network at the start of step t.
		#
		# Example for testing 
		# Nnodes =	30;  q_m =	 4 / 11;  q_c = 5 / 11;
		q_m 	=	parameters[1]
		q_c 	=	parameters[2]
		# The algorithm for generating a DMC graph with n nodes is as follows 
		# (from Larson and Onnela, 2023, pp.1-2; accepted by Phys. Rev. E in July 25 2023):
		for (t_ in 1 : Nnodes) {
			if (t_ == 1){
				# Step 1: Begin with a seed graph containing a single node:
				sim.G 	= 	make_empty_graph(n = 0, directed = FALSE)
				j 		=	1
				sim.G 	= 	add_vertices(sim.G, 1)# V(sim.G) lists nodes/vertices in current network
			} 
			if (t_ > 1){
				# Step 2(a): Duplication:
				#	i.	 Select anchor node uniformly at random:
				anchorNode	=	sample(length(V(sim.G)), 1)# V(sim.G) lists nodes/vertices in current network
				# 	ii.	 Add a new node to the current network graph:
				j 		=	j + 1 # Updated number of nodes at time t_ (also the ID of the new node added)
				sim.G 	= 	add_vertices(sim.G, 1)# V(sim.G) lists nodes/vertices in current network
				# 	iii. Connect the new node to the anchor node's neighbors:
				anchorNodeNeighbors 	=	neighbors(sim.G, anchorNode, mode = "all")
				anchorNodeNeighbors 	=	as.numeric(anchorNodeNeighbors)
				newEdges	=	c(t(matrix(c(rep(j, length(anchorNodeNeighbors)), c(anchorNodeNeighbors)), ncol = 2)))
				sim.G 		= 	add_edges(sim.G, newEdges) # E(sim.G) lists edges in current network
				sim.G		=	simplify(sim.G) # Removes any multiple edges. See: E(sim.G)
				#
				# Step 2(b): Mutation:
				# 	i. "Modify" each of the anchor node's neighbors independently and with probability q_m.
				IndModifyAnchorNodeNeighbors	= 	runif(length(anchorNodeNeighbors), 0, 1) <= q_m
				IndModifyAnchorNodeNeighbors	=	anchorNodeNeighbors[IndModifyAnchorNodeNeighbors]
				# 	ii.	If a neighbor is modified, remove either the neighbor-anchor node edge 
				#     	or the neighbor-new node edge. Which edge is lost is determined by the flip of a fair coin.
				neighborAnchorNodeEdges 	= 	c(t(matrix(c(rep(anchorNode,length(IndModifyAnchorNodeNeighbors)), c(IndModifyAnchorNodeNeighbors)), ncol = 2)))
				neighborNewNodeEdges 		= 	c(t(matrix(c(rep(j, length(IndModifyAnchorNodeNeighbors)), c(IndModifyAnchorNodeNeighbors)), ncol = 2)))
				neighborAnchorNodeEdgesIDs 	= 	get.edge.ids(sim.G, neighborAnchorNodeEdges)
				neighborNewNodeEdgesIDs		= 	get.edge.ids(sim.G, neighborNewNodeEdges)
				z 								= 	ifelse(runif(length(neighborAnchorNodeEdgesIDs)) <= 1 / 2, 1, 2)
				z[(neighborAnchorNodeEdgesIDs       	+ 	neighborNewNodeEdgesIDs) == 0	]	=	0
				z[(neighborAnchorNodeEdgesIDs > 0	)  	& 	(neighborNewNodeEdgesIDs == 0)	]	=	1
				z[(neighborAnchorNodeEdgesIDs == 0	) 	& 	(neighborNewNodeEdgesIDs >  0)	]	=	2
				deleteEdgesIDs 				=	((z == 1) * neighborAnchorNodeEdgesIDs) + ((z == 2) * neighborNewNodeEdgesIDs)
				deleteEdgesIDs 				=	deleteEdgesIDs[deleteEdgesIDs > 0]
				if(length(deleteEdgesIDs) > 0){sim.G = delete_edges(sim.G, deleteEdgesIDs)}#E(sim.G) lists edges of current network graph
				# Step 2(c) Complementation
				# 	i. Connect the new node to the anchor node with probability q_c.
				if (runif(1) < q_c) {sim.G 	=	add_edges(sim.G, c(j, anchorNode))}# E(sim.G) lists nodes/vertices in current network
			}
		}
	}



	if (model == "DMR"){
		# Simulate an undirected network graph from the Duplication-Mutation-Random mutation models (DMR) 
		# (Solé et al. 2002; Pastor-Satorras etal. 2003) with parameters = (q_d, q_n).
		# From Chen & Onnela (2019, SciRep, pp.4-5):
		# Duplication-divergence models. Duplication-divergence models are a popular class of models 
		# used for protein-protein interaction networks. Some examples include the 
		# duplication-mutation-complementation (DMC) (Vázquez et al. 2003, Complexus) and 
		# duplication-mutation-random mutation models (DMR) (Solé et al. 2002, Adv. Complex Syst.; 
		# Pastor-Satorras etal. 2003, J. Theor. Biology). 
		# Given a seed network, both DMC and DMR models grow the network according to their
		# respective generative mechanisms until the requisite number of nodes, n, is reached.
		#
		# In both the DMC and DMR models, a new node is first added at the beginning of each
		# step in network generation. An existing node is chosen uniformly at random for duplication,
		# and an edge is then added between the new node and each neighbor of the chosen node.
		# [The above process represents Steps 1 and 2(a) mentioned below].
		# After this, the two models diverge. 
		#
		# For DMC, for each neighbor of the chosen node, one of the edge between the chosen node
		# and the neighbor or the edge between the new node and the neighbor is randomly chosen
		# and then removed with probability q_mod (or q_m). The step is concluded by adding 
     	# an edge between the chosen node and the new node with probability q_con (or q_c).
		#
		# For DMR, each edge connected to the new node is removed independently with probability
		# q_del (or q_d). The step concludes by adding an edge between the new node and any
		# existing node at the start of step t with probability qnew/n(t) (or q_n / n(t)), 
		# where n(t) is the number of nodes in the network at the start of step t.
		#
		# Example for testing 
		# Nnodes =	30;  q_d =	 4 / 11;  q_n = 5 / 11;
		q_d 	=	parameters[1]
		q_n 	=	parameters[2]
		# The algorithm for generating a DMC graph with n nodes is as follows 
		# (from Larson and Onnela, 2023, pp.1-2; accepted by Phys. Rev. E in July 25 2023):
		for (t_ in 1 : Nnodes) {
			if (t_ == 1){
				# Step 1: Begin with a seed graph containing a single node:
				sim.G 	= 	make_empty_graph(n = 0, directed = FALSE)
				j 		=	1
				sim.G 	= 	add_vertices(sim.G, 1) # V(sim.G) lists nodes/vertices in current network
			}
			if (t_ > 1){
				n_t		=	length(V(sim.G))	# n_t number of nodes in the network at the start of step t.
				# Step 2(a): Duplication:
				#	i.	 Select anchor node uniformly at random:
				anchorNode	=	sample(length(V(sim.G)), 1)# V(sim.G) lists nodes/vertices in current network
				# 	ii.	 Add a new node to the current network graph:
				j 		=	j + 1 # Updated number of nodes at time t_ (also the ID of the new node added)
				sim.G 	= 	add_vertices(sim.G, 1)# V(sim.G) lists nodes/vertices in current network
				# 	iii. Connect the new node to the anchor node's neighbors:
				anchorNodeNeighbors 	=	neighbors(sim.G, anchorNode, mode = "all")
				anchorNodeNeighbors 	=	as.numeric(anchorNodeNeighbors)
				newEdges	=	c(t(matrix(c(rep(j, length(anchorNodeNeighbors)), c(anchorNodeNeighbors)), ncol = 2)))
				sim.G 		= 	add_edges(sim.G, newEdges) # E(sim.G) lists edges in current network
				sim.G		=	simplify(sim.G) # Removes any multiple edges. See: E(sim.G)
				#
				# Step 2(b): Mutation:
				# 	i. Each edge connected to the new node is removed independently with probability q_del (or q_d).
				if (length(newEdges) > 0) {
					incidentNewNode 	=	incident(sim.G, j, mode = "all")
					incidentNewNode 	=	incidentNewNode[runif(length(incidentNewNode), 0, 1) <= q_d]
					sim.G 				=	delete_edges(sim.G, edges = as_ids(incidentNewNode))	}
				#
				# Step 2(c) Random Mutation
				# 	i. 	The step concludes by adding an edge between the new node and any
				# 		existing node at the start of step t with probability qnew/n(t) (or q_n / n(t)), 
				# 		where n(t) (or n_t) is the number of nodes in the network at the start of step t.
				existingNodes 	=	1 : (j - 1)
				newEdges			=	t(matrix(c(rep(j, j - 1), c(existingNodes)), ncol = 2))
				newEdges 			=	c(newEdges[, runif(j - 1, 0, 1) <= (q_n / n_t)])
				sim.G 				=	add_edges(sim.G, newEdges)# E(sim.G) lists the edges of the current network
				sim.G 				=	simplify(sim.G) 
			}
		}
	}



	if (model == "Watts-Strogatz"){# A model for an undirected network graph
		# Model's parameters = (K, beta), with Nnodes >> K >> log(Nnodes) >> 1 and 0 < beta < 1.
		# K is the mean degree. beta is the rewiring probability.
		# https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model
		# Given the desired number of nodes N, the mean degree K (assumed to be an even integer), 
		# and a parameter beta , all satisfying 0 ≤ beta ≤ 1  and N >> K >> ln(N) >> 1 
		# the model constructs an undirected graph with N nodes and (N*K)/2 edges in the following way:
		# Example inputs:   Nnodes = 100; K = 26; beta = 1/3;  with parameters = c(K, beta)
		K  		=	parameters[1] # Mean degree (should be an even integer).
		beta  	=	parameters[2] # Rewiring probablity.
		sim.G 	= 	sample_smallworld(dim = 1, size = Nnodes, nei = K / 2, p = beta, loops = FALSE, multiple = FALSE)
		#V(sim.G) # Lists vertices (nodes)
		#length(V(sim.G)) # Number of vertices (nodes). Equals Nnodes.
		#E(sim.G) # Lists edges
		#length(E(sim.G)) # Number of edges. Equals (Nnodes*K)/2
		#mean(degree(sim.G)) # Mean degree. Equals K
		#
		# Description from 'igraph' R package manual:
		# (see also https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model )
		# dim 			Integer constant, the dimension of the starting lattice.
		# size 		Integer constant, the size of the lattice along each dimension.[Equals number of nodes (Nnodes)]
		# nei 			Integer constant, the neighborhood within which the vertices
		# 				of the lattice will be connected 
		#				['nei' is denoted 'K'. With Nnodes >> K >> log(Nnodes) >> 1. See Wiki.]
		# p 			Real constant between zero and one, the rewiring probability.
		# 				[Often denoted 'beta', with 0 <= beta <= 1. See Wiki.].
		# loops 		Logical scalar, whether loops edges are allowed in the generated graph.
		# multiple	 	Logical scalar, whether multiple edges are allowed int the generated graph.
		#
		# Details
		#	First a lattice is created with the given dim, size and nei arguments. 
		#	Then the edges of the lattice are rewired uniformly randomly with probability p.
		#	Note that this function might create graphs with loops and/or multiple edges. 
		#	You can use simplify() to get rid of these.
		# Value:  A graph object.
	}

	if ((model == "Kretzschmar-Morris") | (model == "KM")){
		# This function simulates an undirected network graph from the Kretzschmar & Morris 
		# (KM) model. The following description of the KM model is taken from Goyal & Onnela 
		# (2020, pp.3, 6-7, in a paper entitled: "Framework For Converting Mechanistic 
		# Network Models To Probabilistic Models").
		# [The KM model is a] ... mechanistic model developed by Kretzschmar and Morris–hereafter 
		# referred to as the KM model (Kretzschmar & Morris, 1996, Mathematical Biosciences; 
       # Morris & Kretzschmar, 1997, AIDS), which played a significant role in identifying 
		# intervention priorities by highlighting the potential impact of concurrency on 
		# epidemic spread in sub-Saharan African  (Morris, Goodreau, & Moody, 2007, 
		# chapter in book Sexually Transmitted Diseases, pages 109–126).
		# The model continues to be the building block of more recent realistic models to study 
		# HIV (Palombi, et al. 2012). As it is believed that HIV epidemic in sub-Saharan Africa 
       # is driven by heterosexual relationships, the model only includes partnerships between 
		# people of the opposite sex, i.e., it is a model of a bipartite graph.
		# The network evolution under the KM model is based on individual-level stochastic rules
       # for partnership formation and dissolution. The population is fixed and the relationships
		# among the population form and dissolve over time. At each time t, an individual can form 
       # new partnerships, dissolve existing partnerships, or both. There are three key components
		# governing the formation and dissolution of relationships:
		# probability of pair formation (p_f), probability of pair separation (p_s), and 
		# a stochastic rule for partner mixing (phi), which can depend on the properties of the nodes.
		# The evolution of a network under the KM model is outlined below:
		# 1. Let g_t denote the network at time t.
		# 2. Repeat the following T_1 times:
		#  	(a) Simulate a Bernoulli process where X = 1 with probability p_f and X = 0 otherwise.
		# 	(b) If X = 1: 
		# 		(i) 	Draw two unconnected individuals at random, one male, i, and one female, j; 
		#		(ii) 	with probability phi(i,j) add edge (i,j) to g, 
		#				otherwise repeat (i) by redrawing two individuals at random.
		# 3. Every connected node pair splits up with probability p_s.
		# 
		# The resulting network following these steps, represents the network at time t + 1, 
		# denoted g_{t + 1}. To use the KM model to simulate an HIV epidemic, one must
		# specify an initial network at time 0, denote this network as g_0. 
		# Once g_0 is specified, the steps outlined above can be used to generate networks 
		# at subsequent times. In the KM model, the network g_0 is generated by starting
		# with an empty bipartite network with n_1 and n_2 nodes representing females and males, 
		# respectively, and then repeating the above steps a large number of times.
		# This procedure is commonly referred to as a burn-in step. After completing this 
		# large number of iterations, the resulting network, g_0, is used at time 0. 
		# The burn-in step ensures that the simulation of the HIV epidemic starts at the 
		# stationary state of the network generation process. 
		# In Section 5, we provide examples on how the MPMC framework can be used to derive 
		# a PMF for the stationary state of the process; therefore, one can sample a 
		# network g_0 from the stationary state instead of using the burn-in step described above. 
		# Note that ... we focus only on the generation of the networks and not on modeling
		# the HIV epidemic on the networks.
		# (Then on p. 7 of Goyal & Onnella 2020, who describe Step 1: Simulate from KM...)
		# We investigate a simple specification of the KM model, pure random mixing, 
		# and we use identical parameter values as the authors of the KM model when 
		# it was first proposed (Kretzschmar & Morris, 1996)
		# In the pure random mixing setting, there exists no preference for nodes to 
		# form edges based on the covariates of the nodes. 
		# The phi function for this setting is the following:
		# phi(i, j) = 1 if k_i < d_m and k_j < d_m; else, phi(i, j) = 0.    (8)
		# where k_i and k_j are the current degrees of nodes i and j.
		# The following parameters values were used in the original model:
		# n1 = n2 = 1000, p_f = 0.01, p_s = 0.005, T_1 = (n1 + n2)/2 - |g|, and d_m = 10.
		# |g| is the number of pairs (edges) at time t 
		# according to pp.176,181 of Kretzschmar & Morris (1996):
		# p.176: "Let X be the number of singles and P the number of pairs at time t; 
		# so N = X + 2P." and then on p.181: "(b) Repeat (a) [T_1 =] n = N/2 - P times."
		# (where N is the number of nodes according to p.167 of Kretzschmar & Morris, 1996:
		# 2. CONCURRENCY IN NETWORKS
		# By a sexual network we mean the total population together with all the
		# partnerships at a given moment. In graph theoretical terms, individuals
		# are referred to as vertices or nodes and a partnership as an edge between
		# two vertices. The network is then the set of all vertices together with
		# the set of all edges. So the sexual network may also contain isolated nodes,
		# that is, individuals who are single. Concurrent partnerships are described
		# by edges that emanate from the same vertex. ... we assume that we have a
		# population of size N with a given fixed partnership constellation.
		#
		# Example for testing:
		# Nnodes1 = 1000; 	# number of males in sexual network
		# Nnodes2 = 1000;  	# number of females in sexual network
		# p_f 	= 0.01; 	p_s = 0.005;  	d_m = 10; 		Time = 1825
		# p_f 	= 	probability of pair formation
		# p_s 	=	probability of pair separation
		# d_m 	=	maximum (current) degree used to define stochastic rule for partner mixing (phi)
		# Time	= 	simulation time (input).
		# Parameters p_f = 0.01; p_s = 0.005; Time = 1825 from Table 1 of Kretzschmar & Morris (1996)
		# d_m = 10 from p.7 of Goyal & Onnela (2020)
		
		# Input parameters of KM model:
		p_f 	=	parameters[1]# probability of pair formation  (p_f)
		p_s 	=	parameters[2]# probability of pair separation (p_s)
		d_m 	=	parameters[3]# maximum (current) degree used to define stochastic rule for partner mixing (phi)

		# Function to grow a network, g_t, at a time t (for t = 0,1,...)
		grow_g_t 	=	function(g_t, unconnectedMales, unconnectedFemales)	{ 
			# 1. Let g_t denote the network at time t.
			Nedges	=	gsize(g_t) # Number of edges in graph g_t
			T_1 	= 	(Nnodes1 + Nnodes2)/2 - Nedges
			for (s in 1 : T_1)	{
			# 2. Repeat the following T_1 times:
			#  	(a) Simulate a Bernoulli process where X = 1 with probability p_f and X = 0 otherwise.
				X		=	runif(1,0,1) <= p_f
			# 	(b) If X = 1: 
			# 		(i) 	Draw two unconnected individuals at random, one male, i, and one female, j; 
			#		(ii) 	with probability phi(i,j) add edge (i,j) to g, 
			#				otherwise repeat (i) by redrawing two individuals at random.
				if (X == 1){
					addEdge 		= 	FALSE
					while ((addEdge == FALSE) & (length(unconnectedMales) > 0) & (length(unconnectedFemales) > 0)){
						drawMale	= 	unconnectedMales[sample(length(unconnectedMales),1)		]# Drawing one male, i, 
						drawFemale	= 	unconnectedFemales[sample(length(unconnectedFemales),1)]# Drawing one female, j
						k_i			=	degree(g_t, v = drawMale, 	mode = "all")# k_i is current degree of node i
						k_j			=	degree(g_t, v = drawFemale,	mode = "all")# k_j is current degree of node j
						phi_ij 	=	ifelse((k_i < d_m) & (k_j < d_m), 1, 0)# stochastic rule for partner mixing (phi)
						addEdge	=	runif(1, 0, 1) <= phi_ij
						if (addEdge == TRUE) {g_t 	= 	add_edges(g_t, c(drawMale, drawFemale))}
						unconnectedMales				=	unconnectedMales		[unconnectedMales 	!= 	drawMale	]
						unconnectedFemales			=	unconnectedFemales	[unconnectedFemales 	!= 	drawFemale	]}}}
			# 3. Every connected node pair splits up with probability p_s.
				if (runif(1, 0, 1) <= p_s) {
					g_t 					= 	make_empty_graph(n = Nnodes1 + Nnodes2, directed = FALSE)
					unconnectedMales		=	1 : Nnodes1
					unconnectedFemales	=	(Nnodes1 + 1) : (Nnodes1 + Nnodes2)}
		return(list(g_t = g_t, unconnectedMales = unconnectedMales, unconnectedFemales = unconnectedFemales))}

		# Now generate an undirected network graph from the KM model, using the simulation strategy of 
		# Kretzschmar & Morris (1996, p.181) (also described in Goyal & Onnela, 2020, p.3), as follows:

		# Run burn-in step to generate a stationary network graph, g_0 (i.e., graph g_t at time t = 0):
		# Start with empty bipartite network with n_1 and n_2 nodes representing females and males (resp.). 
		g_0 					= 	make_empty_graph(n = Nnodes1 + Nnodes2, directed = FALSE) 
		unconnectedMales		=	1 : Nnodes1
		unconnectedFemales	=	(Nnodes1 + 1) : (Nnodes1 + Nnodes2)
		NburninIterations 	=	200
		for (s in 1 : NburninIterations)	{
			Out						= 	grow_g_t(g_0, unconnectedMales, unconnectedFemales) 
			g_0						=	Out$g_t 
			unconnectedMales		=	Out$unconnectedMales
			unconnectedFemales	=	Out$unconnectedFemales		}

		# Given a burned-in network graph g_0 (at time t = 0)
		# simulate an undirected network g_t at chosen (input) time t
		g_t 	=	g_0 
		for (s in 1 : (Time - 1))	{
			Out						= 	grow_g_t(g_t, unconnectedMales, unconnectedFemales) 
			g_t						=	Out$g_t 
			unconnectedMales		=	Out$unconnectedMales
			unconnectedFemales	=	Out$unconnectedFemales		}
		#
		Edgelist	=	as_edgelist(g_t)# # Below: convert Edgelist to a 'network' object:
		Types  	=	ifelse( (1 : (Nnodes1 + Nnodes2)) <= Nnodes1, 0, 1)
		sim.G 		= 	make_bipartite_graph(types = Types, edges = t(Edgelist), directed = FALSE)
	}


return(sim.G)}	# Output the simulated network graph.