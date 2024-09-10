# This code file defines list of candidate network 'terms' (statistics) from the 'ergm' R package, 
# and defines a list of candidate network 'terms' and 
# a function to compute all network statistics from the 'igraph' R package,
# to enable the calculation of toggled network statistics and the MPLE later on. 

#====================================================================================================================================
#  Define list of terms to compute various network statistics from the 'ergm' R package.
#====================================================================================================================================
#------------------------------------------------------------------------------------------------------------------------------------
#  For undirected networks:
#------------------------------------------------------------------------------------------------------------------------------------
oneTermsUndirected_ergm						= 	sort(unique(c(
	'altkstar(lambda = 2, fixed=TRUE)',		'balance',			'concurrent',		'concurrentties',
	'cycle(3)',		'cycle(4)',				'cyclicalties',		'degcrossprod',		
	'degrange(from = 2, to = +Inf)',		'degrange(from = 3, to = +Inf)',
	'degrange(from = 4, to = +Inf)',		'degrange(from = 5, to = +Inf)',
	'degree(2)',		'degree(3)',		'degree(4)',		'degree(5)',
	'dsp(2)',			'dsp(3)',			'dsp(4)',			'dsp(5)',		
	'degree1.5',		'edges',			'isolatededges',	'isolates',
	'kstar(2)',			'kstar(3)',			'kstar(4)',			'kstar(5)',
	'opentriad',		'threetrail',		'transitiveties',	'triangles',		'twopath')))

twoTermsUndirected_ergm						=	c('gwdegree(cutoff = Nnodes)') # term with 2 parameters (the second one being decay)
multipleTermsUndirected_ergm 				=	c('triadcensus')

oneTermsUndirectedValued_ergm				= 	sort(unique(c(
	'atleast(threshold = 1)', 				'atmost(threshold = 1)'		,
	'cyclicalweights(twopath=“min”, combine=“max”, affect=“min”)'		,
	'diff(attr, pow=1, dir=“t-h”, sign.action=“identity”)'				,
	'equalto(value=0, tolerance=0)',		'greaterthan(threshold=0)'	,
	'ininterval(lower=-Inf, upper=+Inf, open=c(TRUE,TRUE))'				,
	'smallerthan(threshold=0)',				'sum(pow=1)'				,
	'transitiveweights(twopath=“min”, combine=“max”, affect=“min”)')))

oneTermsUndirectedBipartite_ergm			=	sort(unique(c(
	oneTermsUndirected_ergm,
	'b1concurrent',
	'b1degrange(from = 2, to = +Inf)',		'b1degrange(from = 3, to = +Inf)',
	'b1degrange(from = 4, to = +Inf)',		'b1degrange(from = 5, to = +Inf)',
	'b1degree(2)',		'b1degree(3)',		'b1degree(4)',			'b1degree(5)',
	'b1dsp(2)',			'b1dsp(3)',			'b1dsp(4)',				'b1dsp(5)',
	'b1mindegree(2)',	'b1mindegree(3)',	'b1mindegree(4)',		'b1mindegree(5)',
	'b1star(2)',		'b1star(3)',		'b1star(4)',			'b1star(5)',
	'b2concurrent',		
	'b2degrange(from = 2, to = +Inf)',		'b2degrange(from = 3, to = +Inf)',
	'b2degrange(from = 4, to = +Inf)',		'b2degrange(from = 5, to = +Inf)',
	'b2degree(2)',		'b2degree(3)',		'b2degree(4)',			'b2degree(5)',
	'b2dsp(2)',			'b2dsp(3)',			'b2dsp(4)',				'b2dsp(5)',
	'b2mindegree(2)',	'b2mindegree(3)',	'b2mindegree(4)',		'b2mindegree(5)',
	'b2star(2)',		'b2star(3)',		'b2star(4)',			'b2star(5)')))

twoTermsUndirectedBipartite_ergm 			=	sort(unique(c(
	'gwb1degree(cutoff = Nnodes)', 		'gwb1dsp(cutoff = Nnodes)',
	'gwb2degree(cutoff = Nnodes)',		'gwb2dsp(cutoff = Nnodes)')))

multipleTermsUndirectedBipartite_ergm 		=	c('triadcensus')

multipleTermsUndirectedBipartiteDyadic_ergm =	sort(c('coincidence(levels=NULL,active=0)'))

oneTermsUndirectedNodalAttr_ergm			=	sort(unique(c(
	'absdiff(attr, pow = 1)', 				'absdiffcat(attr)',		'attrcov(attr, mat)',
	'mm(attrs, levels=NULL, levels2=-1)',	'nodecov(attr)',
	'nodemain(attr, form=“sum”)',			'nodefactor(attr, base=1, levels=-1)',
	'nodematch(attr, diff=FALSE, keep=NULL, levels=NULL)',
	'match(attr, diff=FALSE, keep=NULL, levels=NULL, form=“sum”)',
	'nodemix(attr, base=NULL, b1levels=NULL, b2levels=NULL, levels=NULL, levels2=-1)',
	'smalldiff(attr, cutoff)',	'sociality(attr=NULL, base=1, levels=NULL, nodes=-1)')))

oneTermsUndirectedDyadicCovariate_ergm		=	sort(unique(c(
	'dyadcov(x, attrname=NULL)', 	'hamming(x, cov, attrname=NULL)', 	'localtriangle(x)')))

oneTermsUndirectedBipartiteNodalAttr_ergm 	=	sort(unique(c(
	'b1cov(attr)', 		'b1factor(attr, base=1, levels=-1)',
	'b1nodematch(attr, diff=FALSE, keep=NULL, alpha=1, beta=1, byb2attr=NULL, levels=NULL)',
	'b1sociality',		'b2sociality',
	'b1starmix(k, attr, base=NULL, diff=TRUE)',
	'b1twostar(b1attr, b2attr, base=NULL, b1levels=NULL, b2levels=NULL, levels2=NULL)',
	'b2cov(attr)', 		'b2factor(attr, base=1, levels=-1)',
	'b2nodematch(attr, diff=FALSE, keep=NULL, alpha=1, beta=1, byb2attr=NULL, levels=NULL)',
	'b2starmix(k, attr, base=NULL, diff=TRUE)',
	'b2twostar(b1attr, b2attr, base=NULL, b1levels=NULL, b2levels=NULL, levels2=NULL)')))

TermsUndirected_ergm			=	c(	oneTermsUndirected_ergm, twoTermsUndirected_ergm, multipleTermsUndirected_ergm)
TermsUndirectedBipartite_ergm	=	c(	oneTermsUndirectedBipartite_ergm, twoTermsUndirectedBipartite_ergm, 
										multipleTermsUndirectedBipartite_ergm)


#------------------------------------------------------------------------------------------------------------------------------------
#  For directed networks:
#------------------------------------------------------------------------------------------------------------------------------------
oneTermsDirected_ergm						=	sort(unique(c(
	'asymmetric',	'balance',		'ctriple',		'ctriad',
	'cycle(2)', 	'cycle(3)',		'cycle(4)',		'cyclicalties',	
	'ddsp(1)',		'ddsp(2)',		'ddsp(3)',		'ddsp(4)',		'ddsp(5)',
	'desp(1)',		'desp(2)',		'desp(3)',		'desp(4)',		'desp(5)',
	'dnsp(1)',		'dnsp(2)',		'dnsp(3)',		'dnsp(4)',		'dnsp(5)',
	'edges',			
	'idegrange(from = 2, to = +Inf)','idegrange(from = 3, to = +Inf)',
	'idegrange(from = 4, to = +Inf)','idegrange(from = 5, to = +Inf)',
	'idegree(1)',	'idegree(2)',	'idegree(3)',	'idegree(4)',	'idegree(5)',
	'idegree1.5',	'intransitive',	'isolates',		
	'istar(1)',		'istar(2)',		'istar(3)',		'istar(4)',		'istar(5)',
	'm2star',		'mutual',		'nearsimmelian',
	'nsp(1)',		'nsp(2)',		'nsp(3)',		'nsp(4)',		'nsp(5)',
	'odegrange(from = 2, to = +Inf)',
	'odegree(1)',	'odegree(2)',	'odegree(3)',	'odegree(4)',	'odegree(5)',
	'odegree1.5',
	'ostar(1)',		'ostar(2)',		'ostar(3)',		'ostar(4)',		'ostar(5)',
	'simmelian',	'simmelianties',
	'threetrail(levels = 1)', 		'threetrail(levels = 2)',
	'threetrail(levels = 3)', 		'threetrail(levels = 4)',
	'transitive',	'transitiveties','triangles')))

twoTermsDirected_ergm						=	sort(unique(c(
	'dgwdsp(cutoff = Nnodes)',	'dgwesp(cutoff = Nnodes)',
	'dgwnsp(cutoff = Nnodes)',	'gwdsp(cutoff = Nnodes)',
	'gwesp(cutoff = Nnodes)',	'gwidegree(cutoff = Nnodes)',
	'gwnsp(cutoff = Nnodes)',	'gwodegree(cutoff = Nnodes)')))

multipleTermsDirected_ergm 	=	sort(c('triadcensus'))

oneTermsDirectedValued_ergm			=	sort(unique(c(
	'atleast(threshold = 1)', 			'atmost(threshold = 1)',
	'cyclicalweights(twopath="min", combine="max", affect="min")',
	'equalto(value=0, tolerance=0)',	'greaterthan(threshold=0)',
	'ininterval(lower=-Inf, upper=+Inf, open=c(TRUE,TRUE))',
	'nodecovar(center, transform)',		'nodeicovar(center, transform)',
	'nodeocovar(center, transform)',	'smallerthan(threshold=0)',		'sum(pow=1)',
	'transitiveweights(twopath="min", combine="max", affect="min")')))

oneTermsDirectedNodalAttr_ergm		= 	sort(unique(c(
	'absdiff(attr, pow = 1)',	'absdiffcat(attr)',	'attrcov(attr, mat)',
	'diff(attr, pow=1, dir=“t-h”, sign.action=“identity”)',
	'mm(attrs, levels=NULL, levels2=-1)',	'nodecov(attr)',
	'nodemain(attr, form=“sum”)',
	'nodefactor(attr, base=1, levels=-1)',
	'nodeicov(attr) nodeicov(attr, form=“sum”)',
	'nodeifactor(attr, base=1, levels=-1)',
	'nodematch(attr, diff=FALSE, keep=NULL, levels=NULL)',
	'match(attr, diff=FALSE, keep=NULL, levels=NULL, form=“sum”)',
	'nodemix(attr, base=NULL, b1levels=NULL, b2levels=NULL, levels=NULL, levels2=-1)',
	'nodeocov(attr)',	'nodeofactor(attr, base=1, levels=-1)',
	'receiver(base=1, nodes=-1)',	'sender(base=1, nodes=-1)',
	'smalldiff(attr, cutoff)')))

oneTermsDirectedDyadicCovariate_ergm=	sort(unique(c(
 	'dyadcov(x, attrname=NULL)', 'edgecov(x, attrname=NULL)', 'hamming(x, cov, attrname=NULL)')))

TermsDirected_ergm		<<-		c(oneTermsDirected_ergm, twoTermsDirected_ergm, multipleTermsDirected_ergm)

#------------------------------------------------------------------------------------------------------------------------------------
# Combine the terms for later use in Metropolis sampling updates.
#------------------------------------------------------------------------------------------------------------------------------------
twoTermsAll 	=	c(twoTermsUndirected_ergm, twoTermsDirected_ergm)

#------------------------------------------------------------------------------------------------------------------------------------
# Notes on ergm terms from ergm package:
#------------------------------------------------------------------------------------------------------------------------------------
# For directed networks:
#	'density' equals 'edges' or 'istar(1)' or 'ostar(1)' divided by Nnodes*(Nnodes-1) . 
# 	'edges' is equal to 'ostar(1)', 'istar(1)', and 'nonzero'.
#	'desp' is the same as 'esp'
#	'm2star' is the same as 'twopath'
#	'triangles' equals 'ttriple' + 'ctriple', so at most two of the three terms can be in a model.
#	'nodecov' equals 'nodeicov' plus 'nodeocov'.
# For undirected networks:
# 	'density' equals 'kstar(1)' or 'edges' divided by n(n-1)/2;
#	'edges' equals 'kstar(1)' and 'nonzero'.
#	The 'edgecov' and 'dyadcov' terms are equivalent for undirected networks.
#  'tripercent' can take a while to compute even for a small network (small number of nodes)
#	found following error when running ergm() with variable 'degcor':
# 		Error in if (any(low.drop.theta)) message(paste("Observed statistic(s)",  : 
#		missing value where TRUE/FALSE needed
# For directed or undirected networks:
#	'edges' equals 'nonzero'
# 	'meandeg' (mean degree) is the average number of edges per node 
#  		in the network graph: Average Degree = Total Edges/Total Nodes
#		This term is a constant multiple of both edges and density . 
# 	'triangles' equals 'triangle';
#	'ttriple' equals 'ttriad'
#	'triangle' equals 'ttriple' + 'ctriple' for a directed network, 
#		so at most two of the three terms can be in a model.
#       Therefore, only including 'triangle' ('triangles') and 'ctriple'
# 		for a directed network
#	'cycle(5)' can take a long time to compute even for a network with a non-large number of nodes.
#	In 'localtriangle(x)', the argument x is an undirected network or an symmetric adjacency matrix 
#		that specifies whether the two nodes are in the same neighborhood. 
#	'triadcensus(levels)' is also implemented by 'igraph' R package (and also for directed networks).
# Therefore, the terms 'density', 'istar(1)', 'ostar(1)', 'nonzero', 'meandegree',
#	'twopath', 'triangle', are redundant and thus not considered as model terms above.
# Also excluded are terms which required to specify specific numerical values 
#	(e.g., specific degree values) and thus were difficult to specify/run automatically in a loop.
# altkstar(lambda, fixed=FALSE) did not converge for some reason (thus excluded)
# 	"altkstar(lambda, fixed = FALSE) — alternating k-stars: This term adds one network statistic
# 	to the model equal to a weighted alternating sequence of k-star statistics with weight parameter lambda. 
# 	This is the version given in Snijders et al. (2006). The gwdegree and altkstar terms produce 
# 	mathematically equivalent models, as long as they are used together with the edges [or kstar(1)] term.
# 	See Section 3 and especially equation (13) of Hunter (2007) for details. 
# 	The optional argument fixed indicates whether the scale parameter lambda is to be fit as a 
# 	curved exponential-family model (see Hunter and Handcock 2006). The default is FALSE, which means 
# 	the scale parameter is not fixed and thus the model is a CEF model. This term can only be used 
# 	with undirected networks." (p.8 of Morris et al. (2008))
# 	Robins et al. (2007, p.198), when considering fixed lambda, considered lambda = 2.
# 	Morris, M., Handcock, M. S., & Hunter, D. R. (2008). Specification of Exponential-Family
#		Random Graph Models: Terms and Computational Aspects. 
#		Journal of statistical software, 24(4), 1548. doi: 10.18637/jss.v024.i04
# 	Robins, G., Snijders, T., Wang, P., Handcock, M., & Pattison, P. (2007). Recent developments 
#		in exponential random graph (p*) models for social networks. Social networks, 29(2), 192-215.
# For a model bipartite network, 
#	-	found following error when running ergm() with variable 'degcor':
# 		Error in if (any(low.drop.theta)) message(paste("Observed statistic(s)",  : 
#		missing value where TRUE/FALSE needed
#	-	found following error when running ergm() with variable ‘gwdegree’:
#		Error in `ergm_Init_abort()`:
#		! In term ‘gwdegree’ in package ‘ergm’: Term may not be used with 
#		networks with bipartite > 0.
# ergm package indicated that 'threepath' is inaccurately named and is equivalent to 'threetrail'.
# Therefore, in the above model terms (for undirected networks, we are not including 'threepath' 
# and just including 'threetrail'.
# 'gwdegree' term is not applicable to bipartite networks.
# In term ‘cycle’ in package ‘ergm’: cycles of length less than 3 cannot exist in an undirected network
# and their statistics (i.e., 'cycle(2)') will be omitted from the above model terms.
# No ERGM term can be used for bipartite directed networks.
#if ((isDirectedGraph == TRUE)&(isBipartiteGraph == TRUE)){stop("No ERGM term can be used for bipartite directed networks.")}
# For more details about and other available ergm package model terms, see:
# https://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html#idegree-ergmTerm-036b6e10 
# Define network statistics from igraph:


#====================================================================================================================================
#  Define function to compute various network statistics from the 'igraph' R package.
#====================================================================================================================================
Terms_igraph   	=  unique(c('mean_node_deg_all',				'var_node_deg_all',
							'mean_node_deg_out',				'var_node_deg_out',
							'mean_node_deg_in',					'var_node_deg_in',
							'assortativity_degree',
							'mean_node_betweenness',			'var_node_betweenness',				'n_clusters_edge_betweenness',
							'n_cliques',	'largest_cliques',	'count_max_cliques',				'clique_num',
							'mean_clique_size_counts',			'var_clique_size_counts',
							'mean_node_closeness_all',			'var_node_closeness_all',
							'mean_node_closeness_out',			'var_node_closeness_out',
							'mean_node_closeness_in',			'var_node_closeness_in',
							'mean_node_cocitation',				'var_node_cocitation',
							'mean_node_cohesion',				'var_node_cohesion',
							'mean_node_max_cohesion',			'var_node_max_cohesion',			'n_cohesive_blocks',
							'count_components',					'count_components_weak',			'count_components_strong',
							'mean_ComponentSize',				'var_ComponentSize',
							'mean_constraint',					'var_constraint',
							'count_motifs_of_Size3',
							'mean_node_coreness_all',			'var_node_coreness_all',
							'mean_node_coreness_out',			'var_node_coreness_out',
							'mean_node_coreness_in',			'var_node_coreness_in',
							'automorphism_group_size',			'diameter',
							'mean_node_diversity',				'var_node_diversity',
							'dyad_census_mut',					'dyad_census_asym',					'dyad_census_null',
							'mean_node_eccentricity_all',		'var_node_eccentricity_all',
							'mean_node_eccentricity_out',		'var_node_eccentricity_out',
							'mean_node_eccentricity_in',		'var_node_eccentricity_in',
							'edge_connectivity',				'edge_density',
							'mean_node_centrality',				'var_node_centrality',				'eigen_node_centrality',
							'girth',							'global_efficiency',
							'mean_node_local_efficiency_all',	'var_node_local_efficiency_all',
							'mean_node_local_efficiency_out',	'var_node_local_efficiency_out',
							'mean_node_local_efficiency_in',	'var_node_local_efficiency_in',
							'mean_node_harmonic_centrality_all','var_node_harmonic_centrality_all',
							'mean_node_harmonic_centrality_out','var_node_harmonic_centrality_out',
							'mean_node_harmonic_centrality_in',	'var_node_harmonic_centrality_in',
							'mean_node_hub_score',				'var_node_hub_score',
							'mean_node_authority_score',		'var_node_authority_score',
							'mean_node_KNNmeandegree_all',		'var_node_KNNmeandegree_all',
							'mean_node_KNNmeandegree_out',		'var_node_KNNmeandegree_out',
							'mean_node_KNNmeandegree_in',		'var_node_KNNmeandegree_in',
							'radius_all',						'radius_out',						'radius_in',
							'reciprocity',						'clustering_coef',
							'triad_census_003',		'triad_census_012',		'triad_census_102',		'triad_census_021D',
							'triad_census_021U',	'triad_census_021C',	'triad_census_111D',	'triad_census_111U',
							'triad_census_030T',	'triad_census_030C',	'triad_census_201',		'triad_census_120D',
							'triad_census_120U',	'triad_census_120C',	'triad_census_210',		'triad_census_300',
							'mean_node_triangles',	'var_node_triangles')	)

# This function computes from the input network graph g network statistics (inputTerms).
# The default 'inputTerms = NULL' means to compute all the network statistics, but this can take a long time.
networkStats_igraph 	=	function(g, inputTerms = NULL){
	#---------------------------------------------------------------------------------------------------------------------------------
	# Basic characteristics of network graph g:
	#---------------------------------------------------------------------------------------------------------------------------------
	Nnodes 			=	vcount(g) # equals gorder(g)
	isDirectedGraph	=	is_directed(g)
	isBipartiteGraph	=	is_bipartite(g)
	isSimpleGraph		=	is_simple(g) # Simple graphs are graphs which do not contain loop and multiple edges.
	Weights 			= 	rep(1, ecount(g)) # Assign weight 1 to each edge in the graph g. Any weighted graph is already dichotomized.
	UndirectedNonbipartite	=	(isDirectedGraph == FALSE) & (isBipartiteGraph == FALSE)
	UndirectedBipartite		=	(isDirectedGraph == FALSE) & (isBipartiteGraph == TRUE )
	DidInputTerms			=	!is.null(inputTerms)
	if (DidInputTerms == FALSE)	{inputTerms 	=  Terms_igraph}
	networkStats 		=	matrix(, nrow = 1, ncol = 0) # Initialize.
	#---------------------------------------------------------------------------------------------------------------------------------
	if (('mean_node_deg_all' %in% inputTerms) | ('var_node_deg_all' %in% inputTerms)) {
		out 			= 	degree(g, mode = "all", normalized = FALSE)}# degree for each node in the network graph g.
	if ('mean_node_deg_all' %in% inputTerms) {
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_deg_all'	}
	if ('var_node_deg_all' %in% inputTerms) 	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm=TRUE))^2, na.rm=TRUE)/sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_deg_all'	}
	# For a directed or undirected graph g, mode = 'all' and mode = 'total' give identical results.
	if (('mean_node_deg_out' %in% inputTerms) | ('var_node_deg_out' %in% inputTerms)) {
		out 		= 	degree(g, mode = "out", normalized = FALSE)}# outdegree for each node in the network.
	if ('mean_node_deg_out' %in% inputTerms) {
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_deg_out'	}
	if ('var_node_deg_out' %in% inputTerms) {
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm=TRUE))^2, na.rm=TRUE)/sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_deg_out'	}
	if (('mean_node_deg_in' %in% inputTerms) | ('var_node_deg_in' %in% inputTerms)) {
		out 		= 	degree(g, mode = "in", normalized = FALSE)}# outdegree for each node in the network.
	if ('mean_node_deg_in' %in% inputTerms) {
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))	
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_deg_in'	}
	if ('var_node_deg_in' %in% inputTerms) {
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm=TRUE))^2, na.rm=TRUE)/sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_deg_in'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# alpha_centrality() calculates the alpha centrality of some (or all) vertices in a graph g.
	# Details:  The alpha centrality measure can be considered as a generalization of eigenvector centrality to
	# directed graphs. It was proposed by Bonacich in 2001 (see reference below).
	# The alpha centrality of the vertices in a graph g is defined as the solution of the following matrix equation:
	# x = alpha * t(A)*x + e; where A is the (not necessarily symmetric) adjacency matrix of the graph g, 
	# e is the vector of exogenous sources of status of the vertices and alpha is the relative importance 
	# of the endogenous versus exogenous factors.
	# Value: A numeric vector contaning the centrality scores for the selected vertices.
	# Warning: Singular adjacency matrices cause problems for this algorithm, the routine may fail is certain cases.
	# 	out <- alpha_centrality(g, alpha = 1, loops = FALSE, exo = 1, sparse = TRUE)
	# 	networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
	# 	#colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_alpha1_centrality'
	# 	networkStats	=	cbind(networkStats, sum((out- mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
	# 	#colnames(networkStats)[ncol(networkStats)]	<-	'var_node_alpha1_centrality'
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('assortativity_degree' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, assortativity_degree(g, directed = TRUE)   )
		colnames(networkStats)[ncol(networkStats)]	<-	'assortativity_degree'				}
	#assortativity(g, types1, types2 = NULL, directed = isDirectedGraph)
	#assortativity_nominal(g, types, directed = isDirectedGraph)
	#---------------------------------------------------------------------------------------------------------------------------------
	if (('mean_node_betweenness' %in% inputTerms) | ('var_node_betweenness' %in% inputTerms))	{
		out				=	betweenness(g, directed = isDirectedGraph, normalized = FALSE)	}
	if ('mean_node_betweenness' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_betweenness'		}
	if ('var_node_betweenness' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_betweenness'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('n_clusters_edge_betweenness' %in% inputTerms){
		networkStats	=	cbind(	networkStats, 
									length(cluster_edge_betweenness(g, directed = isDirectedGraph, 
									edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,modularity = TRUE, membership = TRUE)))
		colnames(networkStats)[ncol(networkStats)]	<-	'n_clusters_edge_betweenness'		}
	# centr_betw_tmax(g, directed = isDirectedGraph)
   	# above gives theoretical maximum (unnormalized) graph (g) betweenness centrality score for graphs
	# with given order and other parameters.
	# centr_clo_tmax(g, mode = c("out", "in", "all", "total"))
	# centr_degree_tmax(g, mode = c("all", "out", "in", "total"),loops = FALSE)
	# centr_eigen_tmax(g, directed = isDirectedGraph, scale = TRUE)
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('n_cliques' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, length(cliques(g)))
		colnames(networkStats)[ncol(networkStats)]	<-	'n_cliques'				}
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('largest_cliques' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, length(largest_cliques(g)))
		colnames(networkStats)[ncol(networkStats)]	<-	'largest_cliques'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('count_max_cliques' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, count_max_cliques(g))
		colnames(networkStats)[ncol(networkStats)]	<-	'count_max_cliques'		}
	# max_cliques(g,min=NULL,max=NULL,subset=NULL,file=NULL)# equals count_max_cliques(g)
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('clique_num' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, clique_num(g))# The size of the largest clique(s).
		colnames(networkStats)[ncol(networkStats)]	<-	'clique_num'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	if (('mean_clique_size_counts' %in% inputTerms) | ('var_clique_size_counts' %in% inputTerms))	{
		out = clique_size_counts(g) }# a histogram of clique sizes, between the given minimum and maximum clique size.
	if ('mean_clique_size_counts' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_clique_size_counts'		}
	if ('var_clique_size_counts' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_clique_size_counts'		}
	#---------------------------------------------------------------------------------------------------------------------------------
   	# Closeness centrality measures how many steps is required to access every other vertex from a given vertex.
	if (('mean_node_closeness_all' %in% inputTerms) | ('var_node_closeness_all' %in% inputTerms))	{
		quiet(out <- closeness(g, mode = 'total', normalized = FALSE)) }
	if ('mean_node_closeness_all' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_closeness_all'	}
	if ('var_node_closeness_all' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_closeness_all'	}
	# For a directed or undirected graph (g), mode = 'all' and mode = 'total' give identical results.
	if (('mean_node_closeness_out' %in% inputTerms) | ('var_node_closeness_out' %in% inputTerms))	{
		out		= 	closeness(g, mode = "out", normalized = FALSE)					}
	if ('mean_node_closeness_out' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_closeness_out'		}
	if ('var_node_closeness_out' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_closeness_out'		}
	if (('mean_node_closeness_in' %in% inputTerms) | ('var_node_closeness_in' %in% inputTerms))	{
		out		= 	closeness(g, mode = "in", normalized = FALSE)						}
	if ('mean_node_closeness_in' %in% inputTerms) 	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))				
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_closeness_in'		}
	if ('var_node_closeness_in' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_closeness_in'			}
	#---------------------------------------------------------------------------------------------------------------------------------
	# cocitation(g) simply counts how many types two vertices are cocited.
	# Two vertices are cocited if there is another vertex citing both of them. 
	# cocitation(g) Nnodes-by-Nnodes matrix output. isSymmetric(cocitation(g)) equals TRUE
	if (('mean_node_cocitation' %in% inputTerms) | ('var_node_cocitation' %in% inputTerms))	{
		out		=	rowSums(cocitation(g))	}
	if ('mean_node_cocitation' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_cocitation'		}
	if ('var_node_cocitation' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_cocitation'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	if 	(('mean_node_cohesion' %in% inputTerms) | ('var_node_cohesion' %in% inputTerms)|('mean_node_max_cohesion' %in% inputTerms)|
		('var_node_max_cohesion' %in% inputTerms) |	('n_cohesive_blocks' %in% inputTerms))	{# The input graph (g) 'must be undirected and simple. (See is_simple().)'
		# cohesive_blocks(g) # For this the graph object g 'must be undirected and simple.'(from 'igraph' package manual)
		# cohesion(cohesive_blocks(g)) returns the cohesion of each block.
		# max_cohesion(cohesive_blocks(g)) # gives, for each vertex, the cohesion of its most cohesive block.
		# length(cohesive_blocks(g)) returns a numeric scalar, the number of blocks.
		# The input graph (g) 'must be undirected and simple. (See is_simple().)' From 'igraph' package manual.
		if (isSimpleGraph == TRUE)	{
			out		=	cohesion(cohesive_blocks(g))
			if ('mean_node_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
				colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_cohesion'	}
			if ('var_node_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
				colnames(networkStats)[ncol(networkStats)]	<-	'var_node_cohesion'		}
			if 	(('mean_node_max_cohesion' %in% inputTerms)|('var_node_max_cohesion' %in% inputTerms) |	('n_cohesive_blocks' %in% inputTerms))	{
				out		=	max_cohesion(cohesive_blocks(g))	}
			if 	('mean_node_max_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
				colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_max_cohesion'	}
			if 	('var_node_max_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
				colnames(networkStats)[ncol(networkStats)]	<-	'var_node_max_cohesion'		}
			if 	('n_cohesive_blocks' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, length(cohesive_blocks(g)))
				colnames(networkStats)[ncol(networkStats)]	<-	'n_cohesive_blocks'			}		}
		if (isSimpleGraph == FALSE)	{
			if ('mean_node_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_cohesion'		}
			if ('var_node_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'var_node_cohesion'			}
			if ('mean_node_max_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_max_cohesion'	}
			if ('var_node_max_cohesion' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'var_node_max_cohesion'		}
			if ('n_cohesive_blocks' %in% inputTerms)	{
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'n_cohesive_blocks'			}	}	}
	#---------------------------------------------------------------------------------------------------------------------------------
	#is_connected(g, mode = c("weak", "strong"))
	#	decides whether the graph g is weakly or strongly connected. The null graph g is considered disconnected.
	#components(g, mode = c("weak", "strong"))
	# 	finds the maximal (weakly or strongly) connected components of a graph g.
	#count_components(g, mode = c("weak", "strong"))
	# 	does almost the same as components() but returns only the number of clusters found 
	#	instead of returning the actual clusters.
	#component_distribution() creates a histogram for the maximal connected component sizes.
	#largest_component(, mode = c("weak", "strong"))
	#	returns the largest connected component of a graph g. 
	#	For directed graphs,	optionally the largest weakly or strongly connected component. 
	#	In case of a tie, the first component by vertex ID order is returned. 
	#	Vertex IDs from the original graph g are not retained in the returned graph.
	#The weakly connected components are found by a simple breadth-first search. The strongly connected
	#components are implemented by two consecutive depth-first searches.
	#mode 	Character string, either "weak" or "strong”. 
	#		For directed graphs "weak" implies weakly, "strong" strongly connected components to search. 
	#		It is ignored for undirected graphs.
	if ('count_components' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, count_components(g))
		colnames(networkStats)[ncol(networkStats)]	<-	'count_components'	}
	if ('count_components_weak' %in% inputTerms) 	{
		networkStats	=	cbind(networkStats, count_components(g, mode = "weak"))
		colnames(networkStats)[ncol(networkStats)]	<-	'count_components_weak'		}
	if ('count_components_strong' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, count_components(g, mode = "strong"))
		colnames(networkStats)[ncol(networkStats)]	<-	'count_components_strong'	}
	if (('mean_ComponentSize' %in% inputTerms) | ('var_ComponentSize' %in% inputTerms))	{
		out		=	component_distribution(g, cumulative = FALSE)
		probs	=	out / sum(out)
		vals	=	0 : (length(out) - 1)	}
	if ('mean_ComponentSize' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum(vals * probs))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_ComponentSize'	}
	if ('var_ComponentSize' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_ComponentSize'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	if (('mean_constraint' %in% inputTerms) | ('var_constraint' %in% inputTerms))		{
		out 	=	constraint(g) }# Burt’s constraint for each vertex.
	if ('mean_constraint' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out))				
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_constraint'	}
	if ('var_constraint' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_constraint'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	if ('count_motifs_of_Size3' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, count_motifs(g, size = 3))	# size=3, so rep(0, 3) = rep(0, size)
		colnames(networkStats)[ncol(networkStats)]	<-	'count_motifs_of_Size3'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The k-core of graph g is a maximal subgraph in which each vertex has at least degree k. 
	# The coreness of a vertex is k if it belongs to the k-core but not to the (k+1)-core.
	# coreness(g, mode = c("all", "out", "in"))
	if (('mean_node_coreness_all' %in% inputTerms) | ('var_node_coreness_all' %in% inputTerms)){
		out 	=	coreness(g, mode = "all")	}
	if ('mean_node_coreness_all' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_coreness_all'	}
	if ('var_node_coreness_all' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_coreness_all'		}
	# For a directed or undirected graph g, mode = "all" and mode = "total" give identical results.
	if (('mean_node_coreness_out' %in% inputTerms) | ('var_node_coreness_out' %in% inputTerms)){
		out		= 	coreness(g, mode = "out")	}
	if ('mean_node_coreness_out' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_coreness_out'	}
	if ('var_node_coreness_out' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_coreness_out'		}
	if (('mean_node_coreness_in' %in% inputTerms) | ('var_node_coreness_in' %in% inputTerms))	{
		out		= 	coreness(g, mode = "in")	}
	if ('mean_node_coreness_in' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_coreness_in'			}
	if ('var_node_coreness_in' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_coreness_in'			}
	#---------------------------------------------------------------------------------------------------------------------------------
	# Calculate the number of automorphisms of a graph g, i.e. the number of isomorphisms to itself.
	# Below gives the size of the automorphism group of the input graph g.
	if ('automorphism_group_size' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, color = NULL, as.numeric(count_automorphisms(g)$group_size)	)
		colnames(networkStats)[ncol(networkStats)]	<-	'automorphism_group_size'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The diameter of a graph g is the length of the longest geodesic.
	if ('diameter' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, diameter(g, directed = isDirectedGraph)	)
		colnames(networkStats)[ncol(networkStats)]	<-	'diameter'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	if (('mean_node_diversity' %in% inputTerms) | ('var_node_diversity' %in% inputTerms))	{
		if (isDirectedGraph == FALSE){
			out		=	diversity(g, weights = Weights)
			if ('mean_node_diversity' %in% inputTerms){
				networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
				colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_diversity'	}
			if ('var_node_diversity' %in% inputTerms){
				networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
				colnames(networkStats)[ncol(networkStats)]	<-	'var_node_diversity'	}	}
		if (isDirectedGraph == TRUE){
			if ('mean_node_diversity' %in% inputTerms){
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_diversity'	}
			if ('var_node_diversity' %in% inputTerms){
				networkStats	=	cbind(networkStats, NA)
				colnames(networkStats)[ncol(networkStats)]	<-	'var_node_diversity'	}	}	}
	# Error in diversity(sample_pa(n = 10)) : 
	# At core/properties/basic_properties.c:140 : Diversity measure works with undirected graphs only. Invalid value
	#---------------------------------------------------------------------------------------------------------------------------------
	# Classify dyads in a directed graphs. The relationship between each pair of vertices is measured. 
	# It can be in three states: mutual, asymmetric or non-existent.
	if (('dyad_census_mut' %in% inputTerms) | ('dyad_census_asym' %in% inputTerms) | ('dyad_census_null' %in% inputTerms))	{
		quiet(out		<-	dyad_census(g)) }# Using quiet() because outputs a warning for an undirected graph g
	if ('dyad_census_mut' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, out$mut)
		colnames(networkStats)[ncol(networkStats)]	<-	'dyad_census_mut'	}
	if ('dyad_census_asym' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, out$asym)
		colnames(networkStats)[ncol(networkStats)]	<-	'dyad_census_asym'	}
	if ('dyad_census_null' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, out$null)
		colnames(networkStats)[ncol(networkStats)]	<-	'dyad_census_null'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The eccentricity of a vertex is its shortest path distance from the farthest other node in the graph g.
	# The inverse distance to an unreachable vertex is considered to be zero.
	# eccentricity(g, vids = V(g), mode = c("all", "out", "in", "total"))
	# mode 	Character string, defining the types of the paths used for measuring the distance in directed graphs. 
	# 			"out" follows paths along the edge directions only,
	# 			"in" traverses the edges in reverse, while "all" ignores edge directions. 
	# 			This argument is ignored for undirected graphs.
	if (('mean_node_eccentricity_all' %in% inputTerms) | ('var_node_eccentricity_all' %in% inputTerms)){
		out 	=	eccentricity(g, mode = "all")	}
	if ('mean_node_eccentricity_all' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_eccentricity_all'	}
	if ('var_node_eccentricity_all' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_eccentricity_all'		}
	# For a directed or undirected graph g, mode = "all" and mode = "total" give identical results.
	if (('mean_node_eccentricity_out' %in% inputTerms) | ('var_node_eccentricity_out' %in% inputTerms)){
		out		= 	eccentricity(g, mode = "out")	}
	if ('mean_node_eccentricity_out' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_eccentricity_out'	}
	if ('var_node_eccentricity_out' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_eccentricity_out'		}
	if (('mean_node_eccentricity_in' %in% inputTerms) | ('var_node_eccentricity_in' %in% inputTerms)){
		out		= 	eccentricity(g, mode = "in")	}
	if ('mean_node_eccentricity_in' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_eccentricity_in'		}
	if ('var_node_eccentricity_in' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_eccentricity_in'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The edge connectivity of a graph (g) or two vertices, this is recently also called group adhesion.
	if ('edge_connectivity' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, edge_connectivity(g))
		colnames(networkStats)[ncol(networkStats)]	<-	'edge_connectivity'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The density of a graph (g) is the ratio of the actual number of edges and the largest possible number of
	# edges in the graph, assuming that no multi-edges are present.
	if ('edge_density' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, edge_density(g, loops = FALSE))
		colnames(networkStats)[ncol(networkStats)]	<-	'edge_density'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	if (('mean_node_centrality' %in% inputTerms) | ('var_node_centrality' %in% inputTerms) | ('eigen_node_centrality' %in% inputTerms)){
		quiet(out <- eigen_centrality(g, directed = isDirectedGraph, scale = TRUE, options = arpack_defaults))	}
	# Value:		A named list with components:
	# vector 	A vector containing the centrality scores.
	# value 	The eigenvalue corresponding to the calculated eigenvector, i.e. the centrality scores.
	# options 	A named list, information about the underlying ARPACK computation. See arpack() for the details.
	if ('mean_node_centrality' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, mean(out$vector, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_centrality'	}
	if ('var_node_centrality' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, sum((out$vector - mean(out$vector, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_centrality'	}
	if ('eigen_node_centrality' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, out$value)
		colnames(networkStats)[ncol(networkStats)]	<-	'eigen_node_centrality'	}
	#-------------------------------------------------------------------------------------------------------------
	# The girth of a graph (g) is the length of the shortest circle in it.
	if ('girth' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, girth(g, circle = FALSE)$girth)
		colnames(networkStats)[ncol(networkStats)]	<-	'girth'					}
	#-------------------------------------------------------------------------------------------------------------
	# Calculate the global, average, or local efficiency of a network.
	# For global_efficiency(), the global efficiency of the graph (g) as a single number.
	# For average_local_efficiency(), the average local efficiency of the graph as a single number.
	# For local_efficiency(), the local efficiency of each vertex in a vector.
	if ('global_efficiency' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, global_efficiency(g, directed = isDirectedGraph))
		colnames(networkStats)[ncol(networkStats)]	<-	'global_efficiency'	}
	if ('mean_node_local_efficiency_all' %in% inputTerms){
		networkStats	=	cbind(networkStats, average_local_efficiency(g, directed = isDirectedGraph, mode = "all"))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_local_efficiency_all'}
	if ('var_node_local_efficiency_all' %in% inputTerms)	{
		out 	= 	local_efficiency(g, directed = isDirectedGraph, mode = "all")
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_local_efficiency_all'	}
	# For a directed or undirected graph g, mode = "all" and mode = "total" give identical results.
	if ('mean_node_local_efficiency_out' %in% inputTerms){
		networkStats	=	cbind(networkStats, average_local_efficiency(g, directed = isDirectedGraph, mode = "out"))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_local_efficiency_out'}
	if ('var_node_local_efficiency_out' %in% inputTerms)	{
		out 	=	local_efficiency(g, directed = isDirectedGraph, mode = "out")
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_local_efficiency_out'	}
	if ('mean_node_local_efficiency_in' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, average_local_efficiency(g, directed = isDirectedGraph, mode = "in"))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_local_efficiency_in'	}
	if ('var_node_local_efficiency_in' %in% inputTerms)	{
		out 	=	local_efficiency(g, directed = isDirectedGraph, mode = "in")
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_local_efficiency_in'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The harmonic centrality of a vertex is the mean inverse distance to all other vertices.
	# The inverse distance to an unreachable vertex is considered to be zero.
	# harmonic_centrality(g, vids = V(g), mode = c("out", "in", "all", "total"), normalized = FALSE)
	# mode 	Character string, defining the types of the paths used for measuring the distance in directed graphs. 
	# 			"out" follows paths along the edge directions only,
	# 			"in" traverses the edges in reverse, while "all" ignores edge directions. 
	# 			This argument is ignored for undirected graphs.
	if (('mean_node_harmonic_centrality_all' %in% inputTerms) | ('var_node_harmonic_centrality_all' %in% inputTerms)){
		out 	=	harmonic_centrality(g, normalized = FALSE)	}
	if ('mean_node_harmonic_centrality_all' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_harmonic_centrality_all'	}
	if ('var_node_harmonic_centrality_all' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_harmonic_centrality_all'	}
	# For a directed or undirected graph (g), mode = "all" and mode = "total" give identical results.
	if (('mean_node_harmonic_centrality_out' %in% inputTerms) | ('var_node_harmonic_centrality_out' %in% inputTerms)){
		out		= 	harmonic_centrality(g, mode = "out", normalized = FALSE)	}
	if ('mean_node_harmonic_centrality_out' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_harmonic_centrality_out'		}
	if ('var_node_harmonic_centrality_out' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm=TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_harmonic_centrality_out'		}
	if (('mean_node_harmonic_centrality_in' %in% inputTerms) | ('var_node_harmonic_centrality_in' %in% inputTerms)){
		out		= 	harmonic_centrality(g, mode = "in", normalized = FALSE)		}
	if ('mean_node_harmonic_centrality_in' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_harmonic_centrality_in'		}
	if ('var_node_harmonic_centrality_in' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_harmonic_centrality_in'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The hub scores of the vertices are defined as the principal eigenvector of A*t(A) (matrix product), 
	# where A is the adjacency matrix of the graph (g).
	# Similarly, the authority scores of the vertices are defined as the principal eigenvector of t(A)*A, 
	# where A is the adjacency matrix of the graph (g).
	# For undirected matrices the adjacency matrix is symmetric and the hub scores equal authority scores.
	# "hub_score()" not computed for bipartite graphs (because ARPACK solver can be unstable and result in error messages).
	if (!isBipartiteGraph)	{
		if (('mean_node_hub_score' %in% inputTerms) | ('var_node_hub_score' %in% inputTerms)){
			out				= 	hub_score(g, scale = TRUE)$vector					}
		if ('mean_node_hub_score' %in% inputTerms)	{
			networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
			colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_hub_score'	}
		if ('var_node_hub_score' %in% inputTerms)	{
			networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))	
			colnames(networkStats)[ncol(networkStats)]	<-	'var_node_hub_score'	}
		if (('mean_node_authority_score' %in% inputTerms) | ('var_node_authority_score' %in% inputTerms)){
			out				= 	authority_score(g, scale = TRUE)$vector	}
		if ('mean_node_authority_score' %in% inputTerms){
			networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
			colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_authority_score'	}
		if ('var_node_authority_score' %in% inputTerms){
			networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))
			colnames(networkStats)[ncol(networkStats)]	<-	'var_node_authority_score'	}	}
	if (isBipartiteGraph)	{
		if ('mean_node_hub_score' %in% inputTerms)			{networkStats	=	cbind(networkStats, NA)
			colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_hub_score'	}
		if ('var_node_hub_score' %in% inputTerms)			{networkStats	=	cbind(networkStats, NA)	
			colnames(networkStats)[ncol(networkStats)]	<-	'var_node_hub_score'	}
		if ('mean_node_authority_score' %in% inputTerms)		{networkStats	=	cbind(networkStats, NA)
			colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_authority_score'	}
		if ('var_node_authority_score' %in% inputTerms)		{networkStats	=	cbind(networkStats, NA)
			colnames(networkStats)[ncol(networkStats)]	<-	'var_node_authority_score'	}	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# Average nearest neighbor degree: Calculate the average nearest neighbor degree
	# of the given vertices and the same quantity in the function of vertex degree
	# knn(g, vids = V(g), 	mode = c("all", "out", "in", "total"), 
	# 									neighbor.degree.mode = c("all", "out", "in", "total"))
	# mode 	Character constant to indicate the type of neighbors to consider in directed graphs.
	#			"out" considers out-neighbors, "in" considers in-neighbors and "all" ignores edge directions.
	# neighbor.degree.mode
	# 		The type of degree to average in directed graphs. "out" averages out-degrees, 
	#		"in" averages in-degrees and "all" ignores edge directions for the degree calculation.
	# knn(g, vids = V(g), 	mode = c("all", "out", "in", "total"), 
	# 									neighbor.degree.mode = c("all", "out", "in", "total"))
	if (('mean_node_KNNmeandegree_all' %in% inputTerms) | ('var_node_KNNmeandegree_all' %in% inputTerms)){
		out 	=	knn(g, mode = "all", neighbor.degree.mode = "all")$knn	}
	if ('mean_node_KNNmeandegree_all' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_KNNmeandegree_all'	}
	if ('var_node_KNNmeandegree_all' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_KNNmeandegree_all'	}
	if (('mean_node_KNNmeandegree_out' %in% inputTerms) | ('var_node_KNNmeandegree_out' %in% inputTerms)){
		out		= 	knn(g, mode = "out", neighbor.degree.mode = "out")$knn				}
	if ('mean_node_KNNmeandegree_out' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_KNNmeandegree_out'	}
	if ('var_node_KNNmeandegree_out' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_KNNmeandegree_out'	}
	if (('mean_node_KNNmeandegree_in' %in% inputTerms) | ('var_node_KNNmeandegree_in' %in% inputTerms)){
		out		= 	knn(g, mode = "in", neighbor.degree.mode = "in")$knn				}
	if ('mean_node_KNNmeandegree_in' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_KNNmeandegree_in'	}
	if ('var_node_KNNmeandegree_in' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_KNNmeandegree_in'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	# The eccentricity of a vertex is its shortest path distance from the farthest other node in the graph (g).
	# The smallest eccentricity in a graph (g) is called its radius
	if ('radius_all' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, radius(g, mode = "all"))
		colnames(networkStats)[ncol(networkStats)]	<-	'radius_all'	}
	if ('radius_out' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, radius(g, mode = "out"))
		colnames(networkStats)[ncol(networkStats)]	<-	'radius_out'	}
	if ('radius_in' %in% inputTerms)		{
		networkStats	=	cbind(networkStats, radius(g, mode = "in"))
		colnames(networkStats)[ncol(networkStats)]	<-	'radius_in'		}
	#---------------------------------------------------------------------------------------------------------------------------------
	# Calculates the reciprocity of a directed graph (g).
	# reciprocity(g, ignore.loops = TRUE, mode = c("default", "ratio"))
	if ('reciprocity' %in% inputTerms)	{
		networkStats	=	cbind(networkStats, reciprocity(g, mode = "default"))
		colnames(networkStats)[ncol(networkStats)]	<-	'reciprocity'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# Calculate selected eigenvalues and eigenvectors of a (supposedly sparse) g.
	# spectrum(g, algorithm = c("arpack", ...),which = list(),options = arpack_defaults)
	# Commenting out command line below (unstable computations led to error messages):
	#if (isDirectedGraph == FALSE){networkStats	=	cbind(networkStats, max(spectrum(g)$values, na.rm=TRUE))}
	#colnames(networkStats)[ncol(networkStats)]	<-	'max_eig'
	# "max_eig" not computed (via spectrum() ) because ARPACK solver can be unstable and result in error messages.
	#---------------------------------------------------------------------------------------------------------------------------------
	# Transitivity measures the probability that the adjacent vertices of a vertex are connected.
	# This is sometimes also called the clustering coefficient.
	# transitivity(g, type = c(	"undirected", "global", "globalundirected", "localundirected", "local",
	# 								"average", "localaverage", "localaverageundirected", "barrat", "weighted"),
	# 								vids = NULL, weights = NULL, isolates = c("NaN", "zero") )
	if ('clustering_coef' %in% inputTerms)	{
		if (isDirectedGraph == FALSE)	{
			networkStats	=	cbind(networkStats, transitivity(g, type = "undirected"))	}
		if (isDirectedGraph == TRUE)	{
			networkStats	=	cbind(networkStats, transitivity(as.undirected(g,mode="collapse"),type= "undirected"))}
		colnames(networkStats)[ncol(networkStats)]	<-	'clustering_coef'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# Count the different induced subgraphs of three vertices in a graph (g).
	if (('triad_census_003' %in% inputTerms)|('triad_census_012' %in% inputTerms)|('triad_census_102' %in% inputTerms)|('triad_census_021D' %in% inputTerms) 
	   |('triad_census_021U' %in% inputTerms)|('triad_census_021C' %in% inputTerms)|('triad_census_111D' %in% inputTerms)|('triad_census_111U' %in% inputTerms) 
	   |('triad_census_030T' %in% inputTerms)|('triad_census_030C' %in% inputTerms)|('triad_census_201' %in% inputTerms)|('triad_census_120D' %in% inputTerms) 
	   |('triad_census_120U' %in% inputTerms)|('triad_census_120C' %in% inputTerms)|('triad_census_210' %in% inputTerms)|('triad_census_300' %in% inputTerms)){
		quiet(out <- triad_census(g)) } # Using quiet() because undirected graph g generates a warning message.
	if ('triad_census_003' 	%in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[1], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_003'	}
	if ('triad_census_012' 	%in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[2], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_012'	}
	if ('triad_census_102' 	%in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[3], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_102'	}
	if ('triad_census_021D' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[4], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_021D'	}
	if ('triad_census_021U' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[5], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_021U'	}
	if ('triad_census_021C' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[6], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_021C'	}
	if ('triad_census_111D' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[7], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_111D'	}
	if ('triad_census_111U' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[8], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_111U'	}
	if ('triad_census_030T' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[9], nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_030T'	}
	if ('triad_census_030C' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[10],nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_030C'	}
	if ('triad_census_201' 	%in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[11],nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_201'	}
	if ('triad_census_120D' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[12],nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_120D'	}
	if ('triad_census_120U' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[13],nrow = 1))				
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_120U'	}
	if ('triad_census_120C' %in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[14],nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_120C'	}
	if ('triad_census_210' 	%in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[15],nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_210'	}
	if ('triad_census_300' 	%in% inputTerms){networkStats	=	cbind(networkStats, matrix(out[16],nrow = 1))
										colnames(networkStats)[ncol(networkStats)]	<-	'triad_census_300'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# Count how many triangles each vertex is part of, in a graph (g), or just list the triangles of a graph.
	if (('mean_node_triangles' %in% inputTerms) | ('var_node_triangles' %in% inputTerms)){
		out 	=	count_triangles(g)	}
	if ('mean_node_triangles' %in% inputTerms){
		networkStats	=	cbind(networkStats, mean(out, na.rm = TRUE))
		colnames(networkStats)[ncol(networkStats)]	<-	'mean_node_triangles'	}
	if ('var_node_triangles' %in% inputTerms){
		networkStats	=	cbind(networkStats, sum((out - mean(out, na.rm = TRUE))^2, na.rm = TRUE) / sum(!is.na(c(out))))
		colnames(networkStats)[ncol(networkStats)]	<-	'var_node_triangles'	}
	#---------------------------------------------------------------------------------------------------------------------------------
	# More network statistics from igraph package (for example)
	# compares the community structure between two networks g and g1:
	# compare(communities(g), communities(g1),method = c("vi", "nmi", "split.join", "rand", "adjusted.rand"))
	# motifs(g, size = 3, cut.prob = rep(0, size))
	# Calculate scan statistics on a time series of graphs. This is done by calculating the local scan statistics 
	# for each graph (g) and each vertex, and then normalizing across the vertices and across the time steps:
	# scan_stat(graphs, tau = 1, ell = 0, locality = c("us", "them"), ...)
	# These functions calculates similarity scores for vertices based on their connection patterns.: 
	# similarity(g,vids = V(g),mode = c("all", "out", "in", "total"),loops = FALSE,method = c("jaccard", "dice", "invlogweighted"))
	# The above statistic seems infeasible for a network with a huge number of nodes.
	# Summing up the edge weights of the adjacent edges for each vertex:
	# strength(g, vids = V(g), mode = c("all", "out", "in", "total"), loops = TRUE, weights = NULL)
	# The above network statistic is redundant with degree() because I am
	# dichotomizing all weighted network graphs before analyzing them.
	# Subgraph centrality of a vertex measures the number of subgraphs a vertex participates in, weighting
	# them according to their size.
	# subgraph_centrality(g, diag = FALSE)
	# Above 'measure can only be calculated for small graphs.' (igraph user's manual)
	# vertex_connectivity(g, source = NULL, target = NULL, checks = TRUE)
	# Above statistic is vertex pair version of cohesion(), already calculated above.
	#networkStats 	= 	c(networkStats)
	
	networkStats	=	matrix(networkStats[match(inputTerms, colnames(networkStats))], nrow = 1)
	
return(networkStats)}
# For testing networkStats_igraph(g) (using an undirected graph input g):
# g 	= 	sample_pa(n = 10, directed = FALSE)
# out 	= 	networkStats_igraph(g) 
# all(colnames(out) == termsUndirected_igraph) # TRUE

# For testing networkStats_igraph(g) (using a directed graph g):
# g 	= 	sample_pa(n = 10, directed = TRUE)
# out 	= 	networkStats_igraph(g) 	 
# all(colnames(out) == termsDirected_igraph) # TRUE



# This function computes the MPLE of ergm model
# based on an input network graph (g) in either the 'network' or 'igraph' format,
# and based on input network Terms (statistics) from the 'ergm' package.
# The default 'inputTerms_ergm = NULL' means to compute all the network statistics from the 'ergm' package.
MPLEergm =	function(g, inputTerms_ergm = NULL)	{
	is_network_format	=	is.network(g)# g in 'network' format? (or 'igraph' format?)
	if (is_network_format == TRUE)	{# then record properties of graph g:
		Nnodes 		=	network.size(g)
		isDirected 	=	get.network.attribute(g, "directed")
		isBipartite	=	get.network.attribute(g, "bipartite") > 0		
		if (isBipartite == TRUE)		{
			Nnodes1 =	get.network.attribute(g, "bipartite")# numeric (#'actors'; Nnodes = #'actors' + #'events').
			Nnodes2 =	Nnodes - Nnodes1 }}# list.network.attributes(graph_igraph)
	if (is_network_format == FALSE){# then convert g into 'network' format and record g's properties:
		Nnodes 		=	vcount(g)
		isDirected	=	is_directed(g)
		isBipartite =	is_bipartite(g) == 1
		ecount_g 	= 	ecount(g)
		if (isBipartite == FALSE)	{
			if (ecount_g	> 	0){g = network(as_edgelist(g),directed = isDirected, matrix.type = "edgelist")}
			if (ecount_g	==	0){g = network.initialize(Nnodes, directed = isDirected)}	}
		if (isBipartite == TRUE)	{
			Types 	=	vertex_attr(g, 'type')
			Nnodes1	=	sum(Types == FALSE)
			Nnodes2	=	sum(Types == TRUE)
			if (ecount_g > 	0){g	=	network(as_edgelist(g),directed = isDirected, bipartite = Nnodes1, matrix.type = "edgelist")}
			if (ecount_g ==	0){g	=	network.initialize(Nnodes, directed = isDirected, bipartite = Nnodes1)}	}
			g		=	set.network.attribute(g, 'n', Nnodes)		}
	# Extract relevant terms (ergmTerms) from 'ergm' package:
	UndirectedNonbipartite	=	(isDirected == FALSE) & (isBipartite == FALSE)
	UndirectedBipartite		=	(isDirected == FALSE) & (isBipartite == TRUE )
	DidInputTerms			=	!is.null(inputTerms_ergm)
	if (DidInputTerms == FALSE)	{
		if (UndirectedNonbipartite)	{inputTerms_ergm = TermsUndirected_ergm			}
		if (UndirectedBipartite)	{inputTerms_ergm = TermsUndirectedBipartite_ergm	}
		if (isDirected)				{inputTerms_ergm = TermsDirected_ergm				}		}	
	nTerms 			= 	length(inputTerms_ergm)
	MPLEout			= 	c();
	varMPLEout 		= 	c();
	#TermsNew = c()
	for (k in 1 : nTerms)	{
		modelFormula=	as.formula(paste(c('g ~ offset(edges) +', paste(inputTerms_ergm[k], collapse = ' + ')), collapse = ''))
		quiet(out 	<-	ergm(modelFormula, estimate = "MPLE", offset.coef = log(1 / Nnodes)))
		Coef		=	coef(out)[-1];		varCoef	=	diag(vcov(out))[-1]
		if (length(names(Coef)) > 0 ){names(Coef)[1] = inputTerms_ergm[k]}
		if (length(names(Coef)) == 0){
			names(Coef)[1] = inputTerms_ergm[k]
			if (length(Coef) > 1){names(Coef)[2 : length(Coef)] =  paste(inputTerms_ergm[k], '_', 2 : length(Coef), sep = '')}	}
			MPLEout		= 	c(MPLEout, Coef);
			varMPLEout	=	c(varMPLEout, varCoef)	}
	MPLEout				=	rbind(MPLEout, varMPLEout)
	rownames(MPLEout)	=	c('MPLE', 'varMPLE')
return(list(MPLE = MPLEout, ergmTerms = inputTerms_ergm))	}



# This function computes the MPLE of ergm model
# based on an input network graph (g) in either the 'network' or 'igraph' format,
# and based on input network Terms (statistics) from the 'igraph' package.
# The default 'inputTerms_igraph = NULL' means to compute all the network statistics from the 'igraph' package.
MPLEigraph = function(g, inputTerms_igraph = NULL)	{
	is_network		=	is.network(g)# g in 'network' format? (or 'igraph' format?)
	if (is_network == TRUE)	{# convert g into 'igraph' format and get its properties:
		Nnodes 		=	network.size(g)
		isDirected 	=	get.network.attribute(g, "directed")
		isBipartite	=	get.network.attribute(g, "bipartite") > 0
		if (isBipartite == FALSE){
			g		=	graph_from_edgelist(as.edgelist(g), directed = isDirected)
			g		=	add_vertices(g, Nnodes - vcount(g))}
		if (isBipartite == TRUE){
			Nnodes1 =	get.network.attribute(g, "bipartite")# numeric (#'actors'; Nnodes = #'actors' + #'events').
			Nnodes2 =	Nnodes - Nnodes1 # list.network.attributes(graph_igraph)
			Types  	=	ifelse( (1 : (Nnodes1 + Nnodes2)) <= Nnodes1, 0, 1)
			g		=	t(as.edgelist(g))
			g		=	make_bipartite_graph(types = Types, edges = g, directed = isDirected)}}
	if (is_network == FALSE)	{# if g in 'igraph' format, then just get its properties 
		Nnodes 		=	vcount(g)
		isDirected	=	is_directed(g)
		isBipartite =	is_bipartite(g) == 1
		if (isBipartite == TRUE)	{
			Types 	=	vertex_attr(g, 'type')
			Nnodes1	=	sum(Types == FALSE)
			Nnodes2	=	sum(Types == TRUE)		}		}

	# Extract relevant terms (ergmTerms) from 'igraph' package:
	UndirectedNonbipartite	=	(isDirected == FALSE) & (isBipartite == FALSE)
	UndirectedBipartite		=	(isDirected == FALSE) & (isBipartite == TRUE )
	DidInputTerms			=	!is.null(inputTerms_igraph)
	if (DidInputTerms == FALSE)	{inputTerms_igraph	=  Terms_igraph}
	nTerms					=	length(inputTerms_igraph)

	# Create pairwise indices (IDs for each pair of node in the network:
	if (isBipartite == FALSE)	{allNodePairs	=	t(combn(1 : Nnodes, 2))}
	if (isBipartite == TRUE)	{allNodePairs 	= 	cbind(	sort(rep(1 : Nnodes1, Nnodes2)),
															rep((Nnodes1 + 1) : (Nnodes1 + Nnodes2), Nnodes1) )		}
	if (isDirected == TRUE) 	{allNodePairs 	=	rbind(allNodePairs, cbind(allNodePairs[,2] , allNodePairs[,1]))	}

	# Function for toggling the network statistics computed from the 'igraph' package.
	# (network statistics from the 'ergm' package can already toggled using the ergmMPLE() command)
	toggle	=	function(nodePair){
		# nodePair = matrix(c(1,3),nrow = 1); nodePair = matrix(c(allNodePairs[45,]),nrow = 1) # For testing
		y_ij		=	ifelse(are_adjacent(g, nodePair[1], nodePair[2]) == TRUE, 1, 0)
		if (y_ij == 1){
			gNew 	=	delete_edges(g, get.edge.ids(g, c(nodePair[1], nodePair[2])))
			DELTA	=	networkStats_igraph(g, inputTerms_igraph)   - 	networkStats_igraph(gNew, inputTerms_igraph)}# Difference in network statistics
		if (y_ij == 0){
			gNew	=	add_edges(g, c(nodePair[1], nodePair[2]))
			DELTA	=	networkStats_igraph(gNew, inputTerms_igraph)- 	networkStats_igraph(g, inputTerms_igraph)}# Difference in network statistics						
		DELTA 		= 	matrix(DELTA, nrow = 1)
	return(list(y_ij = y_ij, DELTA = DELTA))}

	# Set up glm() data from 'igraph' package:
	out 	=	apply(allNodePairs, 1, toggle)
	Y		=	matrix(sapply(out, "[[", "y_ij"),ncol = 1)
	if (nTerms > 	1){DELTA 	= 	t(sapply(out, "[[", "DELTA"))					}
	if (nTerms ==	1){DELTA 	= 	matrix(sapply(out, "[[", "DELTA"), ncol = 1)	}
	O		=	rep(log(1 / Nnodes), length(Y)) # Offsets for glm().

	# Compute network statistics:
	MPLEout	= 	c();  	varMPLEout = c();
	for (k in 1 : nTerms)	{
		OK	=	(!all(is.na(DELTA[,k]))) & (!any(is.infinite(DELTA[,k])))
		if (OK){         # Below, '- 1' means exclude intercept from model:
			quiet(out	<-	glm(Y ~ . - 1, data = as.data.frame(DELTA[, k]), family = "binomial"))
			Coef		=	coef(out);	varCoef =	vcov(out)
			names(Coef)	=	inputTerms_igraph[k]		}
		if (!OK){Coef = NA; varCoef = NA; names(Coef) =	inputTerms_igraph[k]}
		MPLEout			=	c(MPLEout,		Coef	)
		varMPLEout		=	c(varMPLEout,	varCoef	)}
	MPLEout			=	rbind(MPLEout,	varMPLEout)
	rownames(MPLEout)	=	c('MPLE', 	'varMPLE')
return(list(MPLE = MPLEout, Terms = inputTerms_igraph))	}