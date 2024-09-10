rm(list = ls())
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
library('ks')	# for kde()
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('drf');	library('drf')	# for variable importance (selection) analysis of reference table.

# Import the results of the copulaABCdrf (multilayer ERGM) analysis of the BostonBomb2013 dataset. 
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Real network data\\cit-HepPh')
load('Price cit-HepPh 2024-01-28 15_58_51.RData')

# Optionally, since there are 3 summary statistics relative to the two Price model parameters,
# use the top two most important predictors (summaries) in the Reference Table.
#drf.forest 	=	drf(X = MPLE.table, Y = theta.table)# using compute.variable.importance is costly to compute
#variable_importance(drf.forest)# Calculate a simple measure of ’importance’ for each feature.
# A simple weighted sum of how many times feature i was split on at each depth in the forest.
# Output (found that the first and third summary statistics are the most important.
#          [,1]
#[1,] 0.14418114
#[2,] 0.06744952
#[3,] 0.78836933
# variableImportance(drf.forest)# (VERY SLOW) compute an mmd-based variable importance for the drf fit.

# Below, using only the most important summary statistics, #1 and #3:
#euclideanDistances	=	as.matrix(pdist(MPLE.table[, c(1,3)], Y = matrix(X.test[c(1,3)], nrow = 1) ))

# Below, using results based on all 3 summaries:
euclideanDistances	=	as.matrix(pdist(MPLE.table, Y = matrix(X.test, nrow = 1) ))


theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .01),]
# Above, retained 1%,2%, or 3% of thetas with summary statistics sy nearest to observed dataset summaries sx.		

postEtheta.RejectionABC		=	colMeans(theta.table.RejectionABC)
postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
post.kde.RejectionABC		=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate
postMODEtheta.RejectionABC	= 	theta.table.RejectionABC[which.max(post.kde.RejectionABC),]		
priorPDFs.RejectionABC		=	apply(cbind(dunif(theta.table.RejectionABC[,1], min = 0.9, max = 1.1), dunif(theta.table.RejectionABC[,2], min = 0, max = .20)),1,prod)
likelihoods.RejectionABC	=	post.kde.RejectionABC / priorPDFs.RejectionABC
MLEtheta.RejectionABC		=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]
postSDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, sd)
postQtheta.RejectionABC		= 	apply(theta.table.RejectionABC, 2, quantile, probs = c(.025, .25, .75, .975))
Results.rejABC				=	rbind(	postEtheta.RejectionABC, postMEDtheta.RejectionABC, 
											postMODEtheta.RejectionABC, MLEtheta.RejectionABC, 
											postSDtheta.RejectionABC, postQtheta.RejectionABC)
colnames(Results.rejABC)	=	parameterNames
rownames(Results.rejABC)	=	c('postE', 'postMED', 'postMODE', 'MLE', 'postSD', '2.5%', '25%', '75%', '97.5%')
round(Results.rejABC, 2)



# Results based on all 3 summary statistics:

# > round(Results.rejABC, 2)
#           k_0    p
# postE    1.01 0.01
# postMED  1.02 0.01
# postMODE 1.02 0.01
# MLE      1.02 0.01
# postSD   0.06 0.00
# 2.5%     0.90 0.00
# 25%      0.97 0.01
# 75%      1.05 0.01
# 97.5%    1.09 0.01


# Results based on the 2 most important summary statistics:
#> round(Results.rejABC, 2)
#          k_0    p
# postE    1.01 0.01
# postMED  1.02 0.01
# postMODE 1.02 0.01
# MLE      1.02 0.01
# postSD   0.06 0.00
# 2.5%     0.90 0.00
# 25%      0.97 0.01
# 75%      1.05 0.01
# 97.5%    1.09 0.01


