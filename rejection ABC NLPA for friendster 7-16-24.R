rm(list = ls())
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
library('ks')	# for kde()
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('drf');	library('drf')	# for variable importance (selection) analysis of reference table.

# Import the results of the copulaABCdrf (multilayer ERGM) analysis of the BostonBomb2013 dataset. 
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Real network data\\friendster')
load('Barabasi-Albert friendster 2024-02-09 22_46_42.RData')

# Optionally, since there are 3 summary statistics relative to the two Price model parameters,
# use the top two most important predictors (summaries) in the Reference Table.
#drf.forest 	=	drf(X = sy.table, Y = theta.table)# using compute.variable.importance is costly to compute
#variable_importance(drf.forest)# Calculate a simple measure of ’importance’ for each feature.
# A simple weighted sum of how many times feature i was split on at each depth in the forest.
# Output (found that the first and second summary statistics are the most important.
#          [,1]
#[1,] 0.6625702
#[2,] 0.2160702
#[3,] 0.1213596
# variableImportance(drf.forest)# (VERY SLOW) compute an mmd-based variable importance for the drf fit.

# Below, ulsing all 3 summary statistics:
# euclideanDistances			=	as.matrix(pdist(sy.table, Y = matrix(X.test, nrow = 1) ))

# Below, using the most important summary statistics, #1 and #3:
 euclideanDistances			=	as.matrix(pdist(sy.table[, c(1,2)], Y = matrix(X.test[c(1,2)], nrow = 1) ))

theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .01),]
# Above, retained 1%,2%, or 3% of thetas with summary statistics sy nearest to observed dataset summaries sx.		

postEtheta.RejectionABC		=	colMeans(theta.table.RejectionABC)
postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
post.kde.RejectionABC		=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate
postMODEtheta.RejectionABC	= 	theta.table.RejectionABC[which.max(post.kde.RejectionABC),]		
priorPDFs.RejectionABC		=	apply(cbind(dunif(theta.table.RejectionABC[,1], min = 0, max = 3), dunif(theta.table.RejectionABC[,2], min = 0, max = .20)),1,prod)
likelihoods.RejectionABC	=	post.kde.RejectionABC / priorPDFs.RejectionABC
MLEtheta.RejectionABC		=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]
postSDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, sd)
postQtheta.RejectionABC		= 	apply(theta.table.RejectionABC, 2, quantile, probs = c(.025, .25, .75, .975))
Results.rejABC				=	rbind(	postEtheta.RejectionABC, postMEDtheta.RejectionABC, 
											postMODEtheta.RejectionABC, MLEtheta.RejectionABC, 
											postSDtheta.RejectionABC, postQtheta.RejectionABC)
colnames(Results.rejABC)	=	parameterNames
rownames(Results.rejABC)	=	c('postE', 'postMED', 'postMODE', 'MLE', 'postSD', '2.5%', '25%', '75%', '97.5%')
round(Results.rejABC, 3)


# Results based on all 3 summary statistics:
# > round(Results.rejABC, 3)
#         alpha     p
# postE    0.686 0.002
# postMED  0.638 0.002
# postMODE 0.421 0.001
# MLE      0.421 0.001
# postSD   0.437 0.002
# 2.5%     0.041 0.000
# 25%      0.333 0.001
# 75%      1.055 0.004
# 97.5%    1.438 0.005
# 
# > round(Results.rejABC, 4)
#           alpha      p
# postE    0.6859 0.0023
# postMED  0.6377 0.0021
# postMODE 0.4213 0.0015
# MLE      0.4213 0.0015
# postSD   0.4365 0.0015
# 2.5%     0.0411 0.0002
# 25%      0.3325 0.0011
# 75%      1.0551 0.0036
# 97.5%    1.4377 0.0049


# Results based on selecting the two most important summary statistics:
# > round(Results.rejABC, 3)
#          alpha     p
# postE    1.209 0.016
# postMED  1.202 0.017
# postMODE 1.252 0.012
# MLE      1.252 0.012
# postSD   0.124 0.008
# 2.5%     1.043 0.002
# 25%      1.136 0.011
# 75%      1.262 0.024
# 97.5%    1.370 0.028