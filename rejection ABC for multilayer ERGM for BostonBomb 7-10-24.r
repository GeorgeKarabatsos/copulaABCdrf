rm(list = ls())
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
library('ks')	# for kde()
install.packages('MASS');			library('MASS')		# For mvrnorm().
install.packages('mvtnorm');		library('mvtnorm')	# Multivariate Normal and t Distributions

# Import the results of the copulaABCdrf (multilayer ERGM) analysis of the BostonBomb2013 dataset. 
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Real network data\\BostonBomb2013')
load(list.files(path = ".", pattern = "\\.RData$"))

euclideanDistances			=	 as.matrix(pdist(sy.table, Y = matrix(X.test, nrow = 1) ))
theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .03),]
# Above, retained 1%,2%, or 3% of thetas with summary statistics sy nearest to observed dataset summaries sx.		

postEtheta.RejectionABC		=	colMeans(theta.table.RejectionABC)
postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
#
# Command below only applied to data (#variables) of dimension up to 6:
# post.kde.RejectionABC		=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate
# Therefore, using product kernel density estimation:
for (k in 1:length(parameterNames)) { # k = 1;	# For testing
	if (k == 1) {post.kde.RejectionABC = matrix(NA, nrow = nrow(theta.table.RejectionABC), ncol =  length(parameterNames))}
		fit.RejectionABC			= 	kde1d(theta.table.RejectionABC[,k])
		post.kde.RejectionABC[,k]	=	dkde1d(theta.table.RejectionABC[,k], fit.RejectionABC)	}
post.kde.RejectionABC		=	apply(post.kde.RejectionABC, 1, prod)
postMODEtheta.RejectionABC	= 	theta.table.RejectionABC[which.max(post.kde.RejectionABC),]		
priorPDFs.RejectionABC		=	dmvnorm(theta.table.RejectionABC, mean = muPrior, sigma = SigmaPrior)
likelihoods.RejectionABC	=	post.kde.RejectionABC / priorPDFs.RejectionABC
MLEtheta.RejectionABC		=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]
postSDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, sd)
postQtheta.RejectionABC		= 	apply(theta.table.RejectionABC, 2, quantile, probs = c(.025, .25, .75, .975))
Results.rejABC				=	rbind(	postEtheta.RejectionABC, postMEDtheta.RejectionABC, 
											postMODEtheta.RejectionABC, MLEtheta.RejectionABC, 
											postSDtheta.RejectionABC, postQtheta.RejectionABC)
colnames(Results.rejABC)	=	parameterNames
rownames(Results.rejABC)	=	c('postE', 'postMED', 'postMODE', 'MLE', 'postSD', '2.5%', '25%', '75%', '97.5%')
t(round(Results.rejABC, 2)) # transposing table is in the paper