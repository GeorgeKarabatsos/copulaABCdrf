rm(list = ls())
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
library('ks')	# for kde()

# Access the output data from the copulaABCdrf analysis of twisted normal model:
setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Twisted normal model')
RDataFiles		=	list.files(path = ".", pattern = "\\.RData$")
load(RDataFiles[1])

# Do the rejection ABC analysis:
euclideanDistances		=	 as.matrix(pdist(sy.table, Y = t(x)))
theta.table.RejectionABC	=	theta.table[euclideanDistances <= quantile(euclideanDistances, probs = .03),]
# Above, retained 3% of thetas with summary statstics sy nearest to observed dataset summaries sx.
postEtheta.RejectionABC	=	colMeans(theta.table.RejectionABC)
postMEDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, median)
b						=	0.1	# Controls prior correlation between theta[1] and theta[2]
priorPDFs.RejectionABC	= 	exp( -((theta.table.RejectionABC[,1]^2) / 200) - ((theta.table.RejectionABC[,2] - (b*(theta.table.RejectionABC[,1]^2)) + (100*b))^2)/2 - rowSums(theta.table.RejectionABC[,3:d]^2) )	
#post.kde300d.RejectionABC	=	kde(x=theta.table.RejectionABC, eval.points=theta.table.RejectionABC, positive=FALSE, unit.interval=FALSE)$estimate
# Note on above commented out command line:
# In order To estimate the posterior mode and MLE, using a different
# approach to multivariate kernel density estimation using kde1d(), below,
# (as the original approach using kde() used in other simulations 
# did not lead to positive definite results)
postSDtheta.RejectionABC	=	apply(theta.table.RejectionABC, 2, sd)

# Do weighted one-sample (two-tailed) Kolmogorov Smirnov tests (Monahan,2011,p.358)
KS.Dn.RejectionABC		=	c()
KS.stat.RejectionABC		=	c()
postQtheta025.RejectionABC	=	c()
postQtheta975.RejectionABC	=	c()
postQtheta25.RejectionABC	=	c()
postQtheta75.RejectionABC	=	c()
post.kde300d.RejectionABC	=	matrix(NA, nrow = nrow(theta.table.RejectionABC), ncol = p)
#
for (k in 1:p) { # k = 1;	# For testing
	theta.with.weights.RejectionABC	=	cbind(theta.table.RejectionABC[,k], rep(1/length(theta.table.RejectionABC[,k]),length(theta.table.RejectionABC[,k])))
   	# above using rep(1/length(theta.table.RejectionABC[,k]),length(theta.table.RejectionABC[,k])) in place of theta.DRFweights.table[,k]
	sorted.theta.with.weights.RejectionABC	=	theta.with.weights.RejectionABC[order(theta.with.weights.RejectionABC[,1]), ]
	# Univariate local-polynomial (degree =2) likelihood kernel density estimation (from MCMC samples)
	# using plug-in bandwidth of Sheather and Jones (1991) :
	fit 			= 	kde1d(Samples.theta[,k])# Samples.theta are samples from exact posterior distribution.
	F 			= 	pkde1d(sorted.theta.with.weights.RejectionABC[,1], fit)#cdf estimates of true marginal posterior distribution (from MCMC samples).
	Fk			=	cumsum(sorted.theta.with.weights.RejectionABC[,2])
	Fk.minus.1	=	c(0, Fk[1 : (length(Fk) - 1) ] )
	Dn 			= 	max(apply(cbind(abs(Fk - F), abs(F - Fk.minus.1)), 1, max))# Kolmogorov-Smirnov statistic
	S.2			=	sum(sorted.theta.with.weights.RejectionABC[,2] ^ 2)
	KSstat		=	(1 / S.2) * Dn
	# The critical values of KSstat to remember are Q(1.358) = .95 and Q(1.628) = .99.(for 2-tailed KS test)
	KS.Dn.RejectionABC[k]		=	Dn
	KS.stat.RejectionABC[k]	=	KSstat
	#
	# Calculate 95%(50%) credible interval:
	postQtheta025.RejectionABC[k]	=	quantile(theta.table.RejectionABC[,k], .025)
	postQtheta975.RejectionABC[k]	=	quantile(theta.table.RejectionABC[,k], .975)
	postQtheta25.RejectionABC[k]	=	quantile(theta.table.RejectionABC[,k], .25)
	postQtheta75.RejectionABC[k]	=	quantile(theta.table.RejectionABC[,k], .75)
	#
	# Estimated posterior densities from rejection ABC:
	fit.RejectionABC				= 	kde1d(theta.table.RejectionABC[,k])
	post.kde300d.RejectionABC[,k]	=	dkde1d(theta.table.RejectionABC[,k], fit.RejectionABC)
	if (k == p)	{post.kde300d.RejectionABC		=	apply(post.kde300d.RejectionABC, 1, prod)}
}

postMODEtheta.RejectionABC	= 	theta.table.RejectionABC[which.max(post.kde300d.RejectionABC),]
likelihoods.RejectionABC	=	post.kde300d.RejectionABC / priorPDFs.RejectionABC
MLEtheta.RejectionABC		=	theta.table.RejectionABC[which.max(likelihoods.RejectionABC),]

ResultsTable				=	rbind(	postEtheta.RejectionABC,
									postMEDtheta.RejectionABC,
									postMODEtheta.RejectionABC,
									MLEtheta.RejectionABC,
									postSDtheta.RejectionABC,
									postQtheta025.RejectionABC,
									postQtheta975.RejectionABC,
									postQtheta25.RejectionABC,
									postQtheta75.RejectionABC,
									KS.Dn.RejectionABC,
									KS.stat.RejectionABC)
rownames(ResultsTable) 	=	c(		'postE.RejABC',
									'postMED.RejABC',
									'postMODE.RejABC',
									'MLE.RejABC',
									'postSD.RejABC',
									'post2.5%.RejABC',
									'post97.5%.RejABC',
									'postQtheta25.RejABC',
									'postQtheta75.RejABC',
									'KS.Dn.RejABC',
									'KS.stat.RejABC')
colnames(ResultsTable)	=	 paste('theta', 1:300, sep = '')

ResultsTable.3to300minmax 				=	cbind(ResultsTable[,1:2], apply(ResultsTable[,3:300],1,min), apply(ResultsTable[,3:300],1,max)) 
colnames(ResultsTable.3to300minmax)	=	c(paste('theta',1:2,sep=''), 'theta3to300min', 'theta3to300max')

round(ResultsTable.3to300minmax, 2)
round(ResultsTable, 2)


# copulaABCdrf posterior quantile (min and max) results for theta_3,...,theta_300
#> round(cbind(apply(postQtheta[,3:300],1,min), apply(postQtheta[,3:300],1,max)),2)
#         min   max
# 2.5%  -2.44 -0.87
# 25%   -0.94 -0.16
# 50%   -0.36  0.31
# 75%    0.17  0.82
# 97.5%  0.78  2.27

# Exact posterior quantile (min and max) results for theta_3,...,theta_300 
# Q.exact	=	apply(Samples.theta[,3:300], 2, quantile, probs = c(0.025, .05, .25, .50, .75, .975))
# round(cbind(apply(Q.exact,1,min), apply(Q.exact,1,max)),2)
#         min   max
# 2.5%  -1.44 -1.32
# 5%    -1.22 -1.12
# 25%   -0.50 -0.45
# 50%   -0.02  0.02
# 75%    0.45  0.51
# 97.5%  1.33  1.44