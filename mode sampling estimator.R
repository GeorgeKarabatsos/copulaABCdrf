# breaks		# values of theta (for each parameter)
# PI.theta	# Marginal posterior CDFs
# theta.table.keep	# Model parameter samples (thetas) from the Reference Table
# posteriorPDFs		=	dStudentcopula(u.keep,df=post.df,scale=post.Correlations) * apply(pi.thetaRT.keep,1,prod)

# plot(breaks[,1], PI.theta[,1])plot(breaks[,1], PI.theta[,1])
Nmode					=	200000
for (j in 1 : Nmode)	{	# j = 1	# for testing
	if (j == 1) {
		posteriorPDFs.thetaRand	=	matrix(NA, Nmode, 1)
		thetasRand				=	matrix(NA, Nmode, d)
		priorPDFs.thetaRand		=	matrix(NA, Nmode, 1)
		uRandom					=	rStudentcopula(n = Nmode, df = post.df, scale =  post.Correlations)
		options(warn = 2) # This option stops the loops whenever a warning happens (for loop checking)
		# options(warn = 0) # Default
	}
	for (k in 1 : d) {	# k = 1	# for testing
		if (k == 1)	{
			thetaRand				=	matrix(NA, 1, d);	
			marginalPDFs.thetaRand	=	matrix(NA, 1, d)
		}
		Distances		=	abs(uRandom[j, k] - PI.theta[, k])
		Ind				=	which(Distances == min(Distances))
		OK				=	length(PI.theta[Ind, k]) == length(PI.theta[Ind - 1, k])
		if (OK) { 
			marginalPDFs.thetaRand0	=	(PI.theta[Ind, k] - PI.theta[Ind - 1, k]) / (breaks[Ind, k] - breaks[Ind - 1, k])
			thetaRand0				=	breaks[Ind, k]
			Ind2						=	which(marginalPDFs.thetaRand0 == max(marginalPDFs.thetaRand0))[1]
			marginalPDFs.thetaRand[k]	=	marginalPDFs.thetaRand0[Ind2]
			thetaRand[k]				=	thetaRand0[Ind2]
		}
		if (!OK) { 
			marginalPDFs.thetaRand[k]	=	0
			thetaRand[k]				=	breaks[Ind[length(Ind)], k]
		}
		if (k == d)	{
			posteriorPDF.thetaRand0	=	dStudentcopula(uRandom[j,], df = post.df, scale = post.Correlations) * prod(marginalPDFs.thetaRand)
			priorPDF0				=	dmvnorm(thetaRand, mean = muPrior, sigma = SigmaPrior)
		}
	}
	thetasRand[j, ]			=	thetaRand
	posteriorPDFs.thetaRand[j]	=	posteriorPDF.thetaRand0
	priorPDFs.thetaRand[j]	=	priorPDF0
	cat("iteration:", j, " \n")
	flush.console()
}
posteriorMode.thetaRand	=	thetasRand[posteriorPDFs.thetaRand == max(posteriorPDFs.thetaRand), ]
likelihoods.thetaRand		=	posteriorPDFs.thetaRand / priorPDFs.thetaRand
MLE.thetaRand				=	thetasRand[likelihoods.thetaRand == max(likelihoods.thetaRand), ]
#posteriorMode.thetaRand
#MLE.thetaRand
cbind(c(postEtheta), posteriorMode.thetaRand, MLE.thetaRand)
plot(c(postEtheta), posteriorMode.thetaRand)
plot(c(postEtheta), MLE.thetaRand)
cor(c(postEtheta), MLE.thetaRand)
cor(c(postEtheta), posteriorMode.thetaRand)

> cbind(c(postEtheta), posteriorMode.thetaRand, MLE.thetaRand)
                  posteriorMode.thetaRand MLE.thetaRand
 [1,] -5.64529658              -2.9386667    -2.9386667
 [2,] -2.85242735               6.3229865     6.3229865
 [3,] -3.85331843               1.1959503     1.1959503
 [4,]  0.04690447              -0.2175613    -0.2175613
 [5,] -1.88353068               3.7887015     3.7887015
 [6,] -2.14589498               3.7674681     3.7674681
 [7,] -6.38946109              -3.1248589    -3.1248589
 [8,] -2.75526310               8.5099818     8.5099818
 [9,]  1.69536097              -3.5196857    -3.5196857
[10,]  0.23055389               1.3810432     1.3810432
[11,] -2.47551040               5.7856672     5.7856672
[12,] -2.03617468               3.2424265     3.2424265
[13,] -5.66631223              -3.0405638    -3.0405638
[14,] -3.21588605               6.7731753     6.7731753
[15,]  2.61623514              -2.5542355    -2.5542355
[16,] -1.13689347               2.9610693     2.9610693
[17,] -2.14775854               8.4858654     8.4858654
[18,] -1.33510632               5.6072737     5.6072737
> cor(c(postEtheta), MLE.thetaRand)
[1] 0.0008941075
> cor(c(postEtheta), posteriorMode.thetaRand)
[1] 0.0008941075
> 

#> colSums(theta.DRFweights.table > 0)
# [1]  432 1044  449 1906 1357  865  444  872  311 1951 1401 1114  446  707  588 1623 1750 1098

c(theta.table[,1],breaks[,1])
pi.theta[order(ord),k]

