rm(list = ls())
install.packages('pdist');	library('pdist') # To compute Euclidean distances for rejection ABC
install.packages('kde1d');	library('kde1d')	# for kernel CDF() used in Kolmogorov Smirnov test
install.packages('MASS');	library('MASS')		# For mvrnorm().
install.packages('mvtnorm');library('mvtnorm')	# Multivariate Normal and t Distributions
install.packages('nvmix');	library('nvmix')		# For t copula.
install.packages('kde1d');	library('kde1d')	# for kernel PDF and CDF()
simStudy = 21

for (simStudy in 1 : 62)	{ # simStudy = 8; simStudy = 10 ; simStudy = 6; simStudy = 9;# For testing.
	if (simStudy == 1)	{
		simConditionNames	=	c(		'PNmix n=100,n_sim=10',
									'PNmix n=100,n_sim=25',
									'PNmix n=100,n_sim=33',
									'PNmix n=100,n_sim=50',
									'PNmix n=100,n_sim=66',
									'PNmix n=100,n_sim=75',
									'PNmix n=100,n_sim=90',
									'PNmix n=100,n_sim=100',
									'GaussMultimodPstn=n_sim=4',
									'TwGaus,n=1',
									'ERGM, n=50 ,n_sim=5',
									'ERGM, n=50 ,n_sim=16',
									'ERGM, n=50 ,n_sim=25',
									'ERGM, n=50 ,n_sim=37',
									'ERGM, n=50 ,n_sim=50',
									'ERGM, n=300,n_sim=30',
									'ERGM, n=300,n_sim=100',
									'ERGM, n=300,n_sim=150',
									'ERGM, n=300,n_sim=225',
									'ERGM, n=300,n_sim=300',
									'ERGMm,n=100,n_sim=10',
									'ERGMm,n=100,n_sim=20',
									'Price,n=50 ,n_sim=5',
									'Price,n=50 ,n_sim=16',
									'Price,n=50 ,n_sim=25',
									'Price,n=50 ,n_sim=37',
									'Price,n=50 ,n_sim=50',
									'Price,n=300,n_sim=30',
									'Price,n=300,n_sim=100',
									'Price,n=300,n_sim=150',
									'Price,n=300,n_sim=225',
									'Price,n=300,n_sim=300',
									'NLPA, n=50 ,n_sim=5',
									'NLPA, n=50 ,n_sim=16',
									'NLPA, n=50 ,n_sim=25',
									'NLPA, n=50 ,n_sim=37',
									'NLPA, n=50 ,n_sim=50',
									'NLPA, n=300,n_sim=30',
									'NLPA, n=300,n_sim=100',
									'NLPA, n=300,n_sim=150',
									'NLPA, n=300,n_sim=225',
									'NLPA, n=300,n_sim=300',
									'DMC,  n=50 ,n_sim=5',
									'DMC,  n=50 ,n_sim=16',
									'DMC,  n=50 ,n_sim=25',
									'DMC,  n=50 ,n_sim=37',
									'DMC,  n=50 ,n_sim=50',
									'DMC,  n=300,n_sim=30',
									'DMC,  n=300,n_sim=100',
									'DMC,  n=300,n_sim=150',
									'DMC,  n=300,n_sim=225',
									'DMC,  n=300,n_sim=300',
									'DMR,  n=50 ,n_sim=5',
									'DMR,  n=50 ,n_sim=16',
									'DMR,  n=50 ,n_sim=25',
									'DMR,  n=50 ,n_sim=37',
									'DMR,  n=50 ,n_sim=50',
									'DMR,  n=300,n_sim=30',
									'DMR,  n=300,n_sim=100',
									'DMR,  n=300,n_sim=150',
									'DMR,  n=300,n_sim=225',
									'DMR,  n=300,n_sim=300')
		postMODEtheta.NEW1.simulations.mean			= 	matrix(NA, nrow = 62, ncol = 18)
		postMODEtheta.NEW1.simulations.sd				= 	matrix(NA, nrow = 62, ncol = 18)
		MLEtheta.NEW1.simulations.mean					= 	matrix(NA, nrow = 62, ncol = 18)
		MLEtheta.NEW1.simulations.sd					= 	matrix(NA, nrow = 62, ncol = 18)
		postMODEtheta.NEW2.simulations.mean			= 	matrix(NA, nrow = 62, ncol = 18)
		postMODEtheta.NEW2.simulations.sd				= 	matrix(NA, nrow = 62, ncol = 18)
		MLEtheta.NEW2.simulations.mean					= 	matrix(NA, nrow = 62, ncol = 18)
		MLEtheta.NEW2.simulations.sd					= 	matrix(NA, nrow = 62, ncol = 18)
		#
		post.df.NEW.simulations.mean					=	matrix(NA, nrow = 62, ncol = 1)
		post.df.NEW.simulations.sd						=	matrix(NA, nrow = 62, ncol = 1)
		post.scale.NEW.simulations.mean				=	matrix(NA, nrow = 62, ncol = 18*18)
		post.scale.NEW.simulations.sd					=	matrix(NA, nrow = 62, ncol = 18*18)
		#
		MAE.postMODEtheta.NEW1.simulations				= 	matrix(NA, nrow = 62, ncol = 18)
		MSE.postMODEtheta.NEW1.simulations				= 	matrix(NA, nrow = 62, ncol = 18)
		MAE.MLEtheta.NEW1.simulations					= 	matrix(NA, nrow = 62, ncol = 18)
		MSE.MLEtheta.NEW1.simulations					= 	matrix(NA, nrow = 62, ncol = 18)
		MAE.postMODEtheta.NEW2.simulations				= 	matrix(NA, nrow = 62, ncol = 18)
		MSE.postMODEtheta.NEW2.simulations				= 	matrix(NA, nrow = 62, ncol = 18)
		MAE.MLEtheta.NEW2.simulations					= 	matrix(NA, nrow = 62, ncol = 18)
		MSE.MLEtheta.NEW2.simulations					= 	matrix(NA, nrow = 62, ncol = 18)
		#
		rownames(postMODEtheta.NEW1.simulations.mean)	=	simConditionNames
		rownames(postMODEtheta.NEW1.simulations.sd)		=	simConditionNames
		rownames(MLEtheta.NEW1.simulations.mean)		=	simConditionNames
		rownames(MLEtheta.NEW1.simulations.sd)			=	simConditionNames
		rownames(MAE.postMODEtheta.NEW1.simulations)	=	simConditionNames
		rownames(MSE.postMODEtheta.NEW1.simulations)	=	simConditionNames
		rownames(MAE.MLEtheta.NEW1.simulations)			=	simConditionNames
		rownames(MSE.MLEtheta.NEW1.simulations)			=	simConditionNames
		rownames(postMODEtheta.NEW2.simulations.mean)	=	simConditionNames
		rownames(postMODEtheta.NEW2.simulations.sd)		=	simConditionNames
		rownames(MLEtheta.NEW2.simulations.mean)		=	simConditionNames
		rownames(MLEtheta.NEW2.simulations.sd)			=	simConditionNames
		rownames(MAE.postMODEtheta.NEW2.simulations)	=	simConditionNames
		rownames(MSE.postMODEtheta.NEW2.simulations)	=	simConditionNames
		rownames(MAE.MLEtheta.NEW2.simulations)			=	simConditionNames
		rownames(MSE.MLEtheta.NEW2.simulations)			=	simConditionNames
		#
		colnames(postMODEtheta.NEW1.simulations.mean)	=	paste('th', 1 : 18, sep = '')
		colnames(postMODEtheta.NEW1.simulations.sd)		=	paste('th', 1 : 18, sep = '')
		colnames(MLEtheta.NEW1.simulations.mean)		=	paste('th', 1 : 18, sep = '')
		colnames(MLEtheta.NEW1.simulations.sd)			=	paste('th', 1 : 18, sep = '')
		colnames(MAE.postMODEtheta.NEW1.simulations)	=	paste('th', 1 : 18, sep = '')
		colnames(MSE.postMODEtheta.NEW1.simulations)	=	paste('th', 1 : 18, sep = '')
		colnames(MAE.MLEtheta.NEW1.simulations)			=	paste('th', 1 : 18, sep = '')
		colnames(MSE.MLEtheta.NEW1.simulations)			=	paste('th', 1 : 18, sep = '')
		#
		colnames(postMODEtheta.NEW2.simulations.mean)	=	paste('th', 1 : 18, sep = '')
		colnames(postMODEtheta.NEW2.simulations.sd)		=	paste('th', 1 : 18, sep = '')
		colnames(MLEtheta.NEW2.simulations.mean)		=	paste('th', 1 : 18, sep = '')
		colnames(MLEtheta.NEW2.simulations.sd)			=	paste('th', 1 : 18, sep = '')
		colnames(MAE.postMODEtheta.NEW2.simulations)	=	paste('th', 1 : 18, sep = '')
		colnames(MSE.postMODEtheta.NEW2.simulations)	=	paste('th', 1 : 18, sep = '')
		colnames(MAE.MLEtheta.NEW2.simulations)			=	paste('th', 1 : 18, sep = '')
		colnames(MSE.MLEtheta.NEW2.simulations)			=	paste('th', 1 : 18, sep = '')
		#
		rownames(post.df.NEW.simulations.mean)			=	simConditionNames
		rownames(post.df.NEW.simulations.sd)			=	simConditionNames
		rownames(post.scale.NEW.simulations.mean)		=	simConditionNames
		rownames(post.scale.NEW.simulations.sd)			=	simConditionNames
		#
		postMODEtheta.OLD.simulations.mean				= 	matrix(NA, nrow = 62, ncol = 18)
		postMODEtheta.OLD.simulations.sd				= 	matrix(NA, nrow = 62, ncol = 18)
		MLEtheta.OLD.simulations.mean					= 	matrix(NA, nrow = 62, ncol = 18)
		MLEtheta.OLD.simulations.sd					= 	matrix(NA, nrow = 62, ncol = 18)
		post.df.OLD.simulations.mean					=	matrix(NA, nrow = 62, ncol = 1)
		post.df.OLD.simulations.sd						=	matrix(NA, nrow = 62, ncol = 1)
		post.scale.OLD.simulations.mean				=	matrix(NA, nrow = 62, ncol = 18*18)
		post.scale.OLD.simulations.sd					=	matrix(NA, nrow = 62, ncol = 18*18)
		MAE.postMODEtheta.OLD.simulations				= 	matrix(NA, nrow = 62, ncol = 18)
		MSE.postMODEtheta.OLD.simulations				= 	matrix(NA, nrow = 62, ncol = 18)
		MAE.MLEtheta.OLD.simulations					= 	matrix(NA, nrow = 62, ncol = 18)
		MSE.MLEtheta.OLD.simulations					= 	matrix(NA, nrow = 62, ncol = 18)
		rownames(postMODEtheta.OLD.simulations.mean)	=	simConditionNames
		rownames(postMODEtheta.OLD.simulations.sd)		=	simConditionNames
		rownames(MLEtheta.OLD.simulations.mean)			=	simConditionNames
		rownames(MLEtheta.OLD.simulations.sd)			=	simConditionNames
		rownames(MAE.postMODEtheta.OLD.simulations)		=	simConditionNames
		rownames(MSE.postMODEtheta.OLD.simulations)		=	simConditionNames
		rownames(MAE.MLEtheta.OLD.simulations)			=	simConditionNames
		rownames(MSE.MLEtheta.OLD.simulations)			=	simConditionNames
		colnames(postMODEtheta.OLD.simulations.mean)	=	paste('th', 1 : 18, sep = '')
		colnames(postMODEtheta.OLD.simulations.sd)		=	paste('th', 1 : 18, sep = '')
		colnames(MLEtheta.OLD.simulations.mean)			=	paste('th', 1 : 18, sep = '')
		colnames(MLEtheta.OLD.simulations.sd)			=	paste('th', 1 : 18, sep = '')
		colnames(MAE.postMODEtheta.OLD.simulations)		=	paste('th', 1 : 18, sep = '')
		colnames(MSE.postMODEtheta.OLD.simulations)		=	paste('th', 1 : 18, sep = '')
		colnames(MAE.MLEtheta.OLD.simulations)			=	paste('th', 1 : 18, sep = '')
		colnames(MSE.MLEtheta.OLD.simulations)			=	paste('th', 1 : 18, sep = '')
		rownames(post.df.OLD.simulations.mean)			=	simConditionNames
		rownames(post.df.OLD.simulations.sd)			=	simConditionNames
		rownames(post.scale.OLD.simulations.mean)		=	simConditionNames
		rownames(post.scale.OLD.simulations.sd)			=	simConditionNames
	}

# Poisson normal mixture model:
	if (simStudy == 1) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 10')	}
	if (simStudy == 2) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 25')	}
	if (simStudy == 3) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 33')	}
	if (simStudy == 4) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 50')	}
	if (simStudy == 5) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 66')	}
	if (simStudy == 6) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 75')	}
	if (simStudy == 7) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 90')	}
	if (simStudy == 8) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Poisson & normScaleMix\\n = 100, n.table = 100')	}
# Gaussian multimodal posterior:
	if (simStudy == 9) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Gaussian Multimodal Posterior')}
# Twisted normal model (300 dimensional):
	if (simStudy == 10){setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Twisted normal model')}
# 	ERGM:
	if (simStudy == 11) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 50 nodes\\n_sim = 5')}
	if (simStudy == 12) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 50 nodes\\n_sim = 16')}
	if (simStudy == 13) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 50 nodes\\n_sim = 25')}
	if (simStudy == 14) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 50 nodes\\n_sim = 37')}
	if (simStudy == 15) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 50 nodes\\n_sim = 50')}
	if (simStudy == 16) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 300 nodes\\n_sim = 30')}
	if (simStudy == 17) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 300 nodes\\n_sim = 100')}
	if (simStudy == 18) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 300 nodes\\n_sim = 150')}
	if (simStudy == 19) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 300 nodes\\n_sim = 225')}
	if (simStudy == 20) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\n = 300 nodes\\n_sim = 300')}
# ERGM multilayer BostonBomb2013 design:
	if (simStudy == 21){setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\Multilayer-valued ERGM\\n = 100 n_sim = 10')}
	if (simStudy == 22){setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\ERGM\\Multilayer-valued ERGM\\n = 100 n_sim = 20')}
#	Price (n = 50):
	if (simStudy == 23) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 50 nodes\\n_sim = 5')}
	if (simStudy == 24) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 50 nodes\\n_sim = 16')}
	if (simStudy == 25) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 50 nodes\\n_sim = 25')}
	if (simStudy == 26) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 50 nodes\\n_sim = 37')}
	if (simStudy == 27) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 50 nodes\\n_sim = 50')}
#	Price (n = 300):
	if (simStudy == 28) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 300 nodes\\n_sim = 30')}
	if (simStudy == 29) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 300 nodes\\n_sim = 100')}
	if (simStudy == 30) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 300 nodes\\n_sim = 150')}
	if (simStudy == 31) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 300 nodes\\n_sim = 225')}
	if (simStudy == 32) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Price\\n = 300 nodes\\n_sim = 300')}
#	NLPA:
	if (simStudy == 33) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 50 nodes\\n_sim = 5')}
	if (simStudy == 34) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 50 nodes\\n_sim = 16')}
	if (simStudy == 35) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 50 nodes\\n_sim = 25')}
	if (simStudy == 36) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 50 nodes\\n_sim = 37')}
	if (simStudy == 37) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 50 nodes\\n_sim = 50')}
	if (simStudy == 38) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 300 nodes\\n_sim = 30')}
	if (simStudy == 39) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 300 nodes\\n_sim = 100')}
	if (simStudy == 40) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 300 nodes\\n_sim = 150')}
	if (simStudy == 41) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 300 nodes\\n_sim = 225')}
	if (simStudy == 42) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\NLPA\\n = 300 nodes\\n_sim = 300')}
#	DMC:
	if (simStudy == 43) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 50 nodes\\n_sim = 5')}
	if (simStudy == 44) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 50 nodes\\n_sim = 16')}
	if (simStudy == 45) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 50 nodes\\n_sim = 25')}
	if (simStudy == 46) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 50 nodes\\n_sim = 37')}
	if (simStudy == 47) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 50 nodes\\n_sim = 50')}
	if (simStudy == 48) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 300 nodes\\n_sim = 30')}
	if (simStudy == 49) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 300 nodes\\n_sim = 100')}
	if (simStudy == 50) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 300 nodes\\n_sim = 150')}
	if (simStudy == 51) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 300 nodes\\n_sim = 225')}
	if (simStudy == 52) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMC\\n = 300 nodes\\n_sim = 300')}
#	DMR:
	if (simStudy == 53) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 50 nodes\\n_sim = 5')}
	if (simStudy == 54) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 50 nodes\\n_sim = 16')}
	if (simStudy == 55) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 50 nodes\\n_sim = 25')}
	if (simStudy == 56) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 50 nodes\\n_sim = 37')}
	if (simStudy == 57) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 50 nodes\\n_sim = 50')}
	if (simStudy == 58) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 300 nodes\\n_sim = 30')}
	if (simStudy == 59) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 300 nodes\\n_sim = 100')}
	if (simStudy == 60) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 300 nodes\\n_sim = 150')}
	if (simStudy == 61) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 300 nodes\\n_sim = 225')}
	if (simStudy == 62) {setwd('C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\DMR\\n = 300 nodes\\n_sim = 300')}

	RDataFiles	=	list.files(path = ".", pattern = "\\.RData$")

	postMODEtheta.OLD.simulations 			=	c()
	postMODEtheta.NEW1.simulations 		=	c()
	postMODEtheta.NEW2.simulations 		=	c()

	MLEtheta.OLD.simulations 				=	c()
	MLEtheta.NEW1.simulations 				=	c()
	MLEtheta.NEW2.simulations 				=	c()

	post.df.OLD.simulations				=	c()
	post.df.NEW.simulations				=	c()
	post.scale.OLD.simulations				=	c()
	post.scale.NEW.simulations				=	c()

	L1error.postMODEtheta.OLD.simulations	= 	c()
	L1error.postMODEtheta.NEW1.simulations	= 	c()
	L1error.postMODEtheta.NEW2.simulations	= 	c()

	L2error.postMODEtheta.OLD.simulations	= 	c()
	L2error.postMODEtheta.NEW1.simulations	= 	c()
	L2error.postMODEtheta.NEW2.simulations	= 	c()

	L1error.MLEtheta.OLD.simulations		= 	c()
	L1error.MLEtheta.NEW1.simulations		= 	c()
	L1error.MLEtheta.NEW2.simulations		= 	c()

	L2error.MLEtheta.OLD.simulations		= 	c()
	L2error.MLEtheta.NEW1.simulations		= 	c()
	L2error.MLEtheta.NEW2.simulations		= 	c()


	# Loop through the 10 replicas for the given network model, n nodes, and n_sim.
	for (fileNum in 1 : length(RDataFiles)){# fileNum = 1;  fileNum = 4 # For testing.

		load(RDataFiles[fileNum])

		# sy = MPLE.table used for ERGM(nonMultiLayer), Price, DMC models.  sy = sy.table used for NLPA, DMR models, and ERGM(MultiLayer).
		if ( 	((simStudy >= 11)  & (simStudy <= 20)) # ERGM(nonMultiLayer) models
			| 	((simStudy >= 23)  & (simStudy <= 32)) # Price model
			|	((simStudy >= 43)  & (simStudy <= 52))	# DMC model
			){sy.table	=	MPLE.table}
		#
		# --------------------------------------------------------------------------------------------
		# Old results are based on empirical density function (histogram-like):
		if (!is.matrix(posteriorMode)){
			postMODEtheta.OLD					=	posteriorMode	
			MLEtheta.OLD						=	MLE						}
		if (is.matrix(posteriorMode)){
			postMODEtheta.OLD					=	posteriorMode	[1,]
			MLEtheta.OLD						=	MLE[,1]					}
		postMODEtheta.OLD.simulations			=	rbind(postMODEtheta.OLD.simulations, 	postMODEtheta.OLD)
		MLEtheta.OLD.simulations				=	rbind(MLEtheta.OLD.simulations, 		MLEtheta.OLD)
		post.df.OLD.simulations				=	rbind(post.df.OLD.simulations, post.df)
		if ((simStudy != 21) &  (simStudy != 22)){
			post.scale.OLD.simulations			=	rbind(post.scale.OLD.simulations, c(post.Correlations))}
		if ((simStudy == 21) || (simStudy == 22)){
			post.scale.OLD.simulations			=	rbind(post.scale.OLD.simulations, c(post.scale))}
		L1error.postMODEtheta.OLD.simulations	=	rbind(L1error.postMODEtheta.OLD.simulations, abs(postMODEtheta.OLD	- truth)	)
		L2error.postMODEtheta.OLD.simulations	=	rbind(L2error.postMODEtheta.OLD.simulations, (postMODEtheta.OLD 	- truth)^2	)
		L1error.MLEtheta.OLD.simulations		=	rbind(L1error.MLEtheta.OLD.simulations, 	 abs(MLEtheta.OLD	- truth)	)
		L2error.MLEtheta.OLD.simulations		=	rbind(L2error.MLEtheta.OLD.simulations, 	 (MLEtheta.OLD		- truth)^2	)
		# --------------------------------------------------------------------------------------------
		#
		# I found that pi.thetaRT.keep (based on OLD results above)produced many zeros and 
		# u produced many zeros and ones leading to all zeros for the posteriorPDFs computations, below.
		# Therefore, I tried using kernel smoothing for marginal densities and inverse cdfs, below, which mitigated this issue.
		#
		# Compute new method for computing posterior mode and MLE (based on kernel univariate marginal posteriors):
		for (k in 1 : length(truth)){
			if (k==1){	pi.thetaRT.kernel	=	matrix(NA, nrow = N, ncol = length(truth))
						u.kernel			=	matrix(NA, nrow = N, ncol = length(truth))	
						fitKernelRT 		= 	c()													}
			#
			# For Poisson Normal mixture:
			if ((simStudy >= 1) & (simStudy <= 8)){if (k == 1){XMIN = 0; XMAX = NaN}; if (k == 2){XMIN = -10; XMAX = 10}}
			#if ((simStudy >= 1) & (simStudy <= 8)){if (k == 1){XMIN = 0; XMAX = NaN}; if (k == 2){XMIN = NaN; XMAX = NaN}}
			#if ((simStudy >= 1) & (simStudy <= 8)){if (k == 1){XMIN = NaN; XMAX = NaN}; if (k == 2){XMIN = NaN; XMAX = NaN}}
			#
			# For Gaussian multimodal posterior:
			if ( simStudy == 9)	{if (k == 2){XMIN = -4; XMAX = 4}; if (k != 2){XMIN = -3; XMAX = 3}}
			#if ( simStudy == 9)	{XMIN = NaN; XMAX = NaN}
			# For twisted normal model:
			if ( simStudy == 10){XMIN = NaN; XMAX = NaN}
			# For ERGM:
			if ((simStudy >= 11) & (simStudy <= 22)){XMIN = NaN; XMAX = NaN}
			# For Price n = 50:
			if ((simStudy >= 23) & (simStudy <= 27)){if (k == 1){XMIN = 0.9;  XMAX = 1.1 }; if (k == 2){XMIN = 0; XMAX = .10}}
			# For Price n = 300:
			if ((simStudy >= 28) & (simStudy <= 32)){if (k == 1){XMIN = 0.9;  XMAX = 1.1 }; if (k == 2){XMIN = 0; XMAX = .20}}
			# NLPA:
			if ((simStudy >= 33) & (simStudy <= 42)){if (k == 1){XMIN = 0  ;  XMAX = 3   }; if (k == 2){XMIN = 0; XMAX = .20}}
			# DMC or DMR:
			if ((simStudy >= 43) & (simStudy <= 62)){if (k == 1){XMIN = 0.15; XMAX = 0.35}; if (k == 2){XMIN = 0; XMAX = 1  }}
			#
			fitKernelRT[[k]] 		= 	kde1d(theta.table[,k], xmin = XMIN, xmax = XMAX, weights = N * theta.DRFweights.table[,k])
			#
			pi.thetaRT.kernel[,k]	=	dkde1d(theta.table[,k], fitKernelRT[[k]])
			u.kernel[,k]			=	pkde1d(theta.table[,k], fitKernelRT[[k]])
		}
		#"multiplied by N / (N + 1) to avoid evaluating the [copula] density at the edges of the unit square."
		# (from p. 499 of Genest & Neslehova, 2007, Astin Bulletin)
		# See also Genest et al. (1995) and Okhrin 2012, Ch.17, p.484, "Fitting High-Dimensional Copulae to Data"
		#u.kernel		=	u.kernel * ( N / ( N + 1 ) ) # Don't use. Can worsen MAE and MSE results according to simulation studies! 
		keep				=	apply((u.kernel > 0) & (u.kernel < 1), 1, all) & (apply(pi.thetaRT.kernel,1,prod)>0) # sum(keep)
		u.keep			=	u.kernel[keep,]
		theta.table.keep	=	theta.table[keep, ]
		pi.thetaRT.keep 	= 	pi.thetaRT.kernel[keep,]
		# Estimate copula parameters (df, correlation matrix) of Meta-t posterior distribution:
		OutCopulaFit  	=	fitStudentcopula(u.keep,fit.method = "EM-MLE", df.bounds = c(0.1, 1000), verbose = TRUE)
		post.df.NEW		=	OutCopulaFit$df
		post.scale.NEW	=	OutCopulaFit$scale
		post.df.NEW.simulations	=	c(post.df.NEW.simulations, 	post.df.NEW			)
		post.scale.NEW.simulations	=	rbind(post.scale.NEW.simulations, 	c(post.scale.NEW)		)
		#
		#
		posteriorPDFs.NEW		=	dStudentcopula(u.keep, df = post.df.NEW, scale = post.scale.NEW) * apply(pi.thetaRT.keep, 1, prod)
		postMODEtheta.NEW1	=	theta.table.keep[which.max(posteriorPDFs.NEW), 	]
		#
		# For Poisson Normal mixture:
		if ((simStudy >= 1) & (simStudy <= 8))	{priorPDFs	=	apply(cbind(dgamma(theta.table.keep[,1], shape = 1/2, rate = 0.1), dunif(theta.table.keep[,2], min = -10, max = 10)), 1, prod)}
		# For Gaussian multimodal posterior:
		if ( simStudy == 9					 )	{priorPDFs	=	apply(	cbind(	dunif(theta.table.keep[,1],min=-3,max=3),dunif(theta.table.keep[,2],min=-4,max=4),
																				dunif(theta.table.keep[,3],min=-3,max=3),dunif(theta.table.keep[,4],min=-3,max=3),
																				dunif(theta.table.keep[,5],min=-3,max=3)), 1, prod)	}
		# For twisted normal model:
		if ( simStudy == 10					 )	{priorPDFs	= 	exp( -((theta.table.keep[,1]^2) / 200) - ((theta.table.keep[,2] - (b*(theta.table.keep[,1]^2)) + (100*b))^2)/2 - rowSums(theta.table.keep[,3:d]^2) )	}
		# For ERGM:
		if ((simStudy >= 11) & (simStudy <= 22)){priorPDFs	=	dmvnorm(theta.table.keep, mean = muPrior, sigma = SigmaPrior)}
		# For Price n = 50:
		if ((simStudy >= 23) & (simStudy <= 27)){priorPDFs	=	apply(cbind(dunif(theta.table.keep[,1], min = 0.9, max = 1.1), 	dunif(theta.table.keep[,2], min = 0, max = .10)),1,prod)}	
		#  For Price n = 300:
		if ((simStudy >= 28) & (simStudy <= 32)){priorPDFs	=	apply(cbind(dunif(theta.table.keep[,1], min = 0.9, max = 1.1), 	dunif(theta.table.keep[,2], min = 0, max = .20)),1,prod)}
		# For NLPA: 
		if ((simStudy >= 33) & (simStudy <= 42)){priorPDFs	=	apply(cbind(dunif(theta.table.keep[,1], min = 0, max = 3), 		dunif(theta.table.keep[,2], min = 0, max = .20)),1,prod)}
		# For DMC or DMR: 
		if ((simStudy >= 43) & (simStudy <= 62)){priorPDFs	=	apply(cbind(dunif(theta.table.keep[,1], min = 0.15, max = .35), 	dunif(theta.table.keep[,2], min = 0, max = 1)),1,prod)}
		#
		likelihoods.NEW		=	posteriorPDFs.NEW / priorPDFs
		MLEtheta.NEW1			=	theta.table.keep[which.max(likelihoods.NEW),	]
		#
		#
		# Take random samples from fitted student copula:
		r.u.kernel		=	rStudentcopula(n = 100000, df = post.df.NEW, scale = post.scale.NEW)
		for (k in 1 : length(truth)){
			if (k==1){	pi.thetaRT.kernel	=	matrix(NA, nrow = 100000, ncol = length(truth))
						r.theta.kernel	=	matrix(NA, nrow = 100000, ncol = length(truth))	}
			r.theta.kernel[,k]	=	qkde1d(r.u.kernel[,k]		, 	fitKernelRT[[k]])
			pi.thetaRT.kernel[,k]	=	dkde1d(r.theta.kernel[,k]	, 	fitKernelRT[[k]])			}
		posteriorPDFs.NEW		=	dStudentcopula(r.u.kernel, df = post.df.NEW, scale = post.scale.NEW) * apply(pi.thetaRT.kernel, 1, prod)
		# Find posterior mode and MLE of theta using the sampled values:
		postMODEtheta.NEW2	=	r.theta.kernel[which.max(posteriorPDFs.NEW), 	]
		#
		# For Poisson Normal mixture:
		if ((simStudy >= 1) & (simStudy <= 8))	{priorPDFs	=	apply(cbind(dgamma(r.theta.kernel[,1], shape = 1/2, rate = 0.1), dunif(r.theta.kernel[,2], min = -10, max = 10)), 1, prod)}
		# For Gaussian multimodal posterior:
		if ( simStudy == 9					 )	{priorPDFs	=	apply(	cbind(	dunif(r.theta.kernel[,1],min=-3,max=3),dunif(r.theta.kernel[,2],min=-4,max=4),
																				dunif(r.theta.kernel[,3],min=-3,max=3),dunif(r.theta.kernel[,4],min=-3,max=3),
																				dunif(r.theta.kernel[,5],min=-3,max=3)), 1, prod)	}
		# For twisted normal model:
		if ( simStudy == 10					 )	{priorPDFs	= 	exp( -((r.theta.kernel[,1]^2) / 200) - ((r.theta.kernel[,2] - (b*(r.theta.kernel[,1]^2)) + (100*b))^2)/2 - rowSums(r.theta.kernel[,3:d]^2) )	}
		# For ERGM:
		if ((simStudy >= 11) & (simStudy <= 22)){priorPDFs	=	dmvnorm(r.theta.kernel, mean = muPrior, sigma = SigmaPrior)}
		# For Price n = 50:
		if ((simStudy >= 23) & (simStudy <= 27)){priorPDFs	=	apply(cbind(dunif(r.theta.kernel[,1], min = 0.9, max = 1.1), 	dunif(r.theta.kernel[,2], min = 0, max = .10)),1,prod)}	
		#  For Price n = 300:
		if ((simStudy >= 28) & (simStudy <= 32)){priorPDFs	=	apply(cbind(dunif(r.theta.kernel[,1], min = 0.9, max = 1.1), 	dunif(r.theta.kernel[,2], min = 0, max = .20)),1,prod)}
		# For NLPA: 
		if ((simStudy >= 33) & (simStudy <= 42)){priorPDFs	=	apply(cbind(dunif(r.theta.kernel[,1], min = 0, max = 3), 		dunif(r.theta.kernel[,2], min = 0, max = .20)),1,prod)}
		# For DMC or DMR: 
		if ((simStudy >= 43) & (simStudy <= 62)){priorPDFs	=	apply(cbind(dunif(r.theta.kernel[,1], min = 0.15, max = .35), dunif(r.theta.kernel[,2], min = 0, max = 1)),1,prod)}
		#
		likelihoods.NEW		=	posteriorPDFs.NEW / priorPDFs
		MLEtheta.NEW2			=	r.theta.kernel[which.max(likelihoods.NEW),	]
		#
		# Compare results:
		# rbind(postMODEtheta.NEW, 	postMODEtheta.OLD)
		# rbind(MLEtheta.NEW, 		MLEtheta.OLD)
		#
		# Store the results:
		postMODEtheta.NEW1.simulations			=	rbind(postMODEtheta.NEW1.simulations, 	postMODEtheta.NEW1)
		MLEtheta.NEW1.simulations				=	rbind(MLEtheta.NEW1.simulations, 		MLEtheta.NEW1)
		#
		postMODEtheta.NEW2.simulations			=	rbind(postMODEtheta.NEW2.simulations, 	postMODEtheta.NEW2)
		MLEtheta.NEW2.simulations				=	rbind(MLEtheta.NEW2.simulations, 		MLEtheta.NEW2)
		#
		L1error.postMODEtheta.NEW1.simulations	=	rbind(L1error.postMODEtheta.NEW1.simulations, 	abs(postMODEtheta.NEW1- truth)	)
		L2error.postMODEtheta.NEW1.simulations	=	rbind(L2error.postMODEtheta.NEW1.simulations, 	(postMODEtheta.NEW1 	- truth)^2	)
		#
		L1error.MLEtheta.NEW1.simulations		=	rbind(L1error.MLEtheta.NEW1.simulations, 		abs(MLEtheta.NEW1		- truth)	)
		L2error.MLEtheta.NEW1.simulations		=	rbind(L2error.MLEtheta.NEW1.simulations, 		(MLEtheta.NEW1		- truth)^2	)
		#
		L1error.postMODEtheta.NEW2.simulations	=	rbind(L1error.postMODEtheta.NEW2.simulations, 	abs(postMODEtheta.NEW2- truth)	)
		L2error.postMODEtheta.NEW2.simulations	=	rbind(L2error.postMODEtheta.NEW2.simulations, 	(postMODEtheta.NEW2 	- truth)^2	)
		#
		L1error.MLEtheta.NEW2.simulations		=	rbind(L1error.MLEtheta.NEW2.simulations, 		abs(MLEtheta.NEW2		- truth)	)
		L2error.MLEtheta.NEW2.simulations		=	rbind(L2error.MLEtheta.NEW2.simulations, 		(MLEtheta.NEW2		- truth)^2	)
		#
		if (fileNum == length(RDataFiles))	{
			if (simStudy != 10)		{	# If not results for twisted normal model.
				postMODEtheta.NEW1.simulations.mean0=	matrix(apply(postMODEtheta.NEW1.simulations, 2, mean), nrow = 1)
				postMODEtheta.NEW1.simulations.sd0	=	matrix(apply(postMODEtheta.NEW1.simulations, 2, sd), nrow = 1)
				MLEtheta.NEW1.simulations.mean0	=	matrix(apply(MLEtheta.NEW1.simulations, 2, mean), nrow = 1)
				MLEtheta.NEW1.simulations.sd0		=	matrix(apply(MLEtheta.NEW1.simulations, 2, sd), nrow = 1)

				postMODEtheta.NEW2.simulations.mean0=	matrix(apply(postMODEtheta.NEW2.simulations, 2, mean), nrow = 1)
				postMODEtheta.NEW2.simulations.sd0	=	matrix(apply(postMODEtheta.NEW2.simulations, 2, sd), nrow = 1)
				MLEtheta.NEW2.simulations.mean0	=	matrix(apply(MLEtheta.NEW2.simulations, 2, mean), nrow = 1)
				MLEtheta.NEW2.simulations.sd0		=	matrix(apply(MLEtheta.NEW2.simulations, 2, sd), nrow = 1)

				post.df.NEW.simulations.mean0		=	matrix(mean(post.df.NEW.simulations), nrow = 1)
				post.df.NEW.simulations.sd0		=	matrix(sd(post.df.NEW.simulations),   nrow = 1)
				post.scale.NEW.simulations.mean0	=	matrix(apply(post.scale.NEW.simulations, 2, mean), nrow = 1)
				post.scale.NEW.simulations.sd0		=	matrix(apply(post.scale.NEW.simulations, 2, sd), nrow = 1)
				#
				MAE.postMODEtheta.NEW1.simulations0=	matrix(apply(L1error.postMODEtheta.NEW1.simulations, 2, mean), nrow = 1) # Mean Absolute Error (MAE)
				MSE.postMODEtheta.NEW1.simulations0=	matrix(apply(L2error.postMODEtheta.NEW1.simulations, 2, mean), nrow = 1) # Mean Squared Error (MSE)
				MAE.MLEtheta.NEW1.simulations0		=	matrix(apply(L1error.MLEtheta.NEW1.simulations, 2, mean), nrow = 1) 		# Mean Absolute Error (MAE)
				MSE.MLEtheta.NEW1.simulations0		=	matrix(apply(L2error.MLEtheta.NEW1.simulations, 2, mean), nrow = 1) 		# Mean Squared Error (MSE)
				#
				MAE.postMODEtheta.NEW2.simulations0=	matrix(apply(L1error.postMODEtheta.NEW2.simulations, 2, mean), nrow = 1) # Mean Absolute Error (MAE)
				MSE.postMODEtheta.NEW2.simulations0=	matrix(apply(L2error.postMODEtheta.NEW2.simulations, 2, mean), nrow = 1) # Mean Squared Error (MSE)
				MAE.MLEtheta.NEW2.simulations0		=	matrix(apply(L1error.MLEtheta.NEW2.simulations, 2, mean), nrow = 1) 		# Mean Absolute Error (MAE)
				MSE.MLEtheta.NEW2.simulations0		=	matrix(apply(L2error.MLEtheta.NEW2.simulations, 2, mean), nrow = 1) 		# Mean Squared Error (MSE)
				#
				postMODEtheta.OLD.simulations.mean0=	matrix(apply(postMODEtheta.OLD.simulations, 2, 	mean, 	na.rm = T), nrow = 1)
				postMODEtheta.OLD.simulations.sd0	=	matrix(apply(postMODEtheta.OLD.simulations, 2, 	sd, 	na.rm = T), nrow = 1)
				MLEtheta.OLD.simulations.mean0		=	matrix(apply(MLEtheta.OLD.simulations, 		2, 	mean, 	na.rm = T), nrow = 1)
				MLEtheta.OLD.simulations.sd0		=	matrix(apply(MLEtheta.OLD.simulations, 		2, 	sd, 	na.rm = T), nrow = 1)
				post.df.OLD.simulations.mean0		=	matrix(mean(post.df.OLD.simulations), nrow = 1)
				post.df.OLD.simulations.sd0		=	matrix(sd(post.df.OLD.simulations),   nrow = 1)
				post.scale.OLD.simulations.mean0	=	matrix(apply(post.scale.OLD.simulations, 2, mean), nrow = 1)
				post.scale.OLD.simulations.sd0		=	matrix(apply(post.scale.OLD.simulations, 2, sd), nrow = 1)
				MAE.postMODEtheta.OLD.simulations0	=	matrix(apply(L1error.postMODEtheta.OLD.simulations, 2, mean, 	na.rm = T), nrow = 1) # Mean Absolute Error (MAE)
				MSE.postMODEtheta.OLD.simulations0	=	matrix(apply(L2error.postMODEtheta.OLD.simulations, 2, mean, 	na.rm = T), nrow = 1) # Mean Squared Error (MSE)
				MAE.MLEtheta.OLD.simulations0		=	matrix(apply(L1error.MLEtheta.OLD.simulations, 2, mean, 	na.rm = T), nrow = 1) 	# Mean Absolute Error (MAE)
				MSE.MLEtheta.OLD.simulations0		=	matrix(apply(L2error.MLEtheta.OLD.simulations, 2, mean, 	na.rm = T), nrow = 1) 	# Mean Squared Error (MSE)
			}
		}
	print(paste("Completed file #", fileNum, " of ", length(RDataFiles), " files, in simulation study #", simStudy, ".", sep='')); flush.console()
	}
	if (simStudy != 10)		{	# If not results for twisted normal model.
		postMODEtheta.NEW1.simulations.mean[simStudy,1:length(truth)]	=	round(postMODEtheta.NEW1.simulations.mean0, 2)
		postMODEtheta.NEW1.simulations.sd[simStudy, 1:length(truth)]	=	round(postMODEtheta.NEW1.simulations.sd0, 2)
		MLEtheta.NEW1.simulations.mean[simStudy, 1:length(truth)]	=	round(MLEtheta.NEW1.simulations.mean0, 2)
		MLEtheta.NEW1.simulations.sd[simStudy, 1:length(truth)]		=	round(MLEtheta.NEW1.simulations.sd0, 2)
		#
		postMODEtheta.NEW2.simulations.mean[simStudy,1:length(truth)]	=	round(postMODEtheta.NEW2.simulations.mean0, 2)
		postMODEtheta.NEW2.simulations.sd[simStudy, 1:length(truth)]	=	round(postMODEtheta.NEW2.simulations.sd0, 2)
		MLEtheta.NEW2.simulations.mean[simStudy, 1:length(truth)]	=	round(MLEtheta.NEW2.simulations.mean0, 2)
		MLEtheta.NEW2.simulations.sd[simStudy, 1:length(truth)]		=	round(MLEtheta.NEW2.simulations.sd0, 2)
		#
		MAE.postMODEtheta.NEW1.simulations[simStudy, 1:length(truth)]	= 	round(MAE.postMODEtheta.NEW1.simulations0, 2)
		MSE.postMODEtheta.NEW1.simulations[simStudy, 1:length(truth)]	= 	round(MSE.postMODEtheta.NEW1.simulations0, 2)
		MAE.MLEtheta.NEW1.simulations[simStudy, 1:length(truth)]		= 	round(MAE.MLEtheta.NEW1.simulations0, 2)
		MSE.MLEtheta.NEW1.simulations[simStudy, 1:length(truth)]		= 	round(MSE.MLEtheta.NEW1.simulations0, 2)	
		#
		MAE.postMODEtheta.NEW2.simulations[simStudy, 1:length(truth)]	= 	round(MAE.postMODEtheta.NEW2.simulations0, 2)
		MSE.postMODEtheta.NEW2.simulations[simStudy, 1:length(truth)]	= 	round(MSE.postMODEtheta.NEW2.simulations0, 2)
		MAE.MLEtheta.NEW2.simulations[simStudy, 1:length(truth)]		= 	round(MAE.MLEtheta.NEW2.simulations0, 2)
		MSE.MLEtheta.NEW2.simulations[simStudy, 1:length(truth)]		= 	round(MSE.MLEtheta.NEW2.simulations0, 2)	
		#
		post.df.NEW.simulations.mean[simStudy, 1]					=	round(post.df.NEW.simulations.mean0, 	2)
		post.df.NEW.simulations.sd[simStudy, 1]						=	round(post.df.NEW.simulations.sd0, 	2)
		post.scale.NEW.simulations.mean[simStudy, 1:length(c(post.scale.NEW))]=	round(post.scale.NEW.simulations.mean0, 2)
		post.scale.NEW.simulations.sd[simStudy,   1:length(c(post.scale.NEW))]=	round(post.scale.NEW.simulations.sd0, 2)
		#
		postMODEtheta.OLD.simulations.mean[simStudy,1:length(truth)]	=	round(postMODEtheta.OLD.simulations.mean0, 2)
		postMODEtheta.OLD.simulations.sd[simStudy, 1:length(truth)]	=	round(postMODEtheta.OLD.simulations.sd0, 2)
		MLEtheta.OLD.simulations.mean[simStudy, 1:length(truth)]		=	round(MLEtheta.OLD.simulations.mean0, 2)
		MLEtheta.OLD.simulations.sd[simStudy, 1:length(truth)]		=	round(MLEtheta.OLD.simulations.sd0, 2)
		#
		MAE.postMODEtheta.OLD.simulations[simStudy, 1:length(truth)]	= 	round(MAE.postMODEtheta.OLD.simulations0, 2)
		MSE.postMODEtheta.OLD.simulations[simStudy, 1:length(truth)]	= 	round(MSE.postMODEtheta.OLD.simulations0, 2)
		MAE.MLEtheta.OLD.simulations[simStudy, 1:length(truth)]		= 	round(MAE.MLEtheta.OLD.simulations0, 2)
		MSE.MLEtheta.OLD.simulations[simStudy, 1:length(truth)]		= 	round(MSE.MLEtheta.OLD.simulations0, 2)	
		#
		post.df.OLD.simulations.mean[simStudy, 1]					=	round(post.df.OLD.simulations.mean0, 2)
		post.df.OLD.simulations.sd[simStudy, 1]						=	round(post.df.OLD.simulations.sd0, 2)
		post.scale.OLD.simulations.mean[simStudy, 1:length(c(post.scale.NEW))]=	round(post.scale.OLD.simulations.mean0, 2)
		post.scale.OLD.simulations.sd[simStudy,   1:length(c(post.scale.NEW))]=	round(post.scale.OLD.simulations.sd0, 2)
	}
	if (simStudy == 10)		{	# If results for twisted normal model.
		postMODEtheta.NEW1.twistedNormal			=	round(cbind(matrix(postMODEtheta.NEW1[1:2],nrow = 1), min(postMODEtheta.NEW1[3:300]), max(postMODEtheta.NEW1[3:300])),2)
		rownames(postMODEtheta.NEW1.twistedNormal)	=	'postMODE.NEW1.twistedNormal'
		colnames(postMODEtheta.NEW1.twistedNormal)	=	c(paste('th', 1 : 2, sep = ''), 'min(theta3...300)', 'max(theta3...300)')
		MLEtheta.NEW1.twistedNormal				=	round(cbind(matrix(MLEtheta.NEW1[1:2],nrow = 1), min(MLEtheta.NEW1[3:300]), max(MLEtheta.NEW1[3:300])),2)
		rownames(MLEtheta.NEW1.twistedNormal)		=	'MLE.NEW1.twistedNormal'
		colnames(MLEtheta.NEW1.twistedNormal)		=	c(paste('th', 1 : 2, sep = ''), 'min(theta3...300)', 'max(theta3...300)')
		#
		postMODEtheta.NEW2.twistedNormal			=	round(cbind(matrix(postMODEtheta.NEW2[1:2],nrow = 1), min(postMODEtheta.NEW2[3:300]), max(postMODEtheta.NEW2[3:300])),2)
		rownames(postMODEtheta.NEW2.twistedNormal)	=	'postMODE.NEW2.twistedNormal'
		colnames(postMODEtheta.NEW2.twistedNormal)	=	c(paste('th', 1 : 2, sep = ''), 'min(theta3...300)', 'max(theta3...300)')
		MLEtheta.NEW2.twistedNormal				=	round(cbind(matrix(MLEtheta.NEW2[1:2],nrow = 1), min(MLEtheta.NEW2[3:300]), max(MLEtheta.NEW2[3:300])),2)
		rownames(MLEtheta.NEW2.twistedNormal)		=	'MLE.NEW2.twistedNormal'
		colnames(MLEtheta.NEW2.twistedNormal)		=	c(paste('th', 1 : 2, sep = ''), 'min(theta3...300)', 'max(theta3...300)')
		#
		post.df.NEW.twistedNormal					=	round(post.df.NEW.simulations, 2)
		post.scale.NEW3300							=	post.scale.NEW[3:300, 3:300]
		post.scale.NEW3300							=	post.scale.NEW3300[row(post.scale.NEW3300) != col(post.scale.NEW3300)]
		post.scale.NEW.twistedNormal				=	rbind(	cbind(post.scale.NEW[1,2], min(post.scale.NEW[1,3:300]), max(post.scale.NEW[1,3:300]) ),
														  	cbind(				 NA, min(post.scale.NEW[2,3:300]), max(post.scale.NEW[2,3:300]) ),
														  	cbind(				 NA, min(post.scale.NEW3300     ), max(post.scale.NEW3300     ) ))
		colnames(post.scale.NEW.twistedNormal)		=	c('theta2', 'min(theta3...300)', 'max(theta3...300)')
		rownames(post.scale.NEW.twistedNormal)		=	c('theta1', 'theta2'           , 'theta3...300'     )
		#
		postMODEtheta.OLD.twistedNormal			=	round(cbind(matrix(postMODEtheta.OLD[1:2],nrow = 1), min(postMODEtheta.OLD[3:300]), max(postMODEtheta.OLD[3:300])),2)
		rownames(postMODEtheta.OLD.twistedNormal)	=	'postMODE.OLD.twistedNormal'
		colnames(postMODEtheta.OLD.twistedNormal)	=	c(paste('th', 1 : 2, sep = ''), 'min(theta3...300)', 'max(theta3...300)')
		MLEtheta.OLD.twistedNormal					=	round(cbind(matrix(MLEtheta.OLD[1:2],nrow = 1), min(MLEtheta.OLD[3:300]), max(MLEtheta.OLD[3:300])),2)
		rownames(MLEtheta.OLD.twistedNormal)		=	'MLE.OLD.twistedNormal'
		colnames(MLEtheta.OLD.twistedNormal)		=	c(paste('th', 1 : 2, sep = ''), 'min(theta3...300)', 'max(theta3...300)')
		post.df.OLD.twistedNormal					=	round(post.df.OLD.simulations, 2)
		post.scale.OLD								=	post.Correlations
		post.scale.OLD3300							=	post.Correlations[3:300, 3:300]
		post.scale.OLD3300							=	post.scale.OLD3300[row(post.scale.OLD3300) != col(post.scale.OLD3300)]
		post.scale.OLD.twistedNormal				=	rbind(	cbind(post.scale.OLD[1,2], min(post.scale.OLD[1,3:300]), max(post.scale.OLD[1,3:300]) ),
														  	cbind(				 NA, min(post.scale.OLD[2,3:300]), max(post.scale.OLD[2,3:300]) ),
														  	cbind(				 NA, min(post.scale.OLD3300     ), max(post.scale.OLD3300     ) ))
		colnames(post.scale.OLD.twistedNormal)		=	c('theta2', 'min(theta3...300)', 'max(theta3...300)')
		rownames(post.scale.OLD.twistedNormal)		=	c('theta1', 'theta2'           , 'theta3...300'     )
	}

	if (simStudy == 62)	{ 
		postMODEtheta.NEW1.simulations.meansd			=	matrix(paste(postMODEtheta.NEW1.simulations.mean, ' (', postMODEtheta.NEW1.simulations.sd, ')', sep=''), nrow = 62, ncol = 18)
		rownames(postMODEtheta.NEW1.simulations.meansd)	=	simConditionNames
		colnames(postMODEtheta.NEW1.simulations.meansd)	=	paste('th', 1 : 18, sep = '')
		#
		MLEtheta.NEW1.simulations.meansd				=	matrix(paste(MLEtheta.NEW1.simulations.mean, ' (', MLEtheta.NEW1.simulations.sd, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MLEtheta.NEW1.simulations.meansd)		=	simConditionNames
		colnames(MLEtheta.NEW1.simulations.meansd)		=	paste('th', 1 : 18, sep = '')
		#
		MAE.MSE.postMODEtheta.NEW1.simulations			=	matrix(paste(MAE.postMODEtheta.NEW1.simulations, ' (', MSE.postMODEtheta.NEW1.simulations, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MAE.MSE.postMODEtheta.NEW1.simulations)=	simConditionNames
		colnames(MAE.MSE.postMODEtheta.NEW1.simulations)=	paste('th', 1 : 18, sep = '')
		#
		MAE.MSE.MLEtheta.NEW1.simulations				=	matrix(paste(MAE.MLEtheta.NEW1.simulations, ' (', MSE.MLEtheta.NEW1.simulations, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MAE.MSE.MLEtheta.NEW1.simulations)		=	simConditionNames
		colnames(MAE.MSE.MLEtheta.NEW1.simulations)		=	paste('th', 1 : 18, sep = '')
		#
		postMODEtheta.NEW2.simulations.meansd			=	matrix(paste(postMODEtheta.NEW2.simulations.mean, ' (', postMODEtheta.NEW2.simulations.sd, ')', sep=''), nrow = 62, ncol = 18)
		rownames(postMODEtheta.NEW2.simulations.meansd)	=	simConditionNames
		colnames(postMODEtheta.NEW2.simulations.meansd)	=	paste('th', 1 : 18, sep = '')
		#
		MLEtheta.NEW2.simulations.meansd				=	matrix(paste(MLEtheta.NEW2.simulations.mean, ' (', MLEtheta.NEW2.simulations.sd, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MLEtheta.NEW2.simulations.meansd)		=	simConditionNames
		colnames(MLEtheta.NEW2.simulations.meansd)		=	paste('th', 1 : 18, sep = '')
		#
		MAE.MSE.postMODEtheta.NEW2.simulations			=	matrix(paste(MAE.postMODEtheta.NEW2.simulations, ' (', MSE.postMODEtheta.NEW2.simulations, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MAE.MSE.postMODEtheta.NEW2.simulations)=	simConditionNames
		colnames(MAE.MSE.postMODEtheta.NEW2.simulations)=	paste('th', 1 : 18, sep = '')
		#
		MAE.MSE.MLEtheta.NEW2.simulations				=	matrix(paste(MAE.MLEtheta.NEW2.simulations, ' (', MSE.MLEtheta.NEW2.simulations, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MAE.MSE.MLEtheta.NEW2.simulations)		=	simConditionNames
		colnames(MAE.MSE.MLEtheta.NEW2.simulations)		=	paste('th', 1 : 18, sep = '')
		#
		post.scale.NEW.simulations.meansd 				= 	matrix(paste(post.scale.NEW.simulations.mean, ' (', post.scale.NEW.simulations.sd, ')',sep=''), nrow = 62)
		rownames(post.scale.NEW.simulations.meansd)		=	simConditionNames
		post.df.NEW.simulations.meansd 					= 	matrix(paste(post.df.NEW.simulations.mean, ' (', post.df.NEW.simulations.sd, ')',sep=''), nrow = 62)
		rownames(post.df.NEW.simulations.meansd)		=	simConditionNames
		#
		postMODEtheta.OLD.simulations.meansd			=	matrix(paste(postMODEtheta.OLD.simulations.mean, ' (', postMODEtheta.OLD.simulations.sd, ')', sep=''), nrow = 62, ncol = 18)
		rownames(postMODEtheta.OLD.simulations.meansd)	=	simConditionNames
		colnames(postMODEtheta.OLD.simulations.meansd)	=	paste('th', 1 : 18, sep = '')
		#
		MLEtheta.OLD.simulations.meansd				=	matrix(paste(MLEtheta.OLD.simulations.mean, ' (', MLEtheta.OLD.simulations.sd, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MLEtheta.OLD.simulations.meansd)		=	simConditionNames
		colnames(MLEtheta.OLD.simulations.meansd)		=	paste('th', 1 : 18, sep = '')
		#
		MAE.MSE.postMODEtheta.OLD.simulations			=	matrix(paste(MAE.postMODEtheta.OLD.simulations, ' (', MSE.postMODEtheta.OLD.simulations, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MAE.MSE.postMODEtheta.OLD.simulations)	=	simConditionNames
		colnames(MAE.MSE.postMODEtheta.OLD.simulations)	=	paste('th', 1 : 18, sep = '')
		#
		MAE.MSE.MLEtheta.OLD.simulations				=	matrix(paste(MAE.MLEtheta.OLD.simulations, ' (', MSE.MLEtheta.OLD.simulations, ')', sep=''), nrow = 62, ncol = 18)
		rownames(MAE.MSE.MLEtheta.OLD.simulations)		=	simConditionNames
		colnames(MAE.MSE.MLEtheta.OLD.simulations)		=	paste('th', 1 : 18, sep = '')
		#
		post.scale.OLD.simulations.meansd 				= 	matrix(paste(post.scale.OLD.simulations.mean, ' (', post.scale.OLD.simulations.sd, ')',sep=''), nrow = 62)
		rownames(post.scale.OLD.simulations.meansd)		=	simConditionNames
		post.df.OLD.simulations.meansd 				= 	matrix(paste(post.df.OLD.simulations.mean, ' (', post.df.OLD.simulations.sd, ')',sep=''), nrow = 62)
		rownames(post.df.OLD.simulations.meansd)		=	simConditionNames
		#
	}

#OUT = rbind(	MAE.postMODEtheta.NEW1.simulations0,
#			MAE.postMODEtheta.NEW2.simulations0,
#			MAE.postMODEtheta.OLD.simulations0,
#			MSE.postMODEtheta.NEW1.simulations0,
#			MSE.postMODEtheta.NEW2.simulations0,
#			MSE.postMODEtheta.OLD.simulations0,
#			MAE.MLEtheta.NEW1.simulations0,
#			MAE.MLEtheta.NEW2.simulations0,
#			MAE.MLEtheta.OLD.simulations0,
#			MSE.MLEtheta.NEW1.simulations0,
#			MSE.MLEtheta.NEW2.simulations0,
#			MSE.MLEtheta.OLD.simulations0)	
#	
#rownames(OUT)		=	c(	'MAE.postMODEtheta.NEW1.simulations0',
#						'MAE.postMODEtheta.NEW2.simulations0',
#						'MAE.postMODEtheta.OLD.simulations0',
#						'MSE.postMODEtheta.NEW1.simulations0',
#						'MSE.postMODEtheta.NEW2.simulations0',
#						'MSE.postMODEtheta.OLD.simulations0',
#						'MAE.MLEtheta.NEW1.simulations0',
#						'MAE.MLEtheta.NEW2.simulations0',
#						'MAE.MLEtheta.OLD.simulations0',
#						'MSE.MLEtheta.NEW1.simulations0',
#						'MSE.MLEtheta.NEW2.simulations0',
#						'MSE.MLEtheta.OLD.simulations0')
#colnames(OUT)		=	paste('th', 1 : length(truth), sep='')
#round(OUT, 2)
}




#
# OLD vs. NEW results:


for (i in 1:62){
	if (i == 1) {	postMODEtheta.NEWvsOLD.simulations.meansd 	= 	c()
					MAE.MSE.postMODEtheta.NEWvsOLD.simulations	=	c()
					MLEtheta.NEWvsOLD.simulations.meansd		=	c()
					MAE.MSE.MLEtheta.NEWvsOLD.simulations		=	c()
					post.df.NEWvsOLD.simulations.meansd			=	c()
					post.scale.NEWvsOLD.simulations.meansd		=	c()		}
		#
		postMODEtheta.NEWvsOLD.simulations.meansd 	= 	rbind(postMODEtheta.NEWvsOLD.simulations.meansd	, 	postMODEtheta.NEW1.simulations.meansd[i,] 	)
		postMODEtheta.NEWvsOLD.simulations.meansd 	= 	rbind(postMODEtheta.NEWvsOLD.simulations.meansd	, 	postMODEtheta.NEW2.simulations.meansd[i,] 	)
		postMODEtheta.NEWvsOLD.simulations.meansd 	= 	rbind(postMODEtheta.NEWvsOLD.simulations.meansd	, 	postMODEtheta.OLD.simulations.meansd[i,] 	)
		#
		MAE.MSE.postMODEtheta.NEWvsOLD.simulations	=	rbind(MAE.MSE.postMODEtheta.NEWvsOLD.simulations,	MAE.MSE.postMODEtheta.NEW1.simulations[i,]	)
		MAE.MSE.postMODEtheta.NEWvsOLD.simulations	=	rbind(MAE.MSE.postMODEtheta.NEWvsOLD.simulations,	MAE.MSE.postMODEtheta.NEW2.simulations[i,]	)
		MAE.MSE.postMODEtheta.NEWvsOLD.simulations	=	rbind(MAE.MSE.postMODEtheta.NEWvsOLD.simulations,	MAE.MSE.postMODEtheta.OLD.simulations[i,]	)
		#
		MLEtheta.NEWvsOLD.simulations.meansd 		= 	rbind(MLEtheta.NEWvsOLD.simulations.meansd		, 	MLEtheta.NEW1.simulations.meansd[i,] 		)
		MLEtheta.NEWvsOLD.simulations.meansd 		= 	rbind(MLEtheta.NEWvsOLD.simulations.meansd		, 	MLEtheta.NEW2.simulations.meansd[i,] 		)
		MLEtheta.NEWvsOLD.simulations.meansd 		= 	rbind(MLEtheta.NEWvsOLD.simulations.meansd		, 	MLEtheta.OLD.simulations.meansd[i,] 		)
		#
		MAE.MSE.MLEtheta.NEWvsOLD.simulations		=	rbind(MAE.MSE.MLEtheta.NEWvsOLD.simulations		,	MAE.MSE.MLEtheta.NEW1.simulations[i,]		)
		MAE.MSE.MLEtheta.NEWvsOLD.simulations		=	rbind(MAE.MSE.MLEtheta.NEWvsOLD.simulations		,	MAE.MSE.MLEtheta.NEW2.simulations[i,]		)
		MAE.MSE.MLEtheta.NEWvsOLD.simulations		=	rbind(MAE.MSE.MLEtheta.NEWvsOLD.simulations		,	MAE.MSE.MLEtheta.OLD.simulations[i,]		)
		#
		post.df.NEWvsOLD.simulations.meansd		=	rbind(post.df.NEWvsOLD.simulations.meansd		, 	post.df.NEW.simulations.meansd[i,]			)
		post.df.NEWvsOLD.simulations.meansd		=	rbind(post.df.NEWvsOLD.simulations.meansd		, 	post.df.OLD.simulations.meansd[i,]			)
		post.scale.NEWvsOLD.simulations.meansd		=	rbind(post.scale.NEWvsOLD.simulations.meansd	, 	post.scale.NEW.simulations.meansd[i,]		)
		post.scale.NEWvsOLD.simulations.meansd		=	rbind(post.scale.NEWvsOLD.simulations.meansd	, 	post.scale.OLD.simulations.meansd[i,]		)

	if (i == 62){
		colnames(postMODEtheta.NEWvsOLD.simulations.meansd)	=	colnames(postMODEtheta.NEW1.simulations.meansd)
		colnames(MAE.MSE.postMODEtheta.NEWvsOLD.simulations)	=	colnames(postMODEtheta.NEW1.simulations.meansd)
		colnames(MLEtheta.NEWvsOLD.simulations.meansd)		=	colnames(postMODEtheta.NEW1.simulations.meansd)
		colnames(MAE.MSE.MLEtheta.NEWvsOLD.simulations)		=	colnames(postMODEtheta.NEW1.simulations.meansd)
		colnames(post.df.NEWvsOLD.simulations.meansd)		=	"d.f."
		newRowNames	=   paste( rep(c("NEW: ","OLD: "), 62), 	simConditionNames[sort(rep(1:62,2))], sep = '')
		rownames(post.df.NEWvsOLD.simulations.meansd)		=	newRowNames
		rownames(post.scale.NEWvsOLD.simulations.meansd)	=	newRowNames
		#
		newRowNames2	= 	paste( rep(c("NEW1: ", "NEW2: ", "OLD: "), 62), simConditionNames[sort(rep(1:62,3))], sep = '')
		rownames(postMODEtheta.NEWvsOLD.simulations.meansd)	=	newRowNames2
		rownames(MAE.MSE.postMODEtheta.NEWvsOLD.simulations)=	newRowNames2
		rownames(MLEtheta.NEWvsOLD.simulations.meansd)		=	newRowNames2
		rownames(MAE.MSE.MLEtheta.NEWvsOLD.simulations)		=	newRowNames2
	}
}

noquote(postMODEtheta.NEWvsOLD.simulations.meansd)
noquote(MAE.MSE.postMODEtheta.NEWvsOLD.simulations)
noquote(MLEtheta.NEWvsOLD.simulations.meansd)
noquote(MAE.MSE.MLEtheta.NEWvsOLD.simulations)
noquote(post.df.NEWvsOLD.simulations.meansd)
noquote(post.scale.NEWvsOLD.simulations.meansd)

#
# OLD vs. NEW results for twisted normal model:
rbind(	postMODEtheta.NEW1.twistedNormal,
		postMODEtheta.NEW2.twistedNormal, 	
		postMODEtheta.OLD.twistedNormal,
		MLEtheta.NEW1.twistedNormal,
		MLEtheta.NEW2.twistedNormal,
 		MLEtheta.OLD.twistedNormal)
rbind(	post.df.NEW.twistedNormal, 
		post.df.OLD.twistedNormal)
rbind(	post.scale.NEW.twistedNormal, 
		post.scale.OLD.twistedNormal)


save.image("C:\\Users\\George Karabatsos\\Desktop\\George\\PAPERS\\copulaABCdrf\\code ABCnetwork\\Results\\Simulation Results\\Comparing Mode MLE estimators (old vs new) 7-25-24.RData") 