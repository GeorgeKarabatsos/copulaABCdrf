# This code file compares Price model and ERGM model using DRF.
rm(list=ls())
# Import Price model analysis output:
load("C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork\\Results\\Real network data\\cit-HepPh\\Price cit-HepPh 2024-01-28 15_58_51.RData")
Price.MPLE.Table	=	MPLE.table
Y.Price	=	matrix(1, nrow=nrow(Price.MPLE.Table), ncol = 1) # model index = 1 for Price model.

# Import ERGM analysis output:
load("C:\\Users\\George Karabatsos\\Desktop\\code ABCnetwork\\Results\\Real network data\\cit-HepPh\\ERGM cit-HepPh 2024-01-28 17_11_04.RData")
ERGM.MPLE.Table	=	MPLE.table
Y.ERGM		=	matrix(0, nrow=nrow(ERGM.MPLE.Table), ncol = 1) # model index = 0 for ERGM.

# Combine the data:
Y			=	rbind(Y.Price, Y.ERGM)
X			=	rbind(Price.MPLE.Table, ERGM.MPLE.Table)

# Geometric median DAC estimate of MPLE:
# MPLEx$p
#> MPLEx$p	
#$p
#gwidegree(cutoff = 28093)           gwidegree.decay                 triangles 
#             -0.8418461                 0.2365466                 1.6904300

MPLEx				=	c(-0.8418461, 0.2365466, 1.6904300)
MPLEterms			=	c('gwidegree', 'gwidegree.decay', 'triangle')
names(MPLEx)		=	MPLEterms


# Train the DRF:
install.packages('drf');	library('drf')
drf.forest 		=	drf(X = X, Y = Y) 
X.test				=	MPLEx
# Compute posterior probability of Price model
Pr.Price			=	c(predict(drf.forest, newdata	= X.test, functional = "mean")$mean)
# Posterior probability of ERGM equals:
Pr.ERGM			=	1 - Pr.Price

cbind(Pr.Price, Pr.ERGM)

# Comparing model posterior probabilities:
#> cbind(Pr.Price, Pr.ERGM)
#      Pr.Price      Pr.ERGM
#[1,] 0.9992924 0.0007075758