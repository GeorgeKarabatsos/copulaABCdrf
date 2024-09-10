
par(mfrow=c(1,5))
matplot(Samples.theta[,1], type = 'l',xlab='MCMC iteration','ylab' = 'theta_1')
matplot(Samples.theta[,2], type = 'l',xlab='MCMC iteration','ylab' = 'theta_2')
matplot(Samples.theta[,3], type = 'l',xlab='MCMC iteration','ylab' = 'theta_3')
matplot(Samples.theta[,4], type = 'l',xlab='MCMC iteration','ylab' = 'theta_4')
matplot(Samples.theta[,5], type = 'l',xlab='MCMC iteration','ylab' = 'theta_5')