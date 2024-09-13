########################################################################################
###Section 1
####Simulate 5 species MLGCP using Multilogreg package
remove(list = ls())
library(Multilogreg)
library(spatstat)
library(RandomFields)
RFoptions(install="no")
library(ggplot2)

#reproducibility
set.seed(1)

#Initialize the MLGCP
# Size of the observation window
n.x <- n.y <- 1
xx=seq(0,n.x,length=100)
yy=seq(0,n.x,length=100)

# Simulating a covariate
set.seed(1)
cov <- as.matrix(RFsimulate(RMexp(var=1,scale=0.05), x=xx, y=yy, grid=TRUE))

# Simulating the background intensity
gamma <- 0.5
background <- as.matrix(RFsimulate(RMgauss(var=1,scale=0.2), x=xx, y=yy, grid=TRUE))*gamma-gamma^2/2

#Set up parameters
beta1 <- c(0.1,0.2,0.3,0.4,0.5)
beta2 <- c(-0.1,-0.2,0,0.1,0.2)
beta2 <- as.matrix(beta2)

# Parameters in the MLGCP
alpha <- matrix(c(0.5,-1,0.5,0,-1,0,0,0.5,0,0.5),nrow=5,byrow=TRUE)
xi    <- c(0.02,0.03)
sigma <- matrix(c(sqrt(0.5),sqrt(0.5),sqrt(0.5),sqrt(0.5),sqrt(0.5)),ncol=1)
phi   <- matrix(c(0.02,0.02,0.03,0.03,0.04),ncol=1)

n.window <- n.x
#Number of points per species
n.points <- c(rep(200,5))

#Simulation of a multivariate LGCP
nsim <-10 # Change to 100 for full version
X <- list()

for (i in 1:nsim){
  print(i)
  X[[i]] <- sim_lgcp_multi(basecov=background,covariate=cov,betas=beta2,alphas=alpha,xi=xi,
                           sigma=sigma,phis=phi, n.window=n.window,n.points=n.points,beta0s=beta1)
  
}

#####################################################################################
###Section 2
####Fitting MLGCP model for the simulated data
nspecies  <- dim(beta2)[1]

#Second Order Analysis
# Globals
len.x=X[[1]]$markedprocess$window$xrange[2]-X[[1]]$markedprocess$window$xrange[1]
len.y=X[[1]]$markedprocess$window$yrange[2]-X[[1]]$markedprocess$window$yrange[1]
lw=(len.x+len.y)/2
latent = 2
Rmax   = 1
rl     = 1000
r      = seq(0,Rmax,length=rl)

# Parameter spaces (upper,lower)
xib = c(lw/20,lw/50)
sib = c(0.75,0.25)
phb = c(lw/20,lw/100)

par.est <- list()
Betahat <- list()
cv <- list()
C <- as.im(t(cov),W=owin())
F <- as.function(as.im(t(cov),W=owin()))

#Initialize the parameters
beta1 <- round(runif(5,0.1,0.5),1)
beta2.1 <- round(runif(2,-0.2,0),1)
beta2.2 <- round(runif(3,0,0.2),1)
beta2 <- c(beta2.1,beta2.2)
beta2 <- as.matrix(beta2)
Beta.init <- cbind(beta1,beta2)

# Parameter estimation
for ( i in 1:nsim){
  print(i)
  cv <- as.matrix(cbind(1,F(X[[i]]$markedprocess$x,X[[i]]$markedprocess$y)),ncol=2)
  Betahat[[i]]   <- FirstOrderCCL(X=X[[i]]$markedprocess,Beta0=Beta.init,covariate = C)
  par.est[[i]] <- PenalizedSecondOrderCCL(X=X[[i]]$markedprocess,covariate = cv,Beta = Betahat[[i]]$betahat,Rmax = Rmax,
                                          xibound=xib,sigmabound = sib,phibound=phb, lat = latent,tol=2e-5,lamb = 0)
}

# Beta estimates
beta0.est <- matrix(NA, nrow=nsim, ncol=6)
beta1.est <- matrix(NA, nrow=nsim, ncol=6)

for ( i in 1:nsim){
  beta0.est[i,] <- Betahat[[i]]$betahat[1,]
  beta1.est[i,] <- Betahat[[i]]$betahat[2,]
}

beta0_est <- apply(beta0.est, 2, mean)
beta1_est <- apply(beta1.est, 2, mean)

# Parameter estimates
mlxi.est <- matrix(NA, nrow=nsim, ncol=2)
mlalpha.est <- array(NA, dim=c(5, 2, nsim))
mlsig.est <- matrix(NA, nrow=nsim, ncol=5)
mlphi.est <- matrix(NA, nrow=nsim, ncol=5)

for(i in 1:nsim){
  mlxi.est[i,]  <- par.est[[i]][[1]]$xi
  mlalpha.est[, , i] <- par.est[[i]][[1]]$alpha[,,1]
  mlsig.est[i,] <- par.est[[i]][[1]]$sigma
  mlphi.est[i,] <- par.est[[i]][[1]]$phi
}

#Average parameter estimates
xi_est <- round(apply(mlxi.est, 2, mean),3)

sig_est <- matrix(round(apply(mlsig.est, 2, mean),3),ncol=1)

phi_est <- matrix(round(apply(mlphi.est, 2, mean),3),ncol=1)

alp_est <- matrix(round(apply(mlalpha.est, c(1, 2), mean),3),ncol=2)

##################################################################################
###Section 3
####Re-sample from the fitted model
sampleMLGCP <- list()
for (i in 1:nsim){
  print(i)
  sampleMLGCP[[i]] <- sim_lgcp_multi(basecov=background,covariate=cov,betas=beta2,alphas=mlalpha.est[,,i],xi=mlxi.est[i,],
                                     sigma=matrix(mlsig.est[i,]),phis=matrix(mlphi.est[i,]), n.window=n.window,n.points=n.points,beta0s=beta1)$markedprocess
}

W <- owin(c(0, 1), c(0, 1))

samplesmlgcppp <- lapply(sampleMLGCP, function(z) ppp(x = z$x,
                                                      y = z$y,
                                                      marks = z$marks))

#########################################################################################
###Section 4 
####Fitting PPJSDM model to the simulated original MLGCP

library(ppjsdm)

nspecies <- 5

#Define initial parameters
short_range <- list(matrix(0.05, nrow=nspecies,ncol=nspecies), matrix(0.2, nrow=nspecies,ncol=nspecies))
model <- list("exponential", "exponential")

#Define matrices and arrays to store the parameter estimates
ppbeta0.est <- matrix(NA, nrow=nsim, ncol=nspecies)
ppbeta.est <- matrix(NA, nrow=nsim, ncol=nspecies)
ppbeta2.est <- matrix(NA, nrow=nsim, ncol=nspecies)
ppalpha.est <- array(NA, dim=c(nspecies, nspecies, nsim))
ppalpha2.est <- array(NA, dim=c(nspecies, nspecies, nsim))
ppaic.est <- matrix(NA,nrow=nsim, ncol=1)
ppbic.est <- matrix(NA,nrow=nsim, ncol=1)

#Combine the simulated covariate and background intensity to make two covariates in PPJSDM setting
cv <- list(as.im(cov,W=owin()),as.im(background,W=owin()))

#Initialize the parameters 
ppwindow=ppjsdm::Rectangle_window(c(0,1),c(0,1))
dummy_factor = 5 
min_dummy = 1000
saturation = 5
fitting_package = "glmnet"
dummy_distribution = "stratified"

for (i in 1:nsim)
{
  print(i)
  configuration <- ppp(x = X[[i]]$markedprocess$x,y = X[[i]]$markedprocess$y,marks = X[[i]]$markedprocess$marks)
  fit_ppjsdm <- ppjsdm::gibbsm(configuration,
                               covariates = cv,
                               window=ppwindow,
                               dummy_factor = dummy_factor, 
                               model=model,
                               min_dummy = min_dummy,
                               saturation = saturation, 
                               fitting_package = fitting_package,
                               dummy_distribution = dummy_distribution,
                               short_range = short_range,
                               debug = TRUE)
  
  ppbeta0.est[i, ] <- fit_ppjsdm$coefficients$beta0
  ppbeta.est[i, ] <- fit_ppjsdm$coefficients$beta[, 1]
  ppbeta2.est[i, ] <- fit_ppjsdm$coefficients$beta[, 2]
  ppalpha.est[, , i] <- fit_ppjsdm$coefficients$alpha[[1]]
  ppalpha2.est[, , i] <- fit_ppjsdm$coefficients$alpha[[2]]
  ppaic.est[i,] <- fit_ppjsdm$aic
  ppbic.est[i,] <- fit_ppjsdm$bic
  cat(i, "\t", ppbeta0.est[i, ], "\n", ppalpha.est[, ,i], "\n", ppaic.est[i,], "\n",ppbic.est[i,],"\n")
}

#Generate the average estimates
ppbeta0 <- apply(ppbeta0.est, 2, mean)

ppbeta1 <- apply(ppbeta.est, 2, mean)
ppbeta2 <- apply(ppbeta2.est, 2, mean)
#Correlation structure
ppalpha <- apply(ppalpha.est, c(1, 2), mean)
ppalpha2 <- apply(ppalpha2.est, c(1, 2), mean)

##########################################################################################
###Section 5
####Simulate samples from the fitted PPJSDM model
point_wiselg.s <- list()

steps <- 1e4
nthreads <- 4

for(i in 1:nsim){
  print(i)
  configuration <- ppp(x = X[[i]]$markedprocess$x,y = X[[i]]$markedprocess$y,marks = X[[i]]$markedprocess$marks)
  point_wiselg.s[[i]] <- ppjsdm::rgibbs( starting_configuration = as.Configuration(configuration),
                                         alpha = list(ppalpha.est[, , i], ppalpha2.est[, , i]),
                                         short_range = short_range, 
                                         beta0 = ppbeta0.est[i,], 
                                         beta = cbind(ppbeta.est[i, ], ppbeta2.est[i, ]),
                                         steps = steps,
                                         model=model,
                                         types = levels(configuration$marks),
                                         nsim=1, 
                                         saturation = saturation,
                                         nthreads = nthreads,
                                         covariates = cv,
                                         debug = TRUE)
}

#Check convergence
ppjsdm::trace_plot(point_wiselg.s[[1]], 1)
ppjsdm::trace_plot(point_wiselg.s[[1]], 2)

#Convert to spatstat format
W <- owin(c(0, 1), c(0, 1))

samplesppp <- lapply(point_wiselg.s, function(z) ppp(x = z$x,
                                                     y = z$y,
                                                     marks = z$types))
#########################################################################################
###Section 6
####Compute the total Kcross functions for the empirical(data) and fittted MLGCP \& PPJSDM models

##Compute the empirical K functions from data
Xprocess <- list()
for(i in 1:nsim){
  Xprocess[[i]] <- X[[1]]$markedprocess
}

W <- owin(c(0, 1), c(0, 1))

samples <- lapply(Xprocess, function(z) ppp(x = z$x,
                                            y = z$y,
                                            marks = z$marks,window=W))

r <- seq(0, 0.1, length=1e3)  

Xen <- vector(mode = "list", length = nspecies)
Ben <- vector(mode = "list", length = nspecies)
mlen <- vector(mode = "list", length = nspecies)
lower <- vector(mode = "list", length = nspecies)
upper <- vector(mode = "list", length = nspecies)
enlower <- vector(mode = "list", length = nspecies)
enupper <- vector(mode = "list", length = nspecies)
en <- vector(mode = "list", length = nspecies)

for(i in 1:nspecies) {
  Xen[[i]] <- vector(mode = "list", length = i)
  mlen[[i]] <- vector(mode = "list", length = i)
  en[[i]] <- vector(mode = "list", length = i)
  lower[[i]] <- vector(mode = "list", length = i)
  upper[[i]] <- vector(mode = "list", length = i)
  enlower[[i]] <- vector(mode = "list", length = i)
  enupper[[i]] <- vector(mode = "list", length = i)
  Ben[[i]] <- vector(mode = "list", length = i)
  
  for(j in 1:nspecies) {
    
    #Estimate baseline K functions
    temp <- pi*(r^2)
    Ben[[i]][[j]] <- temp
    
    nrank <- 10 # This is how many outliers are removed, the larger it is the tighter the band
    # Spatstat default is 1
    
    #Estimate empirical K functions
    temp <- lapply(samples, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    Xen[[i]][[j]] <- Reduce("+", temp) / length(temp)
    
    #Fitted MLGCP K functions
    temp <- lapply(samplesmlgcppp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    lower[[i]][[j]] <- sapply(seq_len(length(temp[[1]])), function(time) {
      vec <- sapply(temp, function(sample) {
        sample[time]
      })
      sort(vec)[1 + nrank]
    })
    
    upper[[i]][[j]] <- sapply(seq_len(length(temp[[1]])), function(time) {
      vec <- sapply(temp, function(sample) {
        sample[time]
      })
      sort(vec)[length(vec) - nrank]
    })
    
    mlen[[i]][[j]] <- sapply(seq_len(length(temp[[1]])), function(time) {
      vec <- sapply(temp, function(sample) {
        sample[time]
      })
      sort(vec)[length(vec) / 2]
    })
    
    #FittedPPJSDM  functions
    temp <- lapply(samplesppp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    #en[[i]][[j]] <- Reduce("+", temp) / length(temp)
    
    enlower[[i]][[j]] <- sapply(seq_len(length(temp[[1]])), function(time) {
      vec <- sapply(temp, function(sample) {
        sample[time]
      })
      sort(vec)[1 + nrank]
    })
    
    enupper[[i]][[j]] <- sapply(seq_len(length(temp[[1]])), function(time) {
      vec <- sapply(temp, function(sample) {
        sample[time]
      })
      sort(vec)[length(vec) - nrank]
    })
    
    en[[i]][[j]] <- sapply(seq_len(length(temp[[1]])), function(time) {
      vec <- sapply(temp, function(sample) {
        sample[time]
      })
      sort(vec)[length(vec) / 2]
    })
  }
}

i=2
j=5

#Plot the computed K functions

ggplot() +
  geom_ribbon(aes(x = r, y = en[[i]][[j]], ymin = enlower[[i]][[j]], ymax = enupper[[i]][[j]]), alpha = 0.2,fill='lightblue') +
  geom_ribbon(aes(x = r, y = mlen[[i]][[j]], ymin = lower[[i]][[j]], ymax = upper[[i]][[j]]), alpha = 0.2,fill='darkseagreen2') +
  geom_line(aes(x = r,y=mlen[[i]][[j]],colour = "MLGCP estimate" ),colour="chartreuse4") +
  geom_line(aes(x = r, y = en[[i]][[j]], colour = "PPJSDM estimate"),colour="deepskyblue3") +
  geom_line(aes(x = r,y= Xen[[i]][[j]],colour = "empirical K function" ),colour="purple") +
  geom_line(aes(x = r, y = pi * r^2, colour = "Baseline K function")) +
  theme_minimal(base_size = 15) +  #xlim(0,0.03) + ylim(0,0.003)+
  theme(legend.position = 'none')  +
  ylab(paste0("K", i, j)) +  xlab("r")  

#save(samplesppp,samples, samplesmlgcppp,background,point_wiselg.s, 
 #    ppbeta0.est, ppbeta.est, ppbeta2.est, ppalpha.est, ppalpha2.est,ppaic.est,
  #   ppbic.est, sampleMLGCP, mlxi.est,mlalpha.est, mlsig.est,mlphi.est, X,
   #  file = "5species-simulation-redo_option2.RData")
#########################################################################################
###Section 7
####MISE Computation
base.cross.K.est <- c(Ben[[1]][[2]],Ben[[1]][[3]],Ben[[1]][[4]],Ben[[1]][[5]],
                      Ben[[2]][[3]],Ben[[2]][[4]],Ben[[2]][[5]],Ben[[3]][[4]],Ben[[3]][[5]],Ben[[4]][[5]])

base.K.est <- c(Ben[[1]][[1]],Ben[[2]][[2]],Ben[[3]][[3]],Ben[[4]][[4]],Ben[[5]][[5]])

data.cross.K.est <- c(Xen[[1]][[2]],Xen[[1]][[3]],Xen[[1]][[4]],Xen[[1]][[5]],
                      Xen[[2]][[3]],Xen[[2]][[4]],Xen[[2]][[5]],Xen[[3]][[4]],Xen[[3]][[5]],Xen[[4]][[5]])
data.K.est <- c(Xen[[1]][[1]],Xen[[2]][[2]],Xen[[3]][[3]],Xen[[4]][[4]],Xen[[5]][[5]])

ppjsdm.cross.K.est <- c(en[[1]][[2]],en[[1]][[3]],en[[1]][[4]],en[[1]][[5]],
                        en[[2]][[3]],en[[2]][[4]],en[[2]][[5]],en[[3]][[4]],en[[3]][[5]],en[[4]][[5]])
ppjsdm.K.est <- c(en[[1]][[1]],en[[2]][[2]],en[[3]][[3]],en[[4]][[4]],en[[5]][[5]])

mlgcp.cross.K.est <- c(mlen[[1]][[2]],mlen[[1]][[3]],mlen[[1]][[4]],mlen[[1]][[5]],
                       mlen[[2]][[3]],mlen[[2]][[4]],mlen[[2]][[5]],mlen[[3]][[4]],mlen[[3]][[5]],mlen[[4]][[5]])
mlgcp.K.est <- c(mlen[[1]][[1]],mlen[[2]][[2]],mlen[[3]][[3]],mlen[[4]][[4]],mlen[[5]][[5]])

for(i in 1:nspecies) {
  for(j in 1:i) {
    if(i != j) {
      data.cross.K.est <- rbind(data.cross.K.est, Xen[[i]][[j]])
      base.cross.K.est <- rbind(base.cross.K.est, Ben[[i]][[j]])
      ppjsdm.cross.K.est <- rbind(ppjsdm.cross.K.est, en[[i]][[j]])
      mlgcp.cross.K.est <- rbind(mlgcp.cross.K.est, mlen[[i]][[j]])
    } else {
      data.K.est <- rbind(data.K.est, Xen[[i]][[i]])
      base.K.est <- rbind(base.K.est, Ben[[i]][[i]])
      ppjsdm.K.est <- rbind(ppjsdm.K.est, en[[i]][[i]])
      mlgcp.K.est <- rbind(mlgcp.K.est, mlen[[i]][[i]])
    }
  }
}
base_mise_btwn <- sapply(seq_len(nrow(base.cross.K.est)), function(i) mean((base.cross.K.est[i,]- data.cross.K.est[i,] )^2))
base_mise_within <- sapply(seq_len(nrow(base.K.est)), function(i) mean((base.K.est[i,]- data.K.est[i,] )^2))

mean(base_mise_btwn)
mean(base_mise_within)
mean(c(base_mise_btwn, base_mise_within))

ppjsdm_mise_btwn <- sapply(seq_len(nrow(data.cross.K.est)), function(i) mean((ppjsdm.cross.K.est[i,]- data.cross.K.est[i,] )^2))
ppjsdm_mise_within <- sapply(seq_len(nrow(ppjsdm.K.est)), function(i) mean((ppjsdm.K.est[i,]- data.K.est[i,] )^2))

mean(ppjsdm_mise_btwn)
mean(ppjsdm_mise_within)
mean(c(ppjsdm_mise_btwn, ppjsdm_mise_within))

mlgcp_mise_btwn <- sapply(seq_len(nrow(data.cross.K.est)), function(i) mean((mlgcp.cross.K.est[i,]- data.cross.K.est[i,] )^2))
mlgcp_mise_within <- sapply(seq_len(nrow(mlgcp.K.est)), function(i) mean((mlgcp.K.est[i,]- data.K.est[i,] )^2))

mean(mlgcp_mise_btwn)
mean(mlgcp_mise_within)
mean(c(mlgcp_mise_btwn,mlgcp_mise_within))

