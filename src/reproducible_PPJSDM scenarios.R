remove(list = ls())
library(ppjsdm)
library(spatstat)
library(RandomFields)
RFoptions(install="no")
library(Multilogreg)
library(ggplot2)

#reproducibility
seed <- 1
set.seed(seed)

#Define the parameters of the PPJSDM model interested to simulate
#interested window
window <- Rectangle_window(c(0, 1), c(0, 1))
#How many simulations are required
nreplications <- 5 # Change to 100 for full version
#Number of species
ntypes <- 2
#number of steps in the MCMC
steps <- 1e4 # Change to 1e5 for full version

#intercept
beta0 <- c(5, 4.1)
#Define the short range interaction coefficients
alpha <- cbind(c(-0.4, -0.1), c(-0.1, 0.3))

#Define a single covariate
covariates <- list(temperature = function(x, y) x - 0.5)
#Coefficent of the covariate
beta <- cbind(c(-1, 1))

#poteintial model
model <- "square_bump"
saturation <- 2
nthreads <- 4

short_range <- matrix(0.05, ntypes, ntypes)

#Simulate from PPJSDM model
samples <- ppjsdm::rgibbs(window = window,
                          beta0 = beta0,
                          alpha = alpha,
                          short_range = short_range,
                          nsim = nreplications,
                          model = model,
                          saturation = saturation,
                          covariates = covariates,
                          beta = beta,
                          drop = FALSE,
                          nthreads = nthreads,
                          steps = steps)

#Check the number of points in each species
samples[[1]][1]
samples[[1]][2]

# Check convergence of number of points for both species
ppjsdm::trace_plot(samples[[1]], 1)
ppjsdm::trace_plot(samples[[1]], 2)

#Fitting MLGCP model to the simulated PPJSDM scenarios

#Define MLGCP parameters
xib = c(0.03,0.01)
sib = c(0.8,0.1)
phb = c(0.03,0.01)

#Number of latent fields
latent = 1
Rmax   = 1
rl     = 1e4
r      = seq(0,Rmax,length=rl)

par.est <- list()
Betahat <- list()

# Size of the observation window
n.x <- n.y <- 1 #Assuming a square window that starts at (0,0)
xx=seq(0,n.x,length=128)
yy=seq(0,n.x,length=128)

#Conver the covariate into an im object
C <- as.im(covariates$temperature, W = owin())
cov <- outer(xx, yy, FUN = function(x, y) covariates$temperature(x, y))

##############################################################################
# Parameter estimation
for ( i in 1:nreplications){
  print(i)
  FC <- as.function(as.im(t(cov),W=owin()))
  cv <- as.matrix(cbind(1,FC(samples[[i]]$x,samples[[i]]$y)),ncol=2)
  configuration <- ppp(x = samples[[i]]$x,y = samples[[i]]$y,marks = samples[[i]]$types)
  
  #Choose initial values for Beta coefficients
  #Use the PPJSDM model to derive accurate estimates of intercept and coefficient of the covariate
  fit <- gibbsm(samples[[i]], short_range = matrix(0, 2, 2),  min_dummy = 1e3, dummy_distribution = "stratified")
  
  beta1<- c(0,0); beta1 <- as.matrix(beta1)
  beta0 <- coef(fit)$beta0
  Beta <- cbind(beta0,beta1)
  
  Betahat[[i]] <- FirstOrderCCL(X=configuration,Beta0=Beta,covariate = C)
  
  par.est[[i]] <- PenalizedSecondOrderCCL(X=configuration,covariate = cv,Beta = Betahat[[i]]$betahat,Rmax = Rmax,
                                          xibound=xib,sigmabound = sib,phibound=phb, lat = latent,tol=1e-5,lamb=0)
}

#matrices and arrays to store Parameter estimates
mlxi.est <- matrix(NA, nrow=nreplications, ncol=latent)
mlalpha.est <- array(NA, dim=c(ntypes, latent, nreplications))
mlsig.est <- matrix(NA, nrow=nreplications, ncol=ntypes)
mlphi.est <- matrix(NA, nrow=nreplications, ncol=ntypes)

for(i in 1:nreplications){
  mlxi.est[i,]  <- round(par.est[[i]][[1]]$xi,3)
  mlalpha.est[, , i] <- round(par.est[[i]][[1]]$alpha[,,1],3)
  mlsig.est[i,] <- round(par.est[[i]][[1]]$sigma,3)
  mlphi.est[i,] <- round(par.est[[i]][[1]]$phi,3)
}

#Average MLGCP parameter estimates
apply(mlxi.est, 2, mean)
apply(mlsig.est, 2, mean)
apply(mlphi.est, 2, mean)
apply(mlalpha.est, c(1, 2), mean)

# Size of the observation window (assume a s  quare window that starts at (0,0))
n.window <- n.x

type <- c('type1','type2')

#Compute accurate n.points for MLGCP simulations form the fitted model
points.est <-matrix(NA, nrow=nreplications,ncol=ntypes)

for(i in 1:length(samples)){
  configuration <- ppp(x = samples[[i]]$x,y = samples[[i]]$y,marks = samples[[i]]$types)
  species1 <-subset(configuration, marks == type[1])
  species2 <-subset(configuration, marks == type[2])
  
  points.est[i,] <- c(species1$n,species2$n)
  
}
n.points <- apply(points.est,2,mean)

#Simulate from the fitted MLGCP model
fit <- list()
intensity <- list()
density <- list()
sampleMLGCP <- list()
for (j in 1:nreplications) {
  print(j)
  basecov <- 0
  configuration <- ppp(x = samples[[j]]$x,y = samples[[j]]$y,marks = samples[[j]]$types)
  
  for ( i in 1:ntypes) {
    #Estimation of basecov/rho_0/background intensity as explained in Hessellund et. al (2020) supplementary doc
    spec <- subset(configuration, marks == type[i])
    fit[[i]] <- ppm(spec ~ 1+C)
    
    #estimate the bandwidth
    bw =  bw.CvL(spec)
    #Intensity of each of the species.
    intensity[[i]] <- predict(fit[[i]], type = 'intensity', locations = spec ) 
    density[[i]] <- density(spec,weights = (1/intensity[[i]]),sigma=bw)
    basecov <- basecov + density[[i]]
  }
  background <- basecov/ntypes
  
  #Simulate from the fitted model
  sampleMLGCP[[j]] <- sim_lgcp_multi(basecov= as.matrix(log(background)),
                                     covariate=cov ,betas=matrix(c(1,1)),alphas=matrix(mlalpha.est[,,1],ncol=1),
                                     xi=mlxi.est[j,],
                                     sigma=matrix(mlsig.est[j,]),phis=matrix(mlphi.est[j,]), n.window=n.window,
                                     n.points=n.points,beta0s=NULL)$markedprocess
}

#Convert the fitted MLGCP samples and original PPJSDM samples in to ppp
W <- owin(c(0, 1), c(0, 1))

samplesmlgcppp <- lapply(sampleMLGCP, function(z) ppp(x = z$x,
                                                      y = z$y,
                                                      marks = z$marks))
samplesdatappp <- lapply(samples, function(cc) ppp(x = cc$x,
                                                   y = cc$y,
                                                   marks = cc$types))

for(i in 1:nreplications){
  levels(samplesdatappp[[i]]$marks) <- c("1", "2")
  
}

r <- seq(0, 0.1, length = 1e3)  

mlen <- vector(mode = "list", length = ntypes)
en <- vector(mode = "list", length = ntypes)

for(i in 1:ntypes) {
  print(i)
  mlen[[i]] <- vector(mode = "list", length = i)
  en[[i]] <- vector(mode = "list", length = i)
  for(j in 1:ntypes) {
    print(j)
    #Fitted MLGCP K functions
    temp <- lapply(samplesmlgcppp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    mlen[[i]][[j]] <- Reduce("+", temp) / length(temp)
    
    #Empirical PPJSDM  functions
    temp <- lapply(samplesdatappp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    en[[j]][[i]] <- Reduce("+", temp) / length(temp)
  }
  
}

ggplot() + 
  geom_line(aes(x = r, y = en[[2]][[2]], colour = "PPJSDM estimate")) +
  geom_line(aes(x = r,y=mlen[[2]][[2]],colour = "MLGCP estimate" )) +
  geom_line(aes(x = r, y = pi * r^2, colour = "Baseline")) +
  theme_minimal(base_size = 15) +  #xlim(0,0.03) + ylim(0,0.003)+
  theme(legend.title = element_blank()) +
  ylab(paste0("K", i, j)) +  xlab("r")  
