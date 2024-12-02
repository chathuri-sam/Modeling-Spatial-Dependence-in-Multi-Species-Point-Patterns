library(ecespa)
library(ggplot2)
library(Multilogreg)
library(ppjsdm)
library(RandomFields)
library(spatstat)

RFoptions(install="no")

set.seed(1)

# Load fixed version of simulation code
source("Sim_MLGCP_adjusted.R")

###MLGCP model fit for the Swamp follows the method used in Hasselhund et al. (2022)
###The following code is adapted by Hasselhund et al. (2022)
### First order analysis ###
nspecies <- length(unique(swamp$sp))
swamp$type <- as.numeric(swamp$sp)
swamp$type <- as.factor(swamp$type)

pointprocess <- ppp(x = swamp$y,y = swamp$x,marks = swamp$type,window = owin(xrange = c(0,200),yrange = c(0,50)) )

window <- ppjsdm::Rectangle_window(c(0, 200), c(0, 50))

xx <- seq(0,200, length=128)
yy <- seq(0,50,length=128)

#Generate artificial covariate
covariates <- list(horizontal = function(x, y) x / 100)

FC <- covariates$horizontal(pointprocess$x,pointprocess$y)
cv <- as.matrix(cbind(1,FC),ncol=2)

ncovs <- dim(cv)[2]

# Initial parameters
Beta0=matrix(runif(nspecies*ncovs,-0.5,0.5),nrow =nspecies)

# Choosing the last type of species as the baseline
Beta0 <- as.matrix(t(scale(Beta0,center=Beta0[nspecies,],scale=FALSE)))

# Parameter estimation

C <- as.im(covariates$horizontal, W = window)

Betahat=FirstOrderCCL(X=pointprocess,Beta0 = Beta0,covariate = C)

# Parameter estimates
colnames(Betahat$betahat)=c("FX", "NS", "NX","TD", "OT")
rownames(Betahat$betahat)=c("Intercept", "Horizontal")

Betahat$betahat

### Second order analysis###

# Globals
len.x= window$x_range[2]- window$x_range[1]
len.y= window$y_range[2]-window$y_range[1]
lw=(len.x+len.y)/2
latent = 2
Rmax   = 200
rl     = 100
r      = seq(0,Rmax,length=rl)

# Parameter spaces (upper,lower)
xib = c(lw,lw/500)
sib = c(5,0.01)
phb = c(lw/10,lw/1000)

# Parameter estimation

par.est <-  PenalizedSecondOrderCCL(X=pointprocess ,covariate = cv,Beta = Betahat$betahat, 
                                    Rmax = Rmax,xibound=xib,sigmabound = sib,phibound=phb, lat = latent,lamb = 2.5,tol=1e-5)

# Parameter estimates
mlalpha.est = matrix(par.est[[1]]$alpha[,,1],ncol=latent)
mlxi.est  = par.est[[1]]$xi
mlsig.est = par.est[[1]]$sigma
mlphi.est = par.est[[1]]$phi

#Compute accurate n.points
species1 <-subset(pointprocess, marks == "1")
species2 <-subset(pointprocess, marks == "2")
species3 <-subset(pointprocess, marks == "3")
species4 <-subset(pointprocess, marks == "4")
species5 <-subset(pointprocess, marks == "5")

n.points <- c(species1$n,species2$n,species3$n,species4$n,species5$n)

cov<-outer(xx, yy, FUN = function(x, y) covariates$horizontal(x, y))

n.window <- owin(xrange = c(0,200),yrange = c(0,50))

#Estimate the background intensity
fit <- list() 
intensity <- list() 
density <- list()
basecov <- 0
for ( i in 1:nspecies) {    
  spec <- subset(pointprocess, marks ==i)   
  fit[[i]] <- ppm(spec ~ 1 + C)
  
  bw =  bw.CvL(spec)   
  #Intensity of each of the species.   
  intensity[[i]] <- predict(fit[[i]], type = 'intensity', locations = spec )   
  density[[i]] <- density(spec,weights = (1/intensity[[i]]) ,sigma=bw)   
  basecov <- basecov + density[[i]]    
} 

#Take the average of all the species intensities for the estimate of basecov background <- basecov/6
background <- basecov/6

plot(background,main="estimated background intensity")

#Derive samples from fitted model
nsim <- 100
X <- list()
for (i in 1:nsim){
  print(i)
  X[[i]] <- sim_lgcp_multi_fixed(basecov=log(background),covariate=cov,
                                 betas= as.matrix(Betahat$betahat[-1,]),
                                 alphas=matrix(mlalpha.est,ncol=latent),xis=mlxi.est,
                                 sigmas=mlsig.est,phis=mlphi.est, n.window=n.window,n.points=n.points,beta0s=NULL)$markedprocess
  
}


samplesmlgcppp <- lapply(X, function(cc) ppp(x = cc$x,
                                             y = cc$y,
                                             marks = cc$marks,window = owin(xrange = c(0,200),yrange = c(0,50)) ))


######################################################################
#Fit for the Swamp data using Saturated Pairwise Interaction Gibbs Point Process (Flint et all (2020))
library(ppjsdm)

#Define the parameters 
short_range <- matrix(5,5,5)
model <- "square_exponential"
ndummy <- 5e3
nthreads <- 4
saturation <- 2

nspecies <- length(unique(swamp$sp))
swamp$type <- as.numeric(swamp$sp)
swamp$type <- as.factor(swamp$type)

configuration <- ppjsdm::Configuration(swamp$y, swamp$x, swamp$type)

back <- as.im(background, W = owin(xrange = c(0,200),yrange = c(0,50)))
C <- as.im(covariates$horizontal, W = owin(xrange = c(0,200),yrange = c(0,50)))

#Fit the PPJSDM model
fit.pp <- ppjsdm::gibbsm(configuration, 
                         window = window, 
                         model = model,
                         covariates = list(C,back),
                         short_range = short_range,
                         fitting_package = "glm",
                         dummy_distribution = "stratified",
                         min_dummy = ndummy,
                         saturation = saturation,
                         nthreads = nthreads,
                         dummy_factor = 1)
plot(fit.pp)
fit.pp

#Compute summary statistics
summ <- summary(fit.pp) 
summ
box_plot(fit.pp, summ = summ,which="within",full_names = c(`1`="FX", `2`="NS", `3`="NX",`4`="OT", `5`="TD"))
plot(potentials(fit.pp, 1, 1))

#Simulate samples from the fitted model
nsim <- 100

samples <- ppjsdm::rgibbs(fit.pp,
                          steps = 1e4, # Number of steps in the Metropolis-Hastings algorithm
                          nsim = nsim,
                          nthreads = 4,
                          debug=TRUE) 


samplesppp <- lapply(samples, function(cc) ppp(x = cc$x,
                                               y = cc$y,
                                               marks = as.factor(as.numeric(cc$types)),window = owin(xrange = c(0,200),yrange = c(0,50)) ))


#Conditional intensity predictions
type <- 5
plot_papangelou(fit.pp, type = type, show = type, use_log = TRUE, drop_type_from_configuration = TRUE,
                window = owin(xrange = c(0,200),yrange = c(0,50)),
                type_description = "Tree Species",
                legend_title = "Log-Cond.int.")

#Compute the AUC
aucs <- sapply(1:5, function(i) {
  conditional_intensity <- plot_papangelou(fit.pp, type = i, show = i, return_papangelou = TRUE,drop_type_from_configuration = TRUE)
  conditional_intensity$v[is.na(conditional_intensity$v)] <- 0.
  conf <- ppp(x = configuration$x,
              y = configuration$y,
              marks = configuration$types,
              window = owin(xrange = c(0,200),yrange = c(0,50)))
  z <- subset(conf, conf$marks == i)
  z <- subset(z, inside.owin(z$x, z$y, conditional_intensity))
  auc(X = z, covariate = as.function(conditional_intensity))
})

aucs

########################################################################
########################################################################
#Compute K function

Xen <- vector(mode = "list", length = nspecies)
Ben <- vector(mode = "list", length = nspecies)
mlen <- vector(mode = "list", length = nspecies)
lower <- vector(mode = "list", length = nspecies)
upper <- vector(mode = "list", length = nspecies)
enlower <- vector(mode = "list", length = nspecies)
enupper <- vector(mode = "list", length = nspecies)
en <- vector(mode = "list", length = nspecies)

i=j=1
temp <- Kcross(pointprocess,i=i,j=j,correction = "best")
r <- temp$r

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
    
    #Estimate empirical K functions
    temp <- Kcross(pointprocess,i=i,j=j,correction = "best")
    Xen[[i]][[j]] <- temp$iso
    #lo <- loess(temp$iso ~ r, span = 0.05)
    #Xen[[i]][[j]] <- lo
    
    #Fitted MLGCP K functions
    temp <- lapply(samplesmlgcppp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    nrank <- 10 # This is how many outliers are removed, the larger it is the tighter the band
    # Spatstat default is 1
    
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

i=1
j=1

#Plot the computed K functions

ggplot() +
  geom_ribbon(aes(x = r, y = en[[i]][[j]], ymin = enlower[[i]][[j]], ymax = enupper[[i]][[j]]), alpha = 0.2,fill='lightblue') +
  geom_ribbon(aes(x = r, y = mlen[[i]][[j]], ymin = lower[[i]][[j]], ymax = upper[[i]][[j]]), alpha = 0.2,fill='darkseagreen2') +
  geom_line(aes(x = r,y=mlen[[i]][[j]],colour = "MLGCP estimate" ),colour="chartreuse4") +
  geom_line(aes(x = r, y = en[[i]][[j]], colour = "PPJSDM estimate"),colour="deepskyblue3") +
  geom_line(aes(x = r,y= Xen[[i]][[j]],colour = "empirical K function" ),colour="#e66101") +
  geom_line(aes(x = r, y = pi * r^2, colour = "Baseline K function"),colour = "red") +
  theme_minimal(base_size = 15) +  #xlim(0,0.03) + ylim(0,0.003)+
  theme(legend.position = 'none')  +
  ylab(paste0("K", "FX", "FX")) +  xlab("r")  

#Plot papangelou conditional intensity
type <- 1

plot_papangelou(fit, type = type, show = type, use_log = TRUE, drop_type_from_configuration = TRUE,
                window = owin(xrange = c(0,200),yrange = c(0,50)),
                type_description = "Tree Species",
                legend_title = "Log-Cond.int.",
                base_size = 17,
                mark_range = 5,
                colours = "black")

#Compute the AUC
aucs <- sapply(1:5, function(i) {
  conditional_intensity <- plot_papangelou(fit, type = i, show = i, return_papangelou = TRUE,drop_type_from_configuration = TRUE)
  conditional_intensity$v[is.na(conditional_intensity$v)] <- 0.
  conf <- ppp(x = configuration$x,
              y = configuration$y,
              marks = configuration$types,
              window = owin(xrange = c(0,200),yrange = c(0,50)))
  z <- subset(conf, conf$marks == i)
  z <- subset(z, inside.owin(z$x, z$y, conditional_intensity))
  auc(X = z, covariate = as.function(conditional_intensity))
})


