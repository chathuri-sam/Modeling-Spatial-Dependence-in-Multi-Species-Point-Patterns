remove(list = ls())
library(ppjsdm)
library(spatstat)
library(RandomFields)
RFoptions(install="no")
library(Multilogreg)
library(ggplot2)

# Reproducibility
set.seed(1)

#define the size of the observation window, for x axis (xx) and y axis (yy)
xx=seq(0,1,length=100)
yy=seq(0,1,length=100)

# Simulate a single covariate (need to be in the matrix form)
cov <- as.matrix(RFsimulate(RMexp(var=1,scale=0.05), x=xx, y=yy, grid=TRUE))

# Simulate the background intensity as given in the paper
gamma <- 0.8
background <- as.matrix(RFsimulate(RMgauss(var=1,scale=0.2), x=xx, y=yy, grid=TRUE))*gamma-gamma^2/2

#Set up parameters beta0 and beta1 for intercept and covariate coefficient
beta0 <- rep(0, 2)
beta1 <- c(0.5,-0.02)
beta1 <- as.matrix(beta1)

#Set up the Parameters in the MLGCP
alpha <- matrix(c(-0.7 ,0.7,-0.5,0.5), 2, 2)
xi <- c(0.5)
sigma <- matrix(c(0.01,0.02))
phi <- matrix(c(0.01,0.01))

#Define the number of simulations needed
nsim = 10 # Change to 100 for full version

#Define the upper limit of the window (assuming a sqaure window of interest)
n.x <- 1 
#equate the n.x to the limits of the interested square window
n.window <- n.x

#define the expected number of points for each species required
n.points <- c(100,100)

# Simulation of a multivariate LGCP (using Mulilogreg package)
X <- list()
for(i in 1:nsim){
  print(i)
  X[[i]] <- sim_lgcp_multi(basecov=background,covariate=cov,betas=beta1,alphas=alpha,xi=xi,
                           sigma=sigma,phis=phi, n.window=n.window,n.points=n.points,beta0s=beta0)$markedprocess
  
}

#Define the matrices and arrays for the parameters for fitting model
#with two different short range interaction distances
nspecies <- 2
ppbeta0.est <- matrix(NA, nrow=nsim, ncol=nspecies)
ppbeta.est <- matrix(NA, nrow=nsim, ncol=nspecies)
ppbeta2.est <- matrix(NA, nrow=nsim, ncol=nspecies)
ppalpha.est <- array(NA, dim=c(nspecies, nspecies, nsim))
ppalpha2.est <- array(NA, dim=c(nspecies, nspecies, nsim))
ppaic.est <- matrix(NA,nrow=nsim, ncol=1)
ppbic.est <- matrix(NA,nrow=nsim, ncol=1)

#Define the potential short range intereaction distances
short_range <- list(matrix(0.05, nrow=2,ncol=2), matrix(0.2, nrow=2,ncol=2))
#Define the potenital model options for the two short range distances
model <- list("exponential", "exponential")

#combine the simulated covariate and background intensity to be used as 2 covariates in PPJSDM model fit
cv <-  list(as.im(cov, W=owin()), as.im(background, W=owin()) )

#Define the potential parameters for PPJSDM model fit
ppwindow=ppjsdm::Rectangle_window(c(0,1),c(0,1))
dummy_factor = 5
min_dummy = 4000
model=model
nthreads = 4
saturation = 4
fitting_package = "glmnet"
dummy_distribution = "stratified"

#Fitting of PPJSDM model
for (i in 1:nsim)
{
  fit_ppjsdm<- ppjsdm::gibbsm(X[[i]],
                              window=ppwindow,
                              dummy_factor = dummy_factor,
                              min_dummy = min_dummy,
                              model=model,
                              covariates = cv ,
                              nthreads = nthreads,
                              saturation = saturation,
                              fitting_package = fitting_package,
                              dummy_distribution = dummy_distribution,
                              short_range = short_range)
  ppbeta0.est[i, ] <- fit_ppjsdm$coefficients$beta0
  ppbeta.est[i, ] <- fit_ppjsdm$coefficients$beta[, 1]
  ppbeta2.est[i, ] <- fit_ppjsdm$coefficients$beta[, 2]
  ppalpha.est[, , i] <- fit_ppjsdm$coefficients$alpha[[1]]
  ppalpha2.est[, , i] <- fit_ppjsdm$coefficients$alpha[[2]]
  ppaic.est[i,] <- fit_ppjsdm$aic
  ppbic.est[i,] <- fit_ppjsdm$bic
  cat(i, "\n")
}

#Define the number of steps required to be used in the MCMC
steps <- 1e4 # Change to 1e5 for full version
lg.s <- list()

#Simulate samples from the fitted PPJSDM model
for (i in 1:nsim) {
  print(i)
  lg.s[[i]] <- ppjsdm::rgibbs( alpha = list(ppalpha.est[, , i], ppalpha2.est[, , i]),
                               short_range = short_range,
                               beta0 = ppbeta0.est[i, ],
                               beta = cbind(ppbeta.est[i, ], ppbeta2.est[i, ]),
                               steps = steps,
                               types = levels(X[[1]]$marks),
                               covariate = cv ,
                               model = model,
                               saturation = saturation,
                               nthreads = nthreads,
                               starting_configuration = as.Configuration(X[[i]]),
                               debug = TRUE)
}


# Check convergence of number of points for both species
ppjsdm::trace_plot(lg.s[[1]], 1)
ppjsdm::trace_plot(lg.s[[1]], 2)

#The number of points for each species in the samples
points_1 <- sapply(lg.s, function(ss) length(ss[1]$x))
points_2 <- sapply(lg.s, function(ss) length(ss[2]$x))

#Convert the samples into point patterns
samplesppp <- lapply(lg.s , function(z) ppp(x = z$x,
                                            y = z$y,
                                            marks = z$types))[points_1 > 0 & points_2 > 0]

#Convert the original MLGCP samples into point patterns as well
samplesdatappp <- lapply(X , function(z) ppp(x = z$x,
                                             y = z$y,
                                             marks = z$marks))

#Define the interested range of r for the K functions
r <- seq(0,0.1, length=1e3)

mlen <- vector(mode = "list", length = nspecies)
en <- vector(mode = "list", length = nspecies)

for(i in 1:nspecies) {
  mlen[[i]] <- vector(mode = "list", length = i)
  en[[i]] <- vector(mode = "list", length = i)
  for(j in 1:nspecies) {
    #empirical MLGCP K functions
    temp <- lapply(samplesdatappp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    mlen[[i]][[j]] <- Reduce("+", temp) / length(temp)
    
    #Fitted PPJSDM  functions
    temp <- lapply(samplesppp, function(cc) {   
      z <- Kcross(cc, i = i, j = j, correction = "best", r = r) # Choose best border correction, other options irrelevant
      z[[which(!(names(z) %in% c("r", "theo")))]] # This selects the K function estimate corresponding to that border correction
    }) 
    
    en[[i]][[j]] <- Reduce("+", temp) / length(temp)
  }
}

#Plot the computed K functions
ggplot() + 
  geom_line(aes(x = r, y = en[[2]][[1]], colour = "PPJSDM estimate")) +
  geom_line(aes(x = r,y=mlen[[2]][[1]],colour = "Empirical MLGCP" )) +
  geom_line(aes(x = r, y = pi * r^2, colour = "Baseline")) +
  theme_minimal(base_size = 15) +  #xlim(0,0.03) + ylim(0,0.003)+
  theme(legend.title = element_blank()) +
  ylab(paste0("K", i, j)) +  xlab("r")  
