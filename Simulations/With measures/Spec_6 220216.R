########### MIXTURE SIMULATIONS EXPERIMENT WITH MEASURES ###########################
#rm(list = ls())

# Here I simulate latent factors and measures before testing how well the proposed method works
### Steps:
#1. Simulate latent factors from RHS of CES function
#2. Invent some 'true' parameters and generate the LHS of CES function
#3. Generate measures using some 'true' loadings
#4. Use the EM algorithm to estimate the joint distribution of measures
#5. Use minimum distance to relate the estimated joint distribution of measures to and 
# estimated joint distribution of latent factors
#6. Draw from the estimated joint distribution of latent factors
#7. Estimate a CES function over these draws
#8. Organise output

# Once this script is run, we have both a set of CES estimates, and the EM output for each run

############ load packages

library(ks)
library("micEconCES")
library("corpcor")
library("mixtools")
#library('mclust')
library("fMultivar")
library(xlsx)
library(nortest)
library("Matrix")
library("LICORS")

#set seed
set.seed(8192)

### Load useful functions

source("Functions/mindistance.R") #minimum distance function

### Set up to loop over different values

# Set list of rho parameters to test over

rhos <- c(-1,-0.5,0,0.5,1)

# Set the number of runs desired for each rho parameter

#numsims <- 1

# Create empty list for storing EM step estimates

EMstep <- list()

# create empty vectors / matrices for storing estimates

rhoestimates <- matrix(nrow = length(rhos), ncol = numsims)
alphaestimates <-
  matrix(nrow = length(rhos), ncol = numsims)
betaestimates <- matrix(nrow = length(rhos), ncol = numsims)

# create empty vectors / matrices for storing estimates taken BEFORE simulation stage,
# ie just to test how well the non-linear least sqaures aspect does at finding parameters

rhoestimatespre <- matrix(nrow = length(rhos), ncol = numsims)
alphaestimatespre <-
  matrix(nrow = length(rhos), ncol = numsims)
betaestimatespre <- matrix(nrow = length(rhos), ncol = numsims)

# Start the clock!
ptm <- proc.time()

iteration_num <- 0 

for (i in 1:length(rhos))
{
  
  for (j in 1:numsims)
  {
    
iteration_num <- iteration_num +1
print("iteration number")
print(iteration_num)

########### 1. Simulate RHS latent factors
# Here I present several different options for the RHS latent factors in our CES function

### Joint logmixnormal Xs

t <- 0.5
mus <- rbind(c(0,1), c(0,2)) # first two entries are first mixture
Sigmas <- rbind(cbind(c(0.1,0.06), c(0.06, 0.08)), cbind(c(0.1,0.05), c(0.05, 0.1)))
props <- c(t, 1-t)

logx <- rmvnorm.mixt(1000, mus=mus, Sigmas=Sigmas, props=props)

########### Simulate error

### no error

### normal
u <- rnorm(1000, mean = 0, sd = 0.02) #for use with lognormal x's

### Cauchy (fat tails)

#u <- rcauchy(1000, location = 0, scale = 0.01)


### Chisquare (skew)

#u <- rchisq(1000,0.1)
#u <- 0.1*(u - mean(u))/var(u)

## mixnormal

#umus <- c(-0.2,0.2)
#usigs <- rbind(0.05,0.05)
#u <- rnorm.mixt(1000, mus = umus, sigmas= usigs, props=props)
#u <- u - mean(u)

#U <- exp(u)
#plot3d(logx[,1],logx[,2],u)



############## 2. Generate LHS variable
    
x <- exp(logx) #take exponential to plug into true production function
    
rho <- rhos[i]
#rho <- 1
alpha <- 0.5
beta <- 0.5
    
if (rho == 0) {
   y <- (exp(alpha * log(x[,1]) + beta * log(x[,2]) + (u)))
  } else {
    y <- ((alpha * (x[,1] ^ rho) + beta * (x[,2] ^ rho)) ^ (1 / rho)) * exp(u)
  }
    
#plot3d(x[,1],x[,2],y)
#plot3d(logx[,1],logx[,2],y)
    
### Estimate CES function on 'true' simulated data
# This should be able to recover estimates of parameters, provided the NLS is working
    
simdata.1 <- cbind(y,x)
simdata1 <- data.frame(x1 = y, x2 = x[,1], x3 = x[,2])
#plot(simdata1)
if (rho>=0)
  {
  cesNlstrue <-nlsLM(
        log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
        simdata1, control = nls.lm.control(maxiter = 1000), start = c(
        rho = 0.5, alpha = 0.7, beta = 0.3
          )
      )
  } else if (rho<0 & rho>=-1) {
    cesNlstrue <-nlsLM(
        log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
        simdata1, control = nls.lm.control(maxiter = 1000), start = c(
        rho = -0.3, alpha = 0.7, beta = 0.3
          )
      )
  } else {
      cesNlstrue <-nlsLM(
        log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
        simdata1, control = nls.lm.control(maxiter = 1000), start = c(
        rho = -1.5, alpha = 0.7, beta = 0.3
          )
      )
  }
    print(cesNlstrue)
    
estimpre <- summary(cesNlstrue)
    
rhoestimatespre[i,j] <- estimpre[[10]][1]
    
alphaestimatespre[i,j] <- estimpre[[10]][2]
    
betaestimatespre[i,j] <- estimpre[[10]][3]
    
# need to turn simulated data back into log form
    
simdata_1 <- log(simdata.1)
    
############## 3.Generate observable measures
# Now I generate the observed measures according to some true relationship to latent factors
    
### Factor loadings 
# I have 3 dedicated measurements for each latent factor, one with high variance, one with medium variance, and one with low variance
    
L <- rbind(cbind(1,0,0),cbind(1.8,0,0),cbind(0.4,0,0),
          cbind(0,1,0),cbind(0,0.3,0),cbind(0,3.1,0),
          cbind(0,0,1),cbind(0,0,1.5),cbind(0,0,0.8))
    
### Errors
# I assume nine types of (independent) errors with low, medium and high variances
    
epsilon1 <- rnorm(1000, mean = 0, sd = 0.0005)
epsilon2 <- rnorm(1000, mean = 0, sd = 0.001)
epsilon3 <- rnorm(1000, mean = 0, sd = 0.0015)
epsilon4 <- rnorm(1000, mean = 0, sd = 0.0005)
epsilon5 <- rnorm(1000, mean = 0, sd = 0.001)
epsilon6 <- rnorm(1000, mean = 0, sd = 0.0015)
epsilon7 <- rnorm(1000, mean = 0, sd = 0.0005)
epsilon8 <- rnorm(1000, mean = 0, sd = 0.001)
epsilon9 <- rnorm(1000, mean = 0, sd = 0.0015)
    
epsilon <- rbind(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6, epsilon7, epsilon8, epsilon9)
    
### Measures
#Combing factor loadings, measurement errors and latent factors to generate dataset of observable measures
# No intercept

#simdata1 <- data.frame(x1 = y, x2 = x[,1], x3 = x[,2])
    
measures_1 <- L%*%t(simdata_1)
for (p in 1:dim(measures_1)[[1]])
  {
  measures_1[p,] <- measures_1[p,] + epsilon[p,]
  }
    
measures_1 <- t(measures_1)
measures <- measures_1 
      
#demean
#for (i in 1:dim(measures_1)[[2]])
#{
#  measures[,i] <- measures_1[,i] - mean(measures_1[,i])
#}

measures <-as.data.frame(measures) 
      
############### 4.Estimate joint distribution of measures
# Here I use the EM algorithm to estimate the joint distribution of Measures   
 
#plot(simdata.2)
    
### Mixture distribution
# Set the starting parameters for the EM algorithm.

lambdastart <- c(0.5,0.5)

sigstart <- list(cov(measures), cov(measures))

kmeansplusplus <- kmeanspp(measures)

mustart <- list()
mustart[[1]] <- kmeansplusplus[[2]][1,1:dim(measures)[[2]]]
mustart[[2]] <- kmeansplusplus[[2]][2,1:dim(measures)[[2]]]

# Initial values for M step: this is a list with (sigma, mean, prop)  
    
#kmeans<-simdata2[complete.cases(simdata2),]
#out<-kmeans(kmeans,2) 
#mean0 <- out[[2]]
#cov.start <- cov(simdata2, use ="complete.obs")   
#sigma0 <- list(cov.start, cov.start)
#prop0 <- c(.5,.5)
#prop0 <- c(mean(out[[1]]-1),1-(mean(out[[1]]-1)))
#mstep.start <- list(mean0, sigma0, prop0)
    
# Number of mixtures              
nM <- 2
    
# Number of factors 
nF <-3
    
# Convergence criterion            
#conv <- 1e-08
    
# Number of measurements
nMeas <- 9

#mstepinit <-list(mustart,sigstart,lambdastart)

## EM from cattan paper

#source("EM_alg_v4cattan.R")
#EMest <- try(EM_alg(mstep.start = mstep.start, y = simdata2, nM = 2, conv = conv))
    #if ('try-error' %in% class(EMest)) next
    
## EM from mixtools package
    
EMest <- try(mvnormalmixEM(measures, lambda = lambdastart, mu = mustart, sigma = sigstart, epsilon =  conv, verb = TRUE, maxit = 1000)) 
if ('try-error' %in% class(EMest)) next
EMest[2:5] # These are the important parts of the EM output

EMstep[[i]] <- EMest[2:5] # save output from EM step
    
### Normal distribution
    
#mu_y <- mean(simdata2[,1])
#mu_x1 <- mean(simdata2[,2])
#mu_x2 <- mean(simdata2[,3])
#mus <- c(mu_y,mu_x1,mu_x2)
    
#sigma <- cov(simdata2)

################### 5.Use minimum distance to relate to joint distribution of latent factors

### Free parameters
# As I normalize one loading for each factor to 1, these are not free parameters to be estimated
# Also, we restrict each measurement to only affect one factor, leaving 2 loadings per factor to be estimated
# I denote all free parameters by 1, those normalized to 2 and all others by 0
# Note that these free parameters only refer to the factor loadings aspect of the problem
    
freeparam <- (cbind(rbind(2,1,1,0,0,0,0,0,0),rbind(0,0,0,2,1,1,0,0,0),rbind(0,0,0,0,0,0,2,1,1))) # nine measurements for three factors
    
# Starting values for min distance step (start all at 1)
# In order, these are for the factor loadings (only the free parameters), 
# the means of each mixture, variances in each mixture, variances of epsilon
    
param.start <- c(rep(1, length(which(freeparam==1))), rep(rep(1,nF),nM), rep(rep(1, .5*nF*(nF+1)), nM), rep(1, nMeas)) 

### Optimise the miniminum distance function
# note that here, we constrain the variances of epsilons to be positive (I have swapped the order relative to old code)
out <- optim(param.start,mindistance, method=c("L-BFGS-B"), 
                 lower=c(rep(-Inf, length(param.start)-nMeas), rep(0, nMeas)), upper=rep(Inf, length(param.start)))
    
# this produces the estimated parameters of the distribution of latent factors
    
### Organize the output so that each component is easily isolated

# Factor Loadings

lambdalong  <- rep(0, nMeas*nF) 
lambdalong[which(freeparam==1)]  <- out$par[1:length(which(freeparam==1))]
lambdalong[which(freeparam==2)]  <- 1
lambda_est  <- matrix(lambdalong, nMeas, nF) #this is now a 9x3 matrix of factor loadings, with the first normalized to one for each factor

# Means of nM mixtures

startMean <- length(which(freeparam==1))+1
endMean <- startMean + (nF*nM) -1
mean_est <- t(matrix(out$par[startMean:endMean], nF, nM)) # This is a matrix of means
    
# Variance-covariance matrices of the nM mixtures 
 
startCov <- endMean + 1
endCov <- startCov + length(rep(rep(1, .5*nF*(nF+1)), nM)) -1
    
Lcovfactor_est    <- list()
covfactor_est    <- list()
  for (m in 1:nM){
    Lcovfactor_est[[m]] <- matrix(0, nF, nF)
    Lcovfactor_est[[m]][lower.tri(Lcovfactor_est[[m]], diag=TRUE)] <- out$par[(startCov + (m-1)*0.5*nF*(nF+1)):(startCov + m*0.5*nF*(nF+1) - 1)] #have used lower.tri (Cattan) instead of LowerTri (Nix)
    covfactor_est[[m]] <- Lcovfactor_est[[m]] %*% t(Lcovfactor_est[[m]])
  }
    
# Variances of the uniquenesses 
    
starteps <- endCov +1
endeps <- length(param.start)
eps_est       <- out$par[starteps:endeps]

# Mixing parameter

prop_hat       <- EMest$lambda

est <- list(prob=prop_hat, eps=eps_est, cov=covfactor_est, mean=mean_est, lambda=lambda_est)
    
#print(est)
all <- list(out, est)
#return(all)
    
################### 6. Simulate from estimated distribution
# Now I use the results of the minimum distance step to simulate from the estimated distribution of factors

### Labeling estimates correctly
    
cov_estfin <- list()

for (m in 1:nM) {
  cov_estfin[[m]] <- est$cov[[m]]
  }
estsig <- cov_estfin
    
estmu <- list()
for (m in 1:nM) {
  estmu[[m]] <- est$mean[m,]
  }
    
estprops <- est$prob
    
### Mixture distribution
    
#To deal with symmetric issues
for (m in 1:nM) {
  estsig[[m]]<-forceSymmetric(estsig[[m]])
  }
#To deal with the positive definite issues
  for (m in 1:nM) {
  estsig[[m]]<-make.positive.definite(estsig[[m]])
  }
    
## Reduce to a single multivariate normal if estimated mixing parameter close to one
# also put variances and means into matrices for mixture draws
    
estsig <- rbind(estsig[[1]],estsig[[2]])
estmu <- rbind(estmu[[1]],estmu[[2]])
    
if (estprops[1]>0.99){
  simdata2 <- mvrnorm(1000, mu = estmu[1,], Sigma = estsig[1:3,1:3])
  }else
  if (estprops[2]>0.99){
  simdata2 <- mvrnorm(1000, mu = estmu[2,], Sigma = estsig[4:6,1:3])
  }else {
  simdata2 <- rmvnorm.mixt(1000, mus=estmu, Sigmas=estsig, props=estprops)
}
  
### Normal distribution
    
#simdata3 <- mvrnorm(1000, mu = mus, Sigma = sigma)

################### 7. Estimate CES production function
# Here I use NLS to estimate the production function over factors    

# need to make not logs again and put in dataframe form
simdata.3 <- data.frame(
        x1 = exp(simdata2[,1]), x2 = exp(simdata2[,2]), x3 = exp(simdata2[,3])
      )
    
if (rho>=0)
  {
  cesNlssim <-try(nlsLM(
      log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
      simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha = 0.7, beta = 0.3)
          )
      )
    } else if (rho<0 & rho>=-1) {cesNlssim <-try(nlsLM(
      log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
      simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = -0.3, alpha = 0.7, beta = 0.3)
          )
      )
    } else {
      cesNlssim <-try(nlsLM(
      log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
      simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = -1.5, alpha = 0.7, beta = 0.3)
          )
      )
    }
    if ('try-error' %in% class(cesNlssim)) next
    print(cesNlssim)
    estim <- summary(cesNlssim)
    
#remove extreme values 
    
if (estim[[10]][1]>=4){
    estim[[10]][1]<-NA
    }
    
#if (estim[[10]][1]<=-4){
  #  estim[[10]][1]<-NA
  #}
    
### Save estimates
    
rhoestimates[i,j] <- estim[[10]][1]
    
alphaestimates[i,j] <- estim[[10]][2]
    
betaestimates[i,j] <- estim[[10]][3]
    
    
  }
}

######################## 8. Organise output

### Gather and average EM outputs across runs

### Create empty vectors to fill with average parameter estimates, both pre and post EM/mindistance step

rhomeanest <- vector("numeric",length(rhos))
rhosdest <- vector("numeric",length(rhos))

alphameanest <- vector("numeric",length(rhos))
alphasdest <- vector("numeric",length(rhos))

betameanest <- vector("numeric",length(rhos))
betasdest <- vector("numeric",length(rhos))

rhomeanestpre <- vector("numeric",length(rhos))
rhosdestpre <- vector("numeric",length(rhos))

alphameanestpre <- vector("numeric",length(rhos))
alphasdestpre <- vector("numeric",length(rhos))

betameanestpre <- vector("numeric",length(rhos))
betasdestpre <- vector("numeric",length(rhos))

### Fill vectors with mean and standard deviation of individual estimates

for (i in 1:length(rhos))
{
  rhomeanest[i] <- mean(rhoestimates[i,], na.rm=TRUE)
  rhosdest[i] <- sd(rhoestimates[i,], na.rm=TRUE)
  
  alphameanest[i] <- mean(alphaestimates[i,], na.rm=TRUE)
  alphasdest[i] <- sd(alphaestimates[i,], na.rm=TRUE)
  
  betameanest[i] <- mean(betaestimates[i,], na.rm=TRUE)
  betasdest[i] <- sd(betaestimates[i,], na.rm=TRUE)
  
  rhomeanestpre[i] <- mean(rhoestimatespre[i,])
  rhosdestpre[i] <- sd(rhoestimatespre[i,])
  
  alphameanestpre[i] <- mean(alphaestimatespre[i,])
  alphasdestpre[i] <- sd(alphaestimatespre[i,])
  
  betameanestpre[i] <- mean(betaestimatespre[i,])
  betasdestpre[i] <- sd(betaestimatespre[i,])
}

### Combine into various lists

rhoest <- cbind(rhos, rhomeanest, rhosdest)

alphaest <- cbind(alpha, alphameanest, alphasdest)

betaest <- cbind(beta, betameanest, betasdest)

rhoestpre <- cbind(rhos, rhomeanestpre, rhosdestpre)

alphaestpre <- cbind(alpha, alphameanestpre, alphasdestpre)

betaestpre <- cbind(beta, betameanestpre, betasdestpre)

estpre <- list(rhoestpre, alphaestpre, betaestpre, rhoestimatespre, alphaestimatespre, betaestimatespre)

estimates <- list(rhoest, alphaest, betaest, rhoestimates, alphaestimates, betaestimates, estpre)

# Stop the clock
proc.time() - ptm


### Generate histograms and perform tests

#hist(log(y))
#hist(rhoestimates[1,])

#pearson.test(log(y))

### Save results

save(list = ls(all.names = TRUE), file = "Full simulations/Output/Spec_6 220216.RData")
