########### MIXTURE SIMULATIONS EXPERIMENTS ###########################
#rm(list = ls())
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


########## Set up to loop over different values

rhos <- c(-1,-0.5,0,0.5,1)
#numsims <- 50
#estimates <- matrix(nrow = (length(rhos)*4), ncol = length(numsims))

# create empty vectors / matrices for storing estimates
rhoestimates <- matrix(nrow = length(rhos), ncol = numsims)
alphaestimates <-
  matrix(nrow = length(rhos), ncol = numsims)
betaestimates <- matrix(nrow = length(rhos), ncol = numsims)

# create empty vectors / matrices for storing estimates taken BEFORE simulation stage

rhoestimatespre <- matrix(nrow = length(rhos), ncol = numsims)
alphaestimatespre <-
  matrix(nrow = length(rhos), ncol = numsims)
betaestimatespre <- matrix(nrow = length(rhos), ncol = numsims)

# Create empty list for storing EM step estimates

EMstep <- list()

# Start the clock!
ptm <- proc.time()

iteration_num <- 0 

for (i in 1:length(rhos))
{
  
EMstep[[i]] <- list()

  for (j in 1:numsims)
  {
    
iteration_num <- iteration_num +1
print("iteration number")
print(iteration_num)

########### Simulate RHS variables

### Joint logmixnormal Xs

t <- 0.7
mus <- rbind(c(1,5), c(4,2)) # first two entries are first mixture
Sigmas <- rbind(cbind(c(1,-0.5), c(-0.5, 1)), cbind(c(0.3,-0.1), c(-0.1, 0.8)))
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

############## Generate LHS variable
    
x <-exp(logx) #take exponential to plug into true production function
    
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
    
############# Estimate CES function on 'true' simulated data
# This should be able to recover estimates of rho
    
    
simdata1 <- data.frame(x1 = y, x2 = x[,1], x3 = x[,2])
#plot(simdata1)
if (rho>=0)
  {
  cesNlstrue <-nlsLM(
    log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
    simdata1, control = nls.lm.control(maxiter = 1000), start = c(
    rho = 0.5, alpha = 0.7, beta = 0.3)
        )
  } else if (rho<0 & rho>=-1) {
  cesNlstrue <-nlsLM(
    log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
    simdata1, control = nls.lm.control(maxiter = 1000), start = c(
    rho = -0.3, alpha = 0.7, beta = 0.3)
        )
  } else {
  cesNlstrue <-nlsLM(
    log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
    simdata1, control = nls.lm.control(maxiter = 1000), start = c(
    rho = -1.5, alpha = 0.7, beta = 0.3)
        )
    }
    print(cesNlstrue)
    
    estimpre <- summary(cesNlstrue)
    rhoestimatespre[i,j] <- estimpre[[10]][1]
    alphaestimatespre[i,j] <- estimpre[[10]][2]
    betaestimatespre[i,j] <- estimpre[[10]][3]
    
############### Estimate joint distribution
    
simdata2 <- cbind(log(y),log(x))
simdata.2 <-data.frame(x1 = simdata2[,1], x2 = simdata2[,2], x3 = simdata2[,3])
#plot(simdata.2)
    
### Mixture distribution

# starting parameters

lambdastart <- c(0.3,0.3,0.3)

sigstart <- list(cov(simdata2), cov(simdata2), cov(simdata2))

kmeansplusplus <- kmeanspp(simdata2, k = 3)

mustart <- list()
mustart[[1]] <- kmeansplusplus[[2]][1,1:dim(simdata2)[[2]]]
mustart[[2]] <- kmeansplusplus[[2]][2,1:dim(simdata2)[[2]]]
mustart[[3]] <- kmeansplusplus[[2]][3,1:dim(simdata2)[[2]]]

# Number of mixtures              
nM <- 3
    
# Number of factors 
nF <-3
    
# Convergence criterion            
#conv <- 1e-08
    
## EM from mixtools package
    
EMest <-  try(mvnormalmixEM(simdata2, lambda = lambdastart, mu = mustart, sigma = sigstart, epsilon =  conv, verb = TRUE, maxit = 200))
if ('try-error' %in% class(EMest)) next
EMest[2:5]

EMstep[[i]][[j]] <- EMest[2:5] # save output from EM step

################### Simulate from estimated distribution
    
## Mixture distribution
    
#To deal with symmetric issues
for (m in 1:nM) {
  EMest[[4]][[m]]<-forceSymmetric(EMest[[4]][[m]])
  }
#To deal with the positive definite issues
  for (m in 1:nM) {
  EMest[[4]][[m]]<-make.positive.definite(EMest[[4]][[m]])
  }
    
estmu <- rbind(EMest[[3]][[1]],EMest[[3]][[2]],EMest[[3]][[3]])
estsig <- rbind(EMest[[4]][[1]],EMest[[4]][[2]],EMest[[4]][[3]])
estprops <- EMest[[2]]
    
# extra bit for case when reduces to a single mvn when only one mixture
    
if (estprops[1]>0.99){
  simdata3 <- mvrnorm(1000, mu = estmu[1,], Sigma = estsig[1:3,1:3])
  }else
  if (estprops[2]>0.99){
  simdata3 <- mvrnorm(1000, mu = estmu[2,], Sigma = estsig[4:6,1:3])
  }else {
  simdata3 <- rmvnorm.mixt(1000, mus=estmu, Sigmas=estsig, props=estprops)
  }

################### Estimate production function
    
# need to make not logs again
    
simdata.3 <-
    data.frame(
    x1 = exp(simdata3[,1]), x2 = exp(simdata3[,2]), x3 = exp(simdata3[,3])
    )
    

if (rho>=0)
{
  cesNlssim <-nlsLM(
    log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
    rho = 0.5, alpha = 0.7, beta = 0.3)
    )
  } else if (rho<0 & rho>=-1) {
    cesNlssim <-nlsLM(
    log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
    rho = -0.3, alpha = 0.7, beta = 0.3)
    )
  } else {
    cesNlssim <-nlsLM(
    log(x1) ~ ((1 / rho) * log(alpha * x2 ^ (rho) + beta * x3 ^ (rho))), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
    rho = -1.5, alpha = 0.7, beta = 0.3)
    )
  }
if ('try-error' %in% class(cesNlssim)) next
print(cesNlssim)
estim <- summary(cesNlssim)
    
#remove extreme values
    
if (estim[[10]][1]>=4){
  estim[[10]][1]<-NA
  }
    
if (estim[[10]][1]<=-4){
  estim[[10]][1]<-NA
  }

################### Save estimates
    
rhoestimates[i,j] <- estim[[10]][1]
alphaestimates[i,j] <- estim[[10]][2]
betaestimates[i,j] <- estim[[10]][3]
    
  }
}

# remove local optima

rhoestimates_wloc <- rhoestimates # create backup of full estimates before removing local optima

for (i in 1:length(rhos)){
  for (j in 1:numsims){
    if (rhos[i]!=0 & (rhoestimates[i,j])^2<=0.2){
      rhoestimates[i,j]<-NA
    }
  }
}

rhomeanest <- vector("numeric",length(rhos))
rhoMSE <- vector("numeric",length(rhos))
rhosdest <- vector("numeric",length(rhos))

alphameanest <- vector("numeric",length(rhos))
alphaMSE <- vector("numeric",length(rhos))
alphasdest <- vector("numeric",length(rhos))

betameanest <- vector("numeric",length(rhos))
betaMSE <- vector("numeric",length(rhos))
betasdest <- vector("numeric",length(rhos))

rhomeanestpre <- vector("numeric",length(rhos))
rhoMSEpre <- vector("numeric",length(rhos))
rhosdestpre <- vector("numeric",length(rhos))

alphameanestpre <- vector("numeric",length(rhos))
alphaMSEpre <- vector("numeric",length(rhos))
alphasdestpre <- vector("numeric",length(rhos))

betameanestpre <- vector("numeric",length(rhos))
betaMSEpre <- vector("numeric",length(rhos))
betasdestpre <- vector("numeric",length(rhos))

for (i in 1:length(rhos))
{
  rhomeanest[i] <- mean(rhoestimates[i,], na.rm=TRUE)
  rhoMSE[i] <- mean(((rhoestimates[i,]-rhos[i])^2), na.rm=TRUE)
  rhosdest[i] <- sd(rhoestimates[i,], na.rm=TRUE)
  
  alphameanest[i] <- mean(alphaestimates[i,], na.rm=TRUE)
  alphaMSE[i] <- mean(((alphaestimates[i,]-alpha)^2), na.rm=TRUE)
  alphasdest[i] <- sd(alphaestimates[i,], na.rm=TRUE)
  
  betameanest[i] <- mean(betaestimates[i,], na.rm=TRUE)
  betaMSE[i] <- mean(((betaestimates[i,]-beta)^2), na.rm=TRUE)
  betasdest[i] <- sd(betaestimates[i,], na.rm=TRUE)
  
  rhomeanestpre[i] <- mean(rhoestimatespre[i,])
  rhoMSEpre[i] <- mean(((rhoestimatespre[i,]-rhos[i])^2), na.rm=TRUE)
  rhosdestpre[i] <- sd(rhoestimatespre[i,])
  
  alphameanestpre[i] <- mean(alphaestimatespre[i,])
  alphaMSEpre[i] <- mean(((alphaestimatespre[i,]-alpha)^2), na.rm=TRUE)
  alphasdestpre[i] <- sd(alphaestimatespre[i,])
  
  betameanestpre[i] <- mean(betaestimatespre[i,])
  betaMSEpre[i] <- mean(((betaestimatespre[i,]-beta)^2), na.rm=TRUE)
  betasdestpre[i] <- sd(betaestimatespre[i,])
}



rhoest <- cbind(rhos, rhomeanest, rhoMSE, rhosdest)

alphaest <- cbind(alpha, alphameanest, alphaMSE, alphasdest)

betaest <- cbind(beta, betameanest, betaMSE, betasdest)

rhoestpre <- cbind(rhos, rhomeanestpre, rhoMSEpre, rhosdestpre)

alphaestpre <- cbind(alpha, alphameanestpre, alphaMSEpre, alphasdestpre)

betaestpre <- cbind(beta, betameanestpre, betaMSEpre, betasdestpre)

# Stop the clock
proc.time() - ptm

### Generate histograms and perform tests

#hist(log(y))
#hist(rhoestimates[1,])

#pearson.test(log(y))

#################### Export results
save(list = ls(all.names = TRUE), file = "Full simulations/No measures/Additional specifications/Spec_C 3mix 260216.RData")
