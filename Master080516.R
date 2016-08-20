############################################ MASTER ###########################
# This R script performs all the main tasks of the analysis with the China data
rm(list = ls())
set.seed(1625)

######################### LOAD PACKAGES ###############################

library("readstata13")
library("LICORS")
library("mixtools")
library("Matrix")
library("corpcor")
library("ks")
library("minpack.lm")
library("xtable")
library("texreg")
library("stargazer")
library("broom")

######################### LOAD USEFUL FUNCTIONS ###############################

source("Functions/mindistance020416.R") #minimum distance function

########################### SET UP BOOTSTRAP ##################################

ptm <- proc.time() # start the clock

bootstrap <- 1 # set = 1 if running bootstrap

#k <- 0 # set this = 0 if not running bootstrap, comment out if running bootstrap

numboot <- 500 # number of bootstraps
boot_est <- matrix(nrow = numboot, ncol = 108) # comment out when not bootstrapping

for (k in 1:numboot) { # comment out when not bootstrapping

########################## LOAD STATA .DTA FILE ###############################
  
if (bootstrap == 0 | k == 1) {
measures <- read.dta13("/Users/Jack/Dropbox/Documents/Oxford/Thesis/China/Output/measuresXXXXXX.dta", # measures file in .dta
                       convert.factors = F)
measuresorig <- measures
measuresfull <- measures
measuresfullmat <- as.matrix(measuresfull)
} else {
  nMeasures <- dim(measuresorig)[2]
  nObs <- dim(measuresorig)[1]
  measuresfull <- measuresorig[sample(nrow(measures),size=nObs,replace=TRUE),]
}

if (k == 1) {
  nMeasures <- dim(measuresorig)[2]
  nObs <- dim(measuresorig)[1]
  measuresfull <- measuresorig[sample(nrow(measuresorig),size=nObs,replace=TRUE),]
}

########################## ORGANISE DATA FOR ANALYSIS #########################

### Combine into relevant factors
cog_3 <- cbind(measuresfull$bayley_mdi3, measuresfull$bayley_pdi3, measuresfull$point3)
cog_2 <- cbind(measuresfull$bayley_mdi2, measuresfull$bayley_pdi2, measuresfull$point2)
cog_1 <- cbind(measuresfull$bayley_mdi1, measuresfull$bayley_pdi1, measuresfull$communicate1)
cog_0 <- cbind(measuresfull$bayley_mdi0, measuresfull$bayley_pdi0)
health_3 <- cbind(measuresfull$hb3, measuresfull$ch_weight3, measuresfull$ch_health_index3, measuresfull$ch_height3) 
health_2 <- cbind(measuresfull$hb2, measuresfull$ch_weight2, measuresfull$ch_health_index2, measuresfull$ch_height2) 
health_1 <- cbind(measuresfull$hb1, measuresfull$ch_weight1, measuresfull$ch_health_index1, measuresfull$ch_height1)
health_0 <- cbind(measuresfull$hb0, measuresfull$ch_weight0, measuresfull$ch_health_index0, measuresfull$ch_height0)
inv_2 <- cbind(measuresfull$numfood2, measuresfull$supplements2, measuresfull$inv_song2)
inv_1 <- cbind(measuresfull$numfood1, measuresfull$supplements1, measuresfull$inv_song1)
inv_0 <- cbind(measuresfull$supplements0, measuresfull$foriegnmilk) 
parhealth <- cbind(measuresfull$hb_mom0, measuresfull$hh_height_dad0) 
wealth <- cbind(measuresfull$durables3, measuresfull$hh_houseprice3, measuresfull$hh_size3, measuresfull$log_hh_income)
pared <- cbind(measuresfull$carer_edu0)

# make sure any variables we treat as observed without measurement error are at end of following lists
measureslist <- list(cog_3, cog_2, cog_1, cog_0, health_3, health_2, health_1, health_0, inv_2, inv_1, inv_0, parhealth, wealth, pared)
measures <- cbind(cog_3, cog_2, cog_1, cog_0, health_3, health_2, health_1, health_0, inv_2, inv_1, inv_0, parhealth, wealth, pared)

# Extract the number of measures per factor
nMeas_cog_3 <- dim(cog_3)[2]
nMeas_cog_2 <- dim(cog_2)[2]
nMeas_cog_1 <- dim(cog_1)[2]
nMeas_cog_0 <- dim(cog_0)[2]
nMeas_health_3 <- dim(health_3)[2]
nMeas_health_2 <- dim(health_2)[2]
nMeas_health_1 <- dim(health_1)[2]
nMeas_health_0 <- dim(health_0)[2]
nMeas_inv_2 <- dim(inv_2)[2]
nMeas_inv_1 <- dim(inv_1)[2]
nMeas_inv_0 <- dim(inv_0)[2]
nMeas_parhealth <- dim(parhealth)[2]
nMeas_wealth <- dim(wealth)[2]
nMeas_pared <- dim(pared)[2]

# Demean measures

measures_d <- measures
for (i in 1:dim(measures)[2]){
  measures_d[,i] <- measures[,i] - mean((measures)[,i])
}

# Normalize variance of each measure to one (since objective function very flat-bottomed otherwise)

for (i in 1:dim(measures)[2]){
  measures_d[,i] <- measures_d[,i] / sqrt(var((measures_d)[,i]))
}

####################### SET KEY PARAMETERS ###################################

# Number of mixtures              
nM <- 2

# Number of factors 
nF <- length(measureslist)

# Convergence criterion for EM alg         
conv <- 0.01

# Number of measurements
nMeas <- dim(measures_d)[2]

# Number of factors measured wihout error, ie those with only one measurements
nMeasnoerror <- 1 #pared

########################## ESTIMATE JOINT DISTRIBUTION OF MEASURES ############
# Here I use the EM algorithm to estimate the joint distribution of Measures   

# Set the starting parameters for the EM algorithm using kmeans++ algorithm

lambdastart <- c(0.5,0.5)
sigstart <- list(cov(measures_d), cov(measures_d))
kmeansplusplus <- try(kmeanspp(measures_d))
if ('try-error' %in% class(kmeansplusplus)) next

mustart <- list()
mustart[[1]] <- kmeansplusplus[[2]][1,1:dim(measures_d)[[2]]]
mustart[[2]] <- kmeansplusplus[[2]][2,1:dim(measures_d)[[2]]]

# EM from mixtools package

EMest <- try(mvnormalmixEM(measures_d, epsilon =  conv, verb = TRUE, maxit = 200, lambda = lambdastart, mu = mustart, sigma = sigstart))
if ('try-error' %in% class(EMest)) next # skip to next bootstrap replication if non-convergence

################### MINIMUM DISTANCE ######################################
# Here I apply the minimum distance step to relate the parameters of the distribution of measures
# to those from the distribution of factors

### Free parameters
# Here I build the matrix demonstrating which parameters in the factor loading matrix are free. 2 = normalised to unity, 1 = free.

# This is a useful function in building the free parameter matrix below
gen_parvec <- function(x){
x_parvec <- vector(mode = "numeric", length = x)
for (i in 1:x) {
  x_parvec[i] <- 1
}
x_parvec[1] <- 2
return(x_parvec)
}

# Build free parameter matrix
parvec <- c(gen_parvec(nMeas_cog_3), rep(0, nMeas), gen_parvec(nMeas_cog_2), rep(0, nMeas), 
            gen_parvec(nMeas_cog_1), rep(0, nMeas), gen_parvec(nMeas_cog_0), rep(0, nMeas),
            gen_parvec(nMeas_health_3), rep(0, nMeas), gen_parvec(nMeas_health_2), rep(0, nMeas),
            gen_parvec(nMeas_health_1), rep(0, nMeas),  gen_parvec(nMeas_health_0), rep(0, nMeas),
            gen_parvec(nMeas_inv_2), rep(0, nMeas), gen_parvec(nMeas_inv_1), rep(0, nMeas), 
            gen_parvec(nMeas_inv_0), rep(0, nMeas),
            gen_parvec(nMeas_parhealth), rep(0, nMeas), gen_parvec(nMeas_wealth),
            rep(0, nMeas), gen_parvec(nMeas_pared)) # add in included variables here
freeparam <- matrix(parvec, nrow = nMeas, ncol = nF)

### Starting values for minimum distance step
# This is a vector of starting values. In order, these are for the factor loadings (only the free parameters), 
# the means of each mixture, variances in each mixture, variances of epsilon

# Set starting epsilons to be 10% of the variance of each measure
varmeas <- apply(measures_d[,1:nMeas - nMeasnoerror], 2, var) # extract the variance
epsilon.start <- (varmeas / 10) # start the error variance at 0.1 of the variance of the measure

# Set full starting values vector
param.start <- c(rep(1, length(which(freeparam==1))), rep(rep(1,nF),nM), rep(rep(1, .5*nF*(nF+1)), nM), epsilon.start) 

### Optimise the miniminum distance function

cov_weight <- 1 # weight put on covariance terms in minimum distance step

out <- try(optim(param.start,mindistance, method=c("L-BFGS-B"), 
             lower=c(rep(-Inf, length(param.start)-(nMeas - nMeasnoerror)), rep(0, nMeas - nMeasnoerror)), 
             upper=rep(Inf, length(param.start))))
if ('try-error' %in% class(out)) next # skip if non-convergence

### Organize the output so that each component is easily isolated

# Factor Loadings
lambdalong  <- rep(0, nMeas*nF) 
lambdalong[which(freeparam==1)]  <- out$par[1:length(which(freeparam==1))]
lambdalong[which(freeparam==2)]  <- 1
lambda_est  <- matrix(lambdalong, nMeas, nF) #this is now a matrix of factor loadings, with the first normalized to one for each factor

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
eps_noerror <- rep(0, nMeasnoerror)
eps_est       <- append(out$par[starteps:endeps], eps_noerror)

# Mixing parameter
prop_hat       <- EMest$lambda

# Bring together all estimated parameters
est <- list(prob=prop_hat, eps=eps_est, cov=covfactor_est, mean=mean_est, lambda=lambda_est)


################### SIMULATE FACTORS #####################################
# Now I use the results of the minimum distance step to simulate from the estimated distribution of factors

### Relabel estimates from minimum distance step

estsig <- est$cov[[1]]
for (m in 2:nM){
  estsig <- rbind(estsig,est$cov[[m]])
}
estmu <- est$mean
estprops <- est$prob

### Draw from distribution
# Use single normal if heavily weighted towards one component

if (estprops[1]>0.99){
  simdata2 <- mvrnorm(1000, mu = estmu[1,], Sigma = estsig[1:3,1:3])
}else
  if (estprops[2]>0.99){
    simdata2 <- mvrnorm(1000, mu = estmu[2,], Sigma = estsig[4:6,1:3])
  }else {
    simdata2 <- rmvnorm.mixt(1000, mus=estmu, Sigmas=estsig, props=estprops) 
  }

### Demean factors

for (i in 1:dim(simdata2)[2]){
  simdata2[,i] <- simdata2[,i] - mean((simdata2)[,i])
}

### Calculate signal to noise ratio 

varfactors                   <- apply(simdata2, 2, var)
signal                   <- est$lambda^2 %*%  as.matrix(varfactors)
noise                    <- eps_est
sigtonoise                   <-  signal/(signal + noise)

### Exponentiate and put in dataframe form

simdata.3 <- data.frame(
  cog_3 = exp(simdata2[,1]), cog_2 = exp(simdata2[,2]), cog_1 = exp(simdata2[,3]), cog_0 = exp(simdata2[,4]),
  health_3 = exp(simdata2[,5]), health_2 = exp(simdata2[,6]), health_1 = exp(simdata2[,7]), health_0 = exp(simdata2[,8]),
  inv_2 = exp(simdata2[,9]), inv_1 = exp(simdata2[,10]), inv_0 = exp(simdata2[,11]),
  parhealth = exp(simdata2[,12]), wealth = exp(simdata2[,13]), pared = exp(simdata2[,14])
)

################### ESTIMATE PRODUCTION FUNCTION IGNORING ENDOGENEITY #################################
# Note that these do not always converge

### Cog

PFcog3 <- try(nlsLM(
  log(cog_3) ~ ((1 / rho) * log(
    alpha1 * cog_2 ^ (rho) + alpha2 * health_2 ^ (rho)
    + alpha3 * inv_2 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1
    )
))
if ('try-error' %in% class(PFcog3)) { 
  next;
} else {
summary(PFcog3) }

PFcog2 <- try(nlsLM(
  log(cog_2) ~ ((1 / rho) * log(
    alpha1 * cog_1 ^ (rho) + alpha2 * health_1 ^ (rho)
    + alpha3 * inv_1 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1
    )
))
if ('try-error' %in% class(PFcog2)) { 
  next;
} else {
  summary(PFcog2) }

PFcog1 <- try(nlsLM(
  log(cog_1) ~ ((1 / rho) * log(
    alpha1 * cog_0 ^ (rho) + alpha2 * health_0 ^ (rho)
    + alpha3 * inv_0 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1
    )
))
if ('try-error' %in% class(PFcog1)) { 
  next;
} else {
  summary(PFcog1) }

### Health

PFhealth3 <- try(nlsLM(
  log(health_3) ~ ((1 / rho) * log(
    alpha1 * cog_2 ^ (rho) + alpha2 * health_2 ^ (rho)
    + alpha3 * inv_2 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1
    )
))
if ('try-error' %in% class(PFhealth3)) { 
  next;
} else {
  summary(PFhealth3) }

PFhealth2 <- try(nlsLM(
  log(health_2) ~ ((1 / rho) * log(
    alpha1 * cog_1 ^ (rho) + alpha2 * health_1 ^ (rho)
    + alpha3 * inv_1 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1
    )
))
if ('try-error' %in% class(PFhealth2)) { 
  next;
} else {
  summary(PFhealth2) }

PFhealth1 <- try(nlsLM(
  log(health_1) ~ ((1 / rho) * log(
    alpha1 * cog_0 ^ (rho) + alpha2 * health_0 ^ (rho)
    + alpha3 * inv_0 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1
    )
))
if ('try-error' %in% class(PFhealth1)) { 
  next;
} else {
  summary(PFhealth1) }

### Cog non-markov

#PFcog3nonmarkov <- try(nlsLM(
#  log(cog_3) ~ ((1 / rho) * log(
#    alpha1 * cog_2 ^ (rho) + alpha2 * cog_1 ^ (rho) + alpha3 * cog_0 ^ (rho)
#    + alpha4 * health_2 ^ (rho) + alpha5 * health_1 ^ (rho) + alpha6 * health_0 ^ (rho)
#    + alpha7 * inv_2 ^ (rho) + alpha8 * parhealth ^ (rho)
#    + alpha9 * pared ^ (rho)
#  )
#  ), data =
#    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
#      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, alpha5 = 0.2, alpha6 = 0.2, alpha7 = 0.2,
#      alpha8 = 0.2, alpha9 = 0.2
#    )
#))
#summary(PFcog3nonmarkov)

#PFcog2nonmarkov <- try(nlsLM(
#  log(cog_2) ~ ((1 / rho) * log(
#    alpha1 * cog_1 ^ (rho) + alpha2 * cog_0 ^ (rho) + alpha3 * health_1 ^ (rho)
#    + alpha4 * health_0 ^ (rho) + alpha5 * inv_1 ^ (rho) + alpha6 * parhealth ^ (rho)
#    + alpha7 * pared ^ (rho)
#  )
#  ), data =
#    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
#      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, alpha5 = 0.2, alpha6 = 0.2, alpha7 = 0.2
#    )
#))
#summary(PFcog2nonmarkov)

### Health non-markov

#PFhealth3nonmarkov <- try(nlsLM(
#  log(health_3) ~ ((1 / rho) * log(
#    alpha1 * cog_2 ^ (rho) + alpha2 * cog_1 ^ (rho) + alpha3 * cog_0 ^ (rho)
#    + alpha4 * health_2 ^ (rho) + alpha5 * health_1 ^ (rho) + alpha6 * health_0 ^ (rho)
#    + alpha7 * inv_2 ^ (rho) + alpha8 * parhealth ^ (rho)
#    + alpha9 * pared ^ (rho)
#  )
#  ), data =
#    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
#      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, alpha5 = 0.2, alpha6 = 0.2, alpha7 = 0.2,
#      alpha8 = 0.2, alpha9 = 0.2
#    )
#))
#summary(PFhealth3nonmarkov)

#PFhealth2nonmarkov <- try(nlsLM(
#  log(health_2) ~ ((1 / rho) * log(
#    alpha1 * cog_1 ^ (rho) + alpha2 * cog_0 ^ (rho) + alpha3 * health_1 ^ (rho)
#    + alpha4 * health_0 ^ (rho) + alpha5 * inv_1 ^ (rho) + alpha6 * parhealth ^ (rho)
#    + alpha7 * pared ^ (rho)
#  )
#  ), data =
#    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
#      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, alpha5 = 0.2, alpha6 = 0.2, alpha7 = 0.2
#    )
#))
#summary(PFhealth2nonmarkov)


################### ESTIMATE PRODUCTION FUNCTION ACCOUNTING FOR ENDOGENEITY #############

### Investment equations

inv_eq0 <-
  lm(
    log(inv_0) ~ log(cog_0) + log(health_0) + log(parhealth) + log(pared) + log(wealth), data = simdata.3
  )
summary(inv_eq0)
resid0 <- inv_eq0$residuals

inv_eq1 <-
  lm(
    log(inv_1) ~ log(cog_1) + log(health_1) + log(parhealth) + log(pared) + log(wealth), data = simdata.3
  )
summary(inv_eq1)
resid1 <- inv_eq1$residuals

inv_eq2 <-
  lm(
    log(inv_2) ~ log(cog_2) + log(health_2) + log(parhealth) + log(pared) + log(wealth), data = simdata.3
  )
summary(inv_eq2)
resid2 <- inv_eq2$residuals

#inv_eq0nonmarkov <-
#  lm(
#    log(inv_0) ~ log(cog_0) + log(health_0) + log(parhealth) + log(pared) + log(wealth), data = simdata.3
#  )
#summary(inv_eq0nonmarkov)
#resid0nonmarkov <- inv_eq0nonmarkov$residuals

#inv_eq1nonmarkov <-
#  lm(
#    log(inv_1) ~ log(cog_1) + log(health_1) +  log(cog_0) + log(health_0) +
#      log(parhealth) + log(pared) + log(wealth), data = simdata.3
#  )
#summary(inv_eq1nonmarkov)
#resid1nonmarkov <- inv_eq1nonmarkov$residuals

#inv_eq2nonmarkov <-
#  lm(
#    log(inv_2) ~ log(cog_2) + log(health_2) + log(cog_1) + log(health_1) +  log(cog_0) + log(health_0) +
#      log(parhealth) + log(pared) + log(wealth), data = simdata.3
#  )
#summary(inv_eq2nonmarkov)
#resid2nonmarkov <- inv_eq2nonmarkov$residuals


### Production Functions

### Cog

PFcog1CF <- try(nlsLM(
  log(cog_1) ~ ((1 / rho) * log(
    alpha1 * cog_0 ^ (rho) + alpha2 * health_0 ^ (rho)
    + alpha3 * inv_0 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A + beta * resid0
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1, beta = 0.1
    )
))
if ('try-error' %in% class(PFcog1CF)) { 
  next;
} else {
  summary(PFcog1CF) }

PFcog2CF <- try(nlsLM(
  log(cog_2) ~ ((1 / rho) * log(
    alpha1 * cog_1 ^ (rho) + alpha2 * health_1 ^ (rho)
    + alpha3 * inv_1 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A + beta * resid1
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1, beta = 0.1
    )
))
if ('try-error' %in% class(PFcog2CF)) { 
  next;
} else {
  summary(PFcog2CF) }

PFcog3CF <- try(nlsLM(
  log(cog_3) ~ ((1 / rho) * log(
    alpha1 * cog_2 ^ (rho) + alpha2 * health_2 ^ (rho)
    + alpha3 * inv_2 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A + beta * resid2
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1, beta = 0.1
    )
))
if ('try-error' %in% class(PFcog3CF)) { 
  next;
} else {
  summary(PFcog3CF) }

### Health

PFhealth1CF <- try(nlsLM(
  log(health_1) ~ ((1 / rho) * log(
    alpha1 * cog_0 ^ (rho) + alpha2 * health_0 ^ (rho)
    + alpha3 * inv_0 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A + beta * resid0
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1, beta = 0.1
    )
))
if ('try-error' %in% class(PFhealth1CF)) { 
  next;
} else {
  summary(PFhealth1CF) }

PFhealth2CF <- try(nlsLM(
  log(health_2) ~ ((1 / rho) * log(
    alpha1 * cog_1 ^ (rho) + alpha2 * health_1 ^ (rho)
    + alpha3 * inv_1 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A + beta * resid1
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1, beta = 0.1
    )
))
if ('try-error' %in% class(PFhealth2CF)) { 
  next;
} else {
  summary(PFhealth2CF) }

PFhealth3CF <- try(nlsLM(
  log(health_3) ~ ((1 / rho) * log(
    alpha1 * cog_2 ^ (rho) + alpha2 * health_2 ^ (rho)
    + alpha3 * inv_2 ^ (rho) + alpha4 * parhealth ^ (rho)
    + (1 - alpha1 - alpha2 - alpha3 - alpha4) * pared ^ (rho)
  ) + A + beta * resid2
  ), data =
    simdata.3, control = nls.lm.control(maxiter = 1000), start = c(
      rho = 0.5, alpha1 = 0.4, alpha2 = 0.3, alpha3 = 0.3, alpha4 = 0.2, A = 1, beta = 0.1
    )
))
if ('try-error' %in% class(PFhealth3CF)) { 
  next;
} else {
  summary(PFhealth3CF) }

###################### EXPORT RESULTS ##########################################

### Collect estimated parameters into vectors

## Production Functions

PFcog3_coeff <- vector()
for (i in 1:6) {
PFcog3_coeff <- c(PFcog3_coeff,summary(PFcog3)$coefficients[i,1])
}
PFcog3_coeff <- c(PFcog3_coeff,summary(PFcog3)$coefficients[6,1]) # Need to get alpha5 from other estimates
PFcog3_coeff[6] <- 1 - sum(PFcog3_coeff[2:5])

PFcog2_coeff <- vector()
for (i in 1:6) {
  PFcog2_coeff <- c(PFcog2_coeff,summary(PFcog2)$coefficients[i,1])
}
PFcog2_coeff <- c(PFcog2_coeff,summary(PFcog2)$coefficients[6,1]) # Need to get alpha5 from other estimates
PFcog2_coeff[6] <- 1 - sum(PFcog2_coeff[2:5])

PFcog1_coeff <- vector()
for (i in 1:6) {
  PFcog1_coeff <- c(PFcog1_coeff,summary(PFcog1)$coefficients[i,1])
}
PFcog1_coeff <- c(PFcog1_coeff,summary(PFcog1)$coefficients[6,1]) # Need to get alpha5 from other estimates
PFcog1_coeff[6] <- 1 - sum(PFcog1_coeff[2:5])

PFhealth3_coeff <- vector()
for (i in 1:6) {
  PFhealth3_coeff <- c(PFhealth3_coeff,summary(PFhealth3)$coefficients[i,1])
}
PFhealth3_coeff <- c(PFhealth3_coeff,summary(PFhealth3)$coefficients[6,1]) # Need to get alpha5 from other estimates
PFhealth3_coeff[6] <- 1 - sum(PFhealth3_coeff[2:5])

PFhealth2_coeff <- vector()
for (i in 1:6) {
  PFhealth2_coeff <- c(PFhealth2_coeff,summary(PFhealth2)$coefficients[i,1])
}
PFhealth2_coeff <- c(PFhealth2_coeff,summary(PFhealth2)$coefficients[6,1]) # Need to get alpha5 from other estimates
PFhealth2_coeff[6] <- 1 - sum(PFhealth2_coeff[2:5])

PFhealth1_coeff <- vector()
for (i in 1:6) {
  PFhealth1_coeff <- c(PFhealth1_coeff,summary(PFhealth1)$coefficients[i,1])
}
PFhealth1_coeff <- c(PFhealth1_coeff,summary(PFhealth1)$coefficients[6,1]) # Need to get alpha5 from other estimates
PFhealth1_coeff[6] <- 1 - sum(PFhealth1_coeff[2:5])

PFcog3CF_coeff <- vector()
for (i in 1:7) {
  PFcog3CF_coeff <- c(PFcog3CF_coeff,summary(PFcog3CF)$coefficients[i,1])
}
PFcog3CF_coeff <- c(PFcog3CF_coeff,summary(PFcog3CF)$coefficients[7,1]) # Need to get alpha5 from other estimates
PFcog3CF_coeff[7] <- PFcog3CF_coeff[6]
PFcog3CF_coeff[6] <- 1 - sum(PFcog3CF_coeff[2:5])

PFcog2CF_coeff <- vector()
for (i in 1:7) {
  PFcog2CF_coeff <- c(PFcog2CF_coeff,summary(PFcog2CF)$coefficients[i,1])
}
PFcog2CF_coeff <- c(PFcog2CF_coeff,summary(PFcog2CF)$coefficients[7,1]) # Need to get alpha5 from other estimates
PFcog2CF_coeff[7] <- PFcog2CF_coeff[6]
PFcog2CF_coeff[6] <- 1 - sum(PFcog2CF_coeff[2:5])

PFcog1CF_coeff <- vector()
for (i in 1:7) {
  PFcog1CF_coeff <- c(PFcog1CF_coeff,summary(PFcog1CF)$coefficients[i,1])
}
PFcog1CF_coeff <- c(PFcog1CF_coeff,summary(PFcog1CF)$coefficients[7,1]) # Need to get alpha5 from other estimates
PFcog1CF_coeff[7] <- PFcog1CF_coeff[6]
PFcog1CF_coeff[6] <- 1 - sum(PFcog1CF_coeff[2:5])

PFhealth3CF_coeff <- vector()
for (i in 1:7) {
  PFhealth3CF_coeff <- c(PFhealth3CF_coeff,summary(PFhealth3CF)$coefficients[i,1])
}
PFhealth3CF_coeff <- c(PFhealth3CF_coeff,summary(PFhealth3CF)$coefficients[7,1]) # Need to get alpha5 from other estimates
PFhealth3CF_coeff[7] <- PFhealth3CF_coeff[6]
PFhealth3CF_coeff[6] <- 1 - sum(PFhealth3CF_coeff[2:5])

PFhealth2CF_coeff <- vector()
for (i in 1:7) {
  PFhealth2CF_coeff <- c(PFhealth2CF_coeff,summary(PFhealth2CF)$coefficients[i,1])
}
PFhealth2CF_coeff <- c(PFhealth2CF_coeff,summary(PFhealth2CF)$coefficients[7,1]) # Need to get alpha5 from other estimates
PFhealth2CF_coeff[7] <- PFhealth2CF_coeff[6]
PFhealth2CF_coeff[6] <- 1 - sum(PFhealth2CF_coeff[2:5])

PFhealth1CF_coeff <- vector()
for (i in 1:7) {
  PFhealth1CF_coeff <- c(PFhealth1CF_coeff,summary(PFhealth1CF)$coefficients[i,1])
}
PFhealth1CF_coeff <- c(PFhealth1CF_coeff,summary(PFhealth1CF)$coefficients[7,1]) # Need to get alpha5 from other estimates
PFhealth1CF_coeff[7] <- PFhealth1CF_coeff[6]
PFhealth1CF_coeff[6] <- 1 - sum(PFhealth1CF_coeff[2:5])

## Investment equations

inv_eq2_coeff <- vector()
for (i in 1:6) {
  inv_eq2_coeff <- c(inv_eq2_coeff,summary(inv_eq2)$coefficients[i,1])
}

inv_eq1_coeff <- vector()
for (i in 1:6) {
  inv_eq1_coeff <- c(inv_eq1_coeff,summary(inv_eq1)$coefficients[i,1])
}

inv_eq0_coeff <- vector()
for (i in 1:6) {
  inv_eq0_coeff <- c(inv_eq0_coeff,summary(inv_eq0)$coefficients[i,1])
}

# Combine into list of coefficients
coef <- c(PFcog3_coeff, PFcog2_coeff, PFcog1_coeff, 
                  PFhealth3_coeff, PFhealth2_coeff, PFhealth1_coeff,
                  PFcog3CF_coeff, PFcog2CF_coeff, PFcog1CF_coeff, 
                  PFhealth3CF_coeff, PFhealth2CF_coeff, PFhealth1CF_coeff,
                  inv_eq2_coeff, inv_eq1_coeff, inv_eq0_coeff)

### Save environment

if (bootstrap == 0) {
coef_est <- coef # This is the set of estimated parameters, to be combined with bootstrap estimates in tables
save(list = ls(all.names = TRUE), file = "Output/masteroutXXXXXX") # This should contain all non-bootstrap estimates
} else {
boot_est[k,] <- coef
}

print ("bootstrap number") #comment these lines out if not bootstrapping
print(k)
} # comment out when not bootstrapping

if (bootstrap == 1) {
CI <- matrix(nrow = 2, ncol = length(coef))
SE <- matrix(nrow = 1, ncol = length(coef))
for (i in 1:length(coef)) {
CI[,i] <- quantile(boot_est[,i], c(.05, .95), na.rm = TRUE) # Calculate confidence intervals
SE[i] <- sd(boot_est[,i], na.rm = TRUE) # Calculate standard error
}
boot <- list(CI,SE,boot_est)
}

if (bootstrap == 1) {
  save(boot, file = "Bootstrap/bootXXXXXX.RData")
}

time <- proc.time() - ptm # stop the clock
