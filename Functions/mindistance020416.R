mindistance      <-function(param) {
  ### This function chooses parameters of the distribution of factors such that the linking equations hold
  
  
  ### First, extract output from EM algorithm
  # These are our 'known' parameters for the distribution of measures.
  
  mean_hat  <- EMest$mu[[1]]
  for (m in 2:nM){
mean_hat      <- rbind(mean_hat,EMest$mu[[m]])
  }
  cov_hat       <- EMest$sigma

  ### Extract parameters initial values vector
  # This pulls the correct entries from the starting parameters provided
  
  # Factor loadings
  
  lambdalong  <- rep(0, nMeas*nF) # gives a long matrix of all possible factor loadings (normalized, freely estimated or other)
  lambdalong[which(freeparam==1)]  <- param[1:length(which(freeparam==1))] # maps the free parameters into lambdalong
  lambdalong[which(freeparam==2)]  <- 1 # puts in the normalized loadings
  lambda_par  <- matrix(lambdalong, nMeas, nF) #this is now a matrix of factor loadings, with the first normalized to one for each factor
  
  # Means of each mixture
  
  startMean <- length(which(freeparam==1))+1 # This is the first mean entry in the parameter matrix, as found using the 'freeparam' matrix of loadings
  endMean <- startMean + (nF*nM) -1
  meanfactor_par <- t(matrix(param[startMean:endMean], nF, nM)) # This is a matrix of estimated means
  
  # Variance-covariance matrices of the nM mixtures 
  
  startCov <- endMean + 1
  endCov <- startCov + length(rep(rep(1, .5*nF*(nF+1)), nM)) -1
  Lcovfactor_par    <- list()
  covfactor_par    <- list()
  for (m in 1:nM){
    Lcovfactor_par[[m]] <- matrix(0, nF, nF)
    Lcovfactor_par[[m]][lower.tri(Lcovfactor_par[[m]], diag=TRUE)] <- param[(startCov + (m-1)*0.5*nF*(nF+1)):(startCov + m*0.5*nF*(nF+1) - 1)] #have used lower.tri (Cattan) instead of LowerTri (Nix)
    covfactor_par[[m]] <- Lcovfactor_par[[m]] %*% t(Lcovfactor_par[[m]])
  }
  
  # Variances of the uniquenesses (adjust for variables measured without error by adding zero to end of vector)
  
  starteps <- endCov +1
  endeps <- length(param)
  eps_parvec  <- param[starteps:endeps]
  eps_noerror <- rep(0, nMeasnoerror)
  eps_parvec <- append(eps_parvec, eps_noerror)

  ### Implied measure parameters from estimated factor parameters
  
  # Implied means 
  
  mean_par      <- matrix(0, nM, nMeas)
  for (m in 1:(nM)){ 
    mean_par[m,] <- lambda_par %*% as.matrix(meanfactor_par[m,])
  }
  
  # Implied covariances
  
  eps_par  <-diag(eps_parvec)
  cov_par  <- list()
  for (m in 1:nM){
    cov_par[[m]] <-lambda_par%*%covfactor_par[[m]]%*%t(lambda_par) + eps_par
  } 
  
  ### Calculate sums of squares
  # _par = parameters derived from inputting estimated distribution of factors
  # _hat = parameters we can consider `true', ie derived from EM on measures data
  
  # Covariances 
  
  covsos     <- (sum((cov_hat[[1]]-cov_par[[1]])^2))
  for (m in 2:nM){
    covsos    <- covsos + sum((cov_hat[[m]]-cov_par[[m]])^2)
  }
  
  # Means 
  
  meansos  <- sum((mean_hat[1,]-mean_par[1,])^2)
  for (m in 2:nM){
    meansos    <- meansos + sum((mean_hat[m,]-mean_par[m,])^2)
  }
    
  # Sum of squares 
  
  
  
  sos <- meansos + (cov_weight) * covsos
  
  return(sos)
}