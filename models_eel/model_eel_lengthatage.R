# Eel length at age model (based on schnute_v9)
# Schnute growth equation

# model code
elaa.code <- nimbleCode({
  
  # length likelihood
  for(i in 1:nobs){
    
    length[i] ~ dnorm(l.mu[i], sd = l.sig)
    l.mu[i] <- ( L1[i]^P[i] + (L2[i]^P[i] - L1[i]^P[i]) * (1-exp(-K[i]*(age[i]-A1))) / (1-exp(-K[i]*(A2-A1))) )^(1/P[i])
    
    L1[i] <- exp(par[mb[i],1])
    P[i] <- par[mb[i],2]
    L2[i] <- exp(par[mb[i],3]) + bsL*sex[i] + btL*temp.sc[i]
    K[i] <- exp(par[mb[i],4]) + bsK*sex[i] + btK*temp.sc[i] + hK[hab[i]]
    
    # impute sex when missing in observations
    sex[i]~dbern(ps[1])
    
  }
  
  ## Priors
  l.sig ~ dexp(0.01)
  
  ps[1] ~ dbeta(1,1)  # prop of males (prob. to be 1)
  ps[2] <- 1-ps[1]
  
  # LKJ prior for the par - correlation matrix
  Ustar[1:npar,1:npar] ~ dlkj_corr_cholesky(1.3, npar) # eta = 1.3
  U[1:npar,1:npar] <- uppertri_mult_diag(Ustar[1:npar, 1:npar], sig.par[1:npar])
  
  # basin specific mnorm prior on par
  for(i in 1:nmb){
    par[i, 1:npar] ~ dmnorm(mu.par[er[i],1:npar], cholesky = U[1:npar, 1:npar], prec_param = 0)
  }
  
  # global pars 
  for(k in 1:ner){
    mu.par[k,1] ~ dnorm(mean = l1.lmean, sd = sqrt(log(1 + (0.5)^2)))
    mu.par[k,2] ~ dnorm(mean = p.mean, sd = 0.5)
    mu.par[k,3] ~ dnorm(mean = l2.lmean, sd = sqrt(log(1 + (0.5)^2)))
    mu.par[k,4] ~ dnorm(mean = k.lmean, sd = sqrt(log(1 + (0.5)^2)))
  }
  
  # basin variation of pars
  sig.par[1] ~ dexp(1/l1.lsd)
  sig.par[2] ~ dexp(1/0.5)
  sig.par[3] ~ dexp(1/l2.lsd)
  sig.par[4] ~ dexp(1/k.lsd)
  
  l1.lsd <- sqrt(log(1 + (0.5)^2)) 
  l1.lmean <- log(100) - 0.5 * l1.lsd^2
  p.mean <- 0
  l2.lsd <- sqrt(log(1 + (200/850)^2)) # sd of the Linf in fishbase = 113
  l2.lmean <- log(850) - 0.5 * l2.lsd^2 # mean Linf from fishbase
  k.lsd <- sqrt(log(1 + (0.5)^2)) 
  k.lmean <- log(0.12) - 0.5 * k.lsd^2 # fishbase median vB k
  
  # slope sex
  bsL ~ dnorm(0, sd = 1)
  bsK ~ dnorm(0, sd = 1)
  
  for(j in 1:nhab){
    hK[j] ~ dnorm(0, sd = 0.1)
  }
  
  # slope temperature on L2 and K
  btL ~ dnorm(0, 1)
  btK ~ dnorm(0, 1)
  
})

# Function creating the Cholesky of the covar. matrix
uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(k in 1:p)
      out[ , k] <- mat[ , k] * vec[k]
    return(out)
    # turn off buildDerivs for the i index
  }, buildDerivs = list(run = list(ignore = c('k')))) 

const <- list(A1 = min(data.schnu$age),
              A2 = max(data.schnu$age),
              temp.sc = data.schnu$temp.sc, 
              hab = data.schnu$hab,
              npar = 4,
              ner = length(unique(data.schnu$er)),
              nobs = nrow(data.schnu),
              nmb = length(unique(data.schnu$main_bas)),
              nhab = length(unique(data.schnu$hab)),
              mb = data.schnu$main_bas,
              er = data.schnu$er
)

# initial values generating function
inits <- function(){
  list(l.sig = rexp(1,5),
       Ustar  = diag(const$npar),
       mu.par = matrix(c(rlnorm(const$ner,log(100), sd = 0.1),
                         rlnorm(const$ner,0,sd = 0.1),
                         rlnorm(const$ner,log(850),sd = 50),
                         rlnorm(const$ner,log(0.13), sd = 0.01)),ncol = 4),
       bsL = rnorm(1,0,0.1),
       bsK = rnorm(1,0,0.1),
       btK = rnorm(1,0,0.1),
       btL = rnorm(1,0,0.1),
       hK = rnorm(const$nhab,0,0.01),
       ps = rnorm(2,0.5,.1),
       sig.par = c(rlnorm(1,0.1,0.01),
                   rlnorm(1,-1,0.01),
                   rlnorm(1,1,0.1),
                   rlnorm(1,-0.1,0.01))
  )}

# build model
elaa.model <- nimbleModel(elaa.code,
                          constants = const,
                          inits=inits(),
                          data = data.schnu %>% select(age,length,sex),
                          buildDerivs = TRUE) 

elaa.model$simulate()
elaa.model$calculate()

# identify nodes to sample 
dataNodes <- elaa.model$getNodeNames(dataOnly = TRUE)
parentNodes <- elaa.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variabsLes...
stnodes <- elaa.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- elaa.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars to id NAs
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",elaa.model$calculate(mvars[i]) ))
}

# compile model
elaa.c <- compileNimble(elaa.model)

# configure and build mcmc/hmc
elaa.confmcmc <- configureHMC(elaa.c, monitors = 
                                c("ps","l.sig","Ustar","U","hK","mu.par","sig.par","btL","btK","sex","bsL","bsK","par","l.mu"),
                              useConjugacy = FALSE, enableWAIC = TRUE)

elaa.mcmc <- buildMCMC(elaa.confmcmc, project = elaa.model)

# compile mcmc
elaa.mcmcc <- compileNimble(elaa.mcmc, project = elaa.model, resetFunctions = TRUE)

# HMC samples
elaa.samples.hmc <- runMCMC(elaa.mcmcc, niter = 20000, nburnin = 12500, nchains = 2, thin=5, WAIC=TRUE, samplesAsCodaMCMC = TRUE) 
# NOTE: There are 115 individual pWAIC values that are greater than 0.4.

saveRDS(elaa.samples, file = paste0(home,"/data/elaa.samples_",Sys.Date(),".RData"))
