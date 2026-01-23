# DIASPARA WP2.2 LHT models - Viktor Thunell
# Yellow eel length at age model and mcmc
# Schnute growth equation

# model code
elaa.code <- nimbleCode({
  
  # length likelihood
  for(i in 1:nobs){
    
    length[i] ~ dnorm(l.mu[i], sd = l.sig)
    l.mu[i] <- ( L1[i]^P[i] + (L2[i]^P[i] - L1[i]^P[i]) * (1-exp(-K[i]*(age[i]-A1))) / (1-exp(-K[i]*(A2-A1))) )^(1/P[i])
    
    L1[i] <- exp(bpar[mb[i],1])
    P[i] <- bpar[mb[i],2]
    L2[i] <- exp(bpar[mb[i],3]) + bsL*sex[i] + btL*temp.sc[i]
    K[i] <- exp(bpar[mb[i],4]) + bsK*sex[i] + btK*temp.sc[i] + hK[hab[i]]
    
    # impute sex when missing in observations
    sex[i]~dbern(ps[1])
  }
  
  ## Priors
  l.sig ~ dexp(0.01)
  
  ps[1] ~ dbeta(1,1)  # prop of males (prob. to be 1)
  ps[2] <- 1-ps[1]
  
  for(j in 1:nmb){
    for(k in 1:4){
      bpar[j,k] ~ dnorm(par[(k-1)*ner+ermb[j]], sd = sig.par[(k-1)*ner+ermb[j]])  #ner=10  
    }
  }
  
  par[1:npar] ~ dmnorm(mu.par[1:npar], prec = tau_p[1:npar, 1:npar])     #npar=40, ecoreg x parameters
  
  # global pars (sd is uncertainty of the pars, assuming CV of 0.5)
  for(k in 1:ner){
    mu.par[k] ~ dnorm(mean = l1.lmean, sd = sqrt(log(1 + (0.5)^2)))
    mu.par[(ner+k)] ~ dnorm(mean = p.mean, sd = 0.5)
    mu.par[(ner*2+k)] ~ dnorm(mean = l2.lmean, sd = sqrt(log(1 + (0.5)^2)))
    mu.par[(ner*3+k)] ~ dnorm(mean = k.lmean, sd = sqrt(log(1 + (0.5)^2)))
  }
  
  #basin variation of pars
  tau_p[1:npar, 1:npar] <- inverse(sigma_p[1:npar, 1:npar])
  
  for(k in 1:npar){
    for(j in 1:npar){
      sigma_p[k,j] <- Rnew[k,j] * sig.par[k] * sig.par[j]
    }
  }
  
  # Prior for correlation matrix (LKJ prior)
  phi[1]  <- eta + (npar - 2)/2
  corY[1] ~ dbeta(phi[1], phi[1])
  r12   <- 2 * corY[1] - 1
  ##
  R[1,1]     <- 1
  R[1,2]     <- r12
  R[2,2]     <- sqrt(1 - r12^2)
  
  R[2:npar,1]   <- 0
  
  for (m in 2:(npar-1)) {
    ## Draw beta random variable
    phi[m] <- phi[(m-1)] - 0.5
    corY[m] ~ dbeta(m / 2, phi[m])
    ## Draw uniformly on a hypersphere
    for (jj in 1:m) {
      corZ[m, jj] ~ dnorm(0, 1)
    }
    scZ[m, 1:m] <- corZ[m, 1:m] / sqrt(inprod(corZ[m, 1:m], corZ[m, 1:m]))
    R[1:m,(m+1)] <- sqrt(corY[m]) * scZ[m,1:m]
    R[(m+1),(m+1)] <- sqrt(1 - corY[m])
    for(jk in (m+1):npar){
      R[jk,m] <- 0
    }
  }  #m
  
  Rnew[1:npar,1:npar] <- t(R[1:npar,1:npar]) %*% R[1:npar,1:npar]
  
  # basin variation of pars
  for(k in 1:ner){
    sig.par[k] ~ dexp(1/l1.lsd)
    sig.par[(ner+k)] ~ dexp(1/0.5)
    sig.par[(ner*2+k)] ~ dexp(1/l2.lsd)
    sig.par[(ner*3+k)] ~ dexp(1/k.lsd)
  }
  
  l1.lsd <- sqrt(log(1 + (0.5)^2)) # CV = 0.5
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

ner = length(unique(data.schnu$er))
ermb = data.schnu %>% distinct(main_bas, er) %>% arrange(main_bas)

const <- list(A1 = min(data.schnu$age),
              A2 = max(data.schnu$age),
              temp.sc = data.schnu$temp.sc, 
              hab = data.schnu$hab,
              ner = ner,
              npar = 4*ner,
              nobs = nrow(data.schnu),
              nmb = length(unique(data.schnu$main_bas)),
              nhab = length(unique(data.schnu$hab)),
              mb = data.schnu$main_bas,
              ermb = ermb$er,
              eta = 2
)

# initial values generating function
inits <- function(){
  list(l.sig = rexp(1,5),
       Ustar  = diag(const$npar),
       mu.par = c(rlnorm(const$ner,log(100), sd = 0.1),
                  rlnorm(const$ner,0,sd = 0.1),
                  rlnorm(const$ner,log(850),sd = 50),
                  rlnorm(const$ner,log(0.13), sd = 0.01)),
       bsL = rnorm(1,0,0.1),
       bsK = rnorm(1,0,0.1),
       btK = rnorm(1,0,0.1),
       btL = rnorm(1,0,0.1),
       hK = rnorm(const$nhab,0,0.01),
       ps = rnorm(2,0.5,.1),
       sig.par = c(rlnorm(const$ner,0.1,0.01),
                   rlnorm(const$ner,-1,0.01),
                   rlnorm(const$ner,1,0.1),
                   rlnorm(const$ner,-0.1,0.01)),
       corZ = matrix(rnorm((const$npar-1)*(const$npar-1), 0, 1),
                     nrow = const$npar-1, ncol = const$npar-1),
       corY = runif((const$npar-1), 0, 1)
  )}

# build model
elaa.model <- nimbleModel(elaa.code,
                             constants = const,
                             inits=inits(),
                             data = data.schnu %>% select(age,length,sex),
                             buildDerivs = TRUE) 

# compile model
elaa.c <- compileNimble(elaa.model)

# configure and build mcmc/hmc
monits = c("ps","l.sig","R","hK","mu.par","sig.par","btL","btK","sex","bsL","bsK","par","bpar","l.mu")
elaa.confmcmc <- configureHMC(elaa.c, monitors = monits,
                                 useConjugacy = FALSE, enableWAIC = TRUE)

elaa.mcmc <- buildMCMC(elaa.confmcmc, project = elaa.model)

# compile mcmc
elaa.mcmcc <- compileNimble(elaa.mcmc, project = elaa.model, resetFunctions = TRUE)

# HMC samples
elaa.samples.hmc <- runMCMC(elaa.mcmcc, niter = 20000, nburnin = 12500, nchains = 2, thin=5, WAIC=TRUE, samplesAsCodaMCMC = TRUE) 

saveRDS(elaa.samples, file = paste0(home,"/models_eel/samples/elaa.samples_",Sys.Date(),".RData"))
