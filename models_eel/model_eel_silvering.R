# DIASPARA WP2.2 LHT models
# Eel length at silvering model and mcmc

esilv.code <- nimbleCode({
  
  # Likelihood
  for(i in 1:nobs){
    length[i] ~ dlnorm(meanlog = s.mu[i], sdlog = s.sig)
    s.mu[i] <- alpha[er[i],year[i]] + bs[er[i]]*sex[i] + bt*temp.sc[i]
  }
  
  s.sig ~ dexp(7)
  
  # first year
  for(l in 1:ner){
    alpha[l,1] ~ dnorm(mu.a.gl[l], sd = sd.a.gl[l]) 
    bs[l] ~ dnorm(mbsl, sd = sdbsl)
  }
  
  # transition of pars from year i to i+1
  for(i in 1:(nyear-1)){
    alpha[1:ner, (i+1)] ~ dmnorm(alpha[1:ner, i], prec = tau_p[1:ner,1:ner])
  }
  
  for(k in 1:ner){
    for(j in 1:ner){
      sigma_p[k,j] <- Rnew[k,j] * sd.a.gl[k] * sd.a.gl[j]  
    }
  }
  
  tau_p[1:ner, 1:ner]<- inverse(sigma_p[1:ner, 1:ner])
  
  #Prior for correlation matrix (LKJ prior)
  phi[1]  <- eta + (ner - 2)/2
  corY[1] ~ dbeta(phi[1], phi[1])
  r12   <- 2 * corY[1] - 1
  ##
  R[1,1]     <- 1
  R[1,2]     <- r12
  R[2,2]     <- sqrt(1 - r12^2)
  
  R[2:ner,1]   <- 0
  
  for (m in 2:(ner-1)) {
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
    for(jk in (m+1):ner){
      R[jk,m] <- 0
    }
  }  #m
  
  Rnew[1:ner,1:ner] <- t(R[1:ner,1:ner]) %*% R[1:ner,1:ner]
  
  for(j in 1:ner){
    mu.a.gl[j] ~ dnorm(msl,sdsl)
    sd.a.gl[j] ~ dgamma(shape = 1.5, rate = 1.5 / sdsl)
  }
  
  # Durif 2005 indiv. mean, https://doi.org/10.1111/j.0022-1112.2005.00662.x
  sdsl <- sqrt(log(1 + (124/658)^2))  
  msl <- log(658) - sdsl^2/2 
  sdbsl <- sqrt(log(1 + (23/393)^2))
  mbsl <- log(393/658) - sdbsl^2/2 #log ratio. Durif 2005, https://doi.org/10.1111/j.0022-1112.2005.00662.x
  
  # prior temp and age slopes
  bt ~ dnorm(0,1)
  
})

# initial values generating function for all nodes
inits <- function() {
  list(s.sig = rexp(1,5),
       alpha = matrix(rnorm(const$ner * const$nyear, log(658), 0.05),
                      nrow = const$ner, ncol = const$nyear),
       bs = rnorm(const$ner, log(393/658), 0.05),
       bt = rnorm(1, 0, 0.1),
       mu.a.gl  = rnorm(const$ner, log(658), 0.05),
       sd.a.gl  = runif(const$ner, 0.2, 0.5),
       corZ = matrix(rnorm((const$ner-1)*(const$ner-1), 0, 1),
                     nrow = const$ner-1, ncol = const$ner-1),
       corY = runif((const$ner-1), 0, 1)
  )
}

nobs <- nrow(data.silv)
ner <- length(unique(data.silv$er))
nyear <- length(unique(data.silv$year))

const <- list(nobs = nobs,
              eta = 2,
              ner = ner,
              nyear = nyear,
              temp.sc = data.silv$temp.sc,
              sex = data.silv$sex,
              year = data.silv$year,
              er = data.silv$er)

# build model
eelsilv.model <- nimbleModel(eelsilv.code,
                             constants = const,
                             inits=inits(),
                             data = data.esilv %>% select(length),buildDerivs = TRUE)

# compile model
silv.c <- compileNimble(silv.model)

monits = c("l.sig","mu.a.gl","sd.a.gl","alpha","bs","bt","s.mu","R")

# configure and build mcmc and add hmc to alpha and sigma nodes
silv.confmcmc <- configureHMC(silv.c, monitors = monits, enableWAIC = TRUE)

silv.mcmc <- buildMCMC(silv.confmcmc, project = silv.model)

# compile mcmc
silv.mcmcc <- compileNimble(silv.mcmc, project = silv.model)

silv.samples <- runMCMC(silv.mcmcc, niter = 6000, nburnin = 3000, nchains = 2, thin = 2, WAIC=TRUE, samplesAsCodaMCMC = TRUE)

saveRDS(silv.samples, file = paste0(home,"/data/silv.samples_",Sys.Date(),".RData"))
