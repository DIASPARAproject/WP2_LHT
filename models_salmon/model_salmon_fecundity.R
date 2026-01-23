# DIASPARA WP2.2 LHT models - Viktor Thunell
# Salmon fecundity model and mcmc

sfecu.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    
    n.eggs[i] ~ dnegbin(prob = p[i], size = r)
    p[i] <- r/(r+mu[i])
    log(mu[i]) <- a[su[i]] + b[su[i]] * length.lc[i] 
    
  }
  
  # priors
  r ~ dgamma(shape = 1, rate = 1)
  
  for(j in 1:nsu){
    a[j] ~ dnorm(mu.a, sd = sd.a)
    b[j] ~ dnorm(mu.b, sd = sd.b)
  }
  
  mu.a ~ dnorm(0, 1)
  sd.a ~ dexp(0.5)
  mu.b ~ dnorm(1, 1)
  sd.b ~ dexp(0.5)
  
})

nsu <- length(unique(data.fec$su))

#initial values generating function
inits <- function(){
  list(a = rnorm(nsu,0,1),
       mu.a = rnorm(,0,1),
       sd.a = rexp(1,0.5),
       b = rnorm(nsu,1,1),
       mu.b = rnorm(1,1,1),
       sd.b = rexp(1,0.5),
       r = rgamma(1,1,1))
}

# build model
sfecu.model <- nimbleModel(sfecu.code,
                           constants = list(nsu = nsu,
                                            nobs = nrow(data.sfecu),
                                            su = data.sfecu$su),
                           inits=inits(),
                           data = data.sfecu %>% select(n.eggs, length.lc),
                           buildDerivs = TRUE)

# configure hmc
sfecu.confhmc <- configureHMC(sfecu.model,
                            monitors = c("a","b","mu.a","mu.b","sd.a","sd.b","p","r","mu"),
                            enableWAIC = TRUE)

# build mcmc (use buidlHMC() when not using configureHMC())
sfecu.hmc <- buildMCMC(sfecu.confhmc)

# compile model
sfecu.c <- compileNimble(sfecu.model)

# compile mcmc  and specify the project model
sfecu.hmcc <- compileNimble(sfecu.hmc)

# hmc samples
sfecu.samples <- runMCMC(sfecu.hmcc, niter = 5000, nburnin = 3000, nchains = 3, WAIC=TRUE, samplesAsCodaMCMC = TRUE)

# Save samples
saveRDS(sfecu.samples, file = paste0(home,"/models_salmon/samples/sfecu.samples_",Sys.Date(),".RData"))
