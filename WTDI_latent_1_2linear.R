


# load data
load("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/analysis/dipper_Data.RData")
library(jagsUI)

####################################################################################################################
###
### run the model for the manuscript
###
####################################################################################################################
sink("wtdi.jags")
cat("
    model {


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # individual random effects
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    for (i in 1:nI){
      eta.star[i] ~ dnorm(0, 1)
      eta[i] = eta.star[i] * sigma.eta
    }
    sigma.eta ~ dgamma(1,1)
    tau.eta = 1/(sigma.eta * sigma.eta)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # detection probability
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    mu.p ~ dnorm(2,1)
    sigma.p ~ dgamma(1,1)
    tau.p <- pow(sigma.p, -2)
    
    for (t in 1:(n.years-1)){
      theta[t] ~ dnorm(mu.p, tau.p)
      logit(p[t]) = theta[t] 
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # annual random effects
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    sigma.kappa ~ dgamma(1,1)
    tau.kappa <- pow(sigma.kappa, -2)
    sigma.omega ~ dgamma(1,1)
    tau.omega <- pow(sigma.omega, -2)
    
    for (t in 1:(n.years-1)){
      kappa.star[t] ~ dnorm(0, 1)     # survival random effect
      kappa[t] = kappa.star[t] * sigma.kappa
    }
    for (t in 1:n.years){
      omega.star[t] ~ dnorm(0, 1)     # fecundity random effect
      omega[t] = omega.star[t] * sigma.omega
    }


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # age priors
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    ###
    ### survival parameters
    ###
    alpha.phi ~ dlogis(0,1)
    beta.phi[1] ~ dlogis(0,1)
    beta.phi[2] ~ dnorm(0,0.01)
    alpha.xi ~ dnorm(0,0.1)    
    beta.xi[1] ~ dnorm(0,0.01)
    beta.xi[2] ~ dnorm(0,0.01)  
    gamma.phi ~ dnorm(0, 0.01)
    gamma.xi = 1  
    
    ###
    ### priors for unknown ages...
    ###
    for(j in 1:7){
      pi.star[j] ~ dgamma(1, 1)
      pi[j] <- pi.star[j]/sum(pi.star[1:7])
    }
 
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # age at first entry into the marked adult population
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:nI){
      a[i] ~ dcat(pi[1:7])
    }



    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # survival likelihood
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:nI){
      z[i,f[i]] <- 1
      for (t in f[i]:(n.years-1)){
        w[i,t] = step(a[i] - f[i] + t - 1.5)
        log(h[i,t]) = alpha.phi * (1 - w[i,t]) + 
                      beta.phi[1] * w[i,t] +
                      beta.phi[2] * w[i,t] * (a[i] - f[i] + t - 2) + 
                      gamma.phi * eta[i] + 
                      kappa[t]
        phi[i,t] = exp(-h[i,t])
        z[i,t+1] ~ dbern(z[i,t] * phi[i,t])
        y[i,t+1,1] ~ dbern(z[i,t+1] * p[t])
      }
    }
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # fecundity
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    for (i in n.fec){
      for (t in fec[i,1:fec.limit[i]]){
        log(xi[i,t]) = alpha.xi + beta.xi[1] * (a[i] - f[i] + t - 1) +
                           beta.xi[2] * (a[i] - f[i] + t - 1) * (a[i] - f[i] + t - 1) + 
                           gamma.xi * eta[i] + omega[t]
        y[i,t,2] ~ dpois(xi[i,t])
      }
    }
    





  }
  ",fill = TRUE)
sink()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# format data and inits, run model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
age.index <- c(1:6,rep(7,n.years+7))
jags.data <- list(y = y,
                  n.years = n.years, nI = dim(y)[1],
                  f = f, fec = fec, fec.limit = fec.limit,  
                  n.fec = n.fec, a = a, z = z)

parameters <- c('alpha.phi','beta.phi','gamma.phi',
                'alpha.xi','beta.xi',
                'p','pi','pi.star','rho','mu.p',
                'Sigma','kappa','omega', 'pi',
                'sigma.phi','sigma.lambda','sigma.p',
                'sigma.kappa','sigma.omega','sigma.eta')

inits <- function(){list(mu.p = 2)}  

nc <- 3
nt <- 10
ni <- 250000
nb <- 100000

Sys.time()
m.dp <- jags(jags.data, inits, parameters, "wtdi.jags", parallel = T,
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m.dp, digits = 3)


save.image("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/results/250k_wtdi_latent_nonlinear.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# preliminary plots
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###
### median phi and fecundity
###
par(mfrow = c(2,1), oma = rep(0,4), family = 'serif', mar = c(5.1,5.1,2.1,2.1))
plot(c(exp(-exp(m.dp$q50$alpha.phi)), exp(-exp(m.dp$q50$beta.phi[1] + m.dp$q50$beta.phi[2]*0:7))),
     ylim = c(0,1), type = 'l', lty = 1, cex.lab = 2,
     xlab = 'Female Age', ylab = 'Survival')

plot(exp(m.dp$q50$alpha.xi + m.dp$q50$beta.xi[1] * 0:8 + m.dp$q50$beta.xi[2] * 0:8 * 0:8), 
     ylim = c(0,1), type = 'l', lty = 1, cex.lab = 2,
     xlab = 'Female Age', ylab = 'Fecundity')

plot(m.dp$sims.list$sigma.eta)
