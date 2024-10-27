
# load data
load("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/analysis/flycatcher_Data.RData")
library(jagsUI)

####################################################################################################################
###
### run the model for the manuscript
###
####################################################################################################################
sink("pifl.jags")
cat("
    model {


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # individual latent variable
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

    for (k in 1:2){
      mu.p[k] ~ dnorm(0,0.1)
      sigma.p[k] ~ dunif(0,2)
      tau.p[k] <- pow(sigma.p[k], -2)

      for (t in 1:(n.years-1)){
        theta[t,k] ~ dnorm(mu.p[k], tau.p[k])
        logit(p[t,k]) <- theta[t,k]
    }}

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
    # age at first capture and latent quality
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
        log(h[i,t]) = beta.phi[1] +
                      beta.phi[2] * (a[i] - f[i] + t - 2) + 
                      gamma.phi * eta[i] + 
                      kappa[t]
        phi[i,t] = exp(-h[i,t])
        z[i,t+1] ~ dbern(z[i,t] * phi[i,t])
        y[i,t+1,1] ~ dbern(z[i,t+1] * p[t,s[i]])
      }
    }
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # fecundity
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    for (i in 1:nI){
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
jags.data <- list(y = y, s = s,
                  n.years = n.years, nI = dim(y)[1],
                  f = f, fec = fec, fec.limit = fec.limit,  
                  a = a, z = z)

parameters <- c('alpha.phi','beta.phi','gamma.phi',
                'alpha.xi','beta.xi',
                'p','pi','pi.star','rho','mu.p',
                'Sigma','kappa','omega', 'pi',
                'sigma.phi','sigma.lambda','sigma.p',
                'sigma.kappa','sigma.omega','sigma.eta')

inits <- function(){list(mu.p = rep(1,2))}  

nc <- 3
nt <- 10
ni <- 250000
nb <- 100000

Sys.time()
m.fl <- jags(jags.data, inits, parameters, "pifl.jags", parallel = T,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()



save.image("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/results/revision/250k_pifl_latent_linear.RData")


m.fl$Rhat
print(m.fl, digits = 3)


plot(m.fl$sims.list$beta.phi[,1], type = 'l')
plot(m.fl$sims.list$beta.phi[,2], type = 'l')
plot(m.fl$sims.list$alpha.phi, type = 'l')
plot(m.fl$sims.list$gamma.phi, type = 'l')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# preliminary plots
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###
### median phi and fecundity
###
par(mfrow = c(2,1), oma = rep(0,4), family = 'serif', mar = c(5.1,5.1,2.1,2.1))
plot(c(exp(-exp(m.fl$q50$alpha.phi)), exp(-exp(m.fl$q50$beta.phi[1] + m.fl$q50$beta.phi[2]*0:7))),
     ylim = c(0,1), type = 'l', lty = 1, cex.lab = 2,
     xlab = 'Female Age', ylab = 'Survival')

plot(exp(m.fl$q50$alpha.xi + m.fl$q50$beta.xi[1] * 0:8 + m.fl$q50$beta.xi[2] * 0:8 * 0:8), 
     ylim = c(0,1), type = 'l', lty = 1, cex.lab = 2,
     xlab = 'Female Age', ylab = 'Fecundity')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# partial plots
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
n.iter <- length(m.fl$sims.list$alpha.phi)
x <- seq(-0.5,0.5,length.out = res)

p.fl.phi <- array(NA, dim = c(n.iter, res))
p.fl.xi <- array(NA, dim = c(n.iter, res))
q.fl.phi <- matrix(NA, n.iter, 3)
q.fl.xi <- matrix(NA, n.iter, 3)
for (j in 1:res){
  p.fl.phi[,j] <- exp(-exp(m.fl$sims.list$alpha.phi + m.fl$sims.list$gamma.phi * x[j]))
  p.fl.xi[,j] <- exp(m.fl$sims.list$alpha.xi + 1 * x[j])
  q.fl.phi[j,] <- quantile(p.fl.phi[,j], c(0.025,0.5,0.975))
  q.fl.xi[j,] <- quantile(p.fl.xi[,j], c(0.025,0.5,0.975))  
}


par(mfrow = c(1,1), oma = rep(0,4), family = 'serif', mar = c(5.1,5.1,2.1,2.1))
smoothScatter(p.fl.phi ~ p.fl.xi, 
              nrpoints = 0, cex.lab = 2,
              xlab = expression(xi),
              ylab = expression(phi))
lines(q.fl.phi[,2] ~ q.fl.xi[,2], lwd = 2, col = 'white')
lines(q.fl.phi[,1] ~ q.fl.xi[,2], lwd = 2, lty = 2, col = 'white')
lines(q.fl.phi[,3] ~ q.fl.xi[,2], lwd = 2, lty = 2, col = 'white')



