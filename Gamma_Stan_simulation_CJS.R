# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(MASS)
library(LaplacesDemon)
library(MCMCvis)
library(parallel)
library(rstan)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of chains/cores, and number of simulations
# number of iterations, length of burn-in, thinning rate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nc <- 4
n.sims <- 1000
ni <- 10000
nb <- 5000
nt <- 5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of releases and occasions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Releases = c(25)                            # releases
Occasions = c(22)                           # intermediate/long-term study

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# arrays to store simulation parameters and parameter estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
truth <- array(NA, dim = c(n.sims, 7))
results.mvn <- array(NA, dim = c(n.sims, 7, 10))
results.sem <- array(NA, dim = c(n.sims, 6, 10))
mvn.run <- rep(NA, n.sims)
sem.run <- rep(NA, n.sims)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stan models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cjs_SEM <- stan_model(file = 'C:/Users/thomas.riecke/Desktop/Simulation_runs/SEM_hazard/cjs_SEM_gamma.stan')
cjs_MVN <- stan_model(file = 'C:/Users/thomas.riecke/Desktop/Simulation_runs/SEM_hazard/cjs_MVN_gamma.stan')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# begin simulation for loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (ii in 1:n.sims){

      
      print(paste0('Simulation ', ii,', Releases ',
                   Releases,' Occasions ',Occasions))
      print(Sys.time())
      n.years <- Occasions
      
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # intercepts (constants)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mu <- rep(NA, 2)
      mu[1] <- -0.5
      mu[2] <- -0.5
      # exp(-exp(mu[1])) + exp(mu[2])
      
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # sd and correlation (constants)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sig <- c(0.4,0.4)
      rho <- -0.6
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # sd and correlation (sampling a range)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # sig <- NULL
      # sig[1] <- rgamma(1, 5, 10)
      # sig[2] <- rgamma(1, 5, 10)
      # rho <- runif(1, 0.25, 0.9)
      
      # variance-covariance matrix
      Sigma <- matrix(c(sig[1] * sig[1],        sig[1] * sig[2] * rho,
                        sig[2] * sig[1] * rho,  sig[2] * sig[2]),
                      byrow = T, nrow = 2)
      
      # number of releases per year, total number of individuals, first release vector
      rel = rep(Releases, n.years - 2)            # releases per year
      n.ind = sum(rel)                            # number of individuals
      f = rep(seq(1, n.years - 2), times = rel)   # first release
      
      
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # individual survival and fecundity and detection
      # set phi and gamma as matrices and add random effects to induce 
      # temporal variation in demographic parameters
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      epsilon = mvrnorm(n.ind, c(0,0), Sigma)  # heterogeneity parameters
      
      h = exp(mu[1] + epsilon[,1])
      phi = exp(-h)        # survival probability
      gamma = exp(mu[2] + epsilon[,2])         # expected fecundity
 
      plot(phi ~ gamma)
      # hist(rbeta(10000,85,15), main = NULL, breaks = 500)
      p <- 0.8
      
      # fill 'truth' array with true values used to simulate data
      truth[ii,] <- c(mu, sig, rho, p, Sigma[1,2])
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # simulate latent states, capture histories, and fecundity
      # z: latent state (1 = alive; 2 = dead)
      # y: 3d array to store capture history [,,1] and fecundity data [,,2]
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      z <- matrix(NA, length(f), n.years)
      y <- array(NA, dim = c(length(f), n.years, 2))  
      
      for (i in 1:n.ind){
        
        z[i,f[i]] = 1
        y[i,f[i],1] = 1
        
        for (t in (f[i]+1):n.years){
          z[i,t] <- rbinom(1, z[i,t-1], phi[i])
          y[i,t,1] <- rbinom(1, z[i,t], p)
        }
        
        # sample fecundity if observed during breeding season
        for (t in f[i]:n.years){
          if (y[i,t,1] == 1 & !is.na(y[i,t,1])){
            y[i,t,2] <- rpois(1, gamma[i])
          }
        }
        
      }

      get.l <- function(x){max(which(x == 1))}
      l <- apply(y[,,1], 1, get.l)
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # model input
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      table(y[,,1])      
      table(y[,,2])
      # y[is.na(y)] <- -999
      # y[,,1][y[,,1] == 0] <- 2
      y[,,1][is.na(y[,,1])] <- -9
      y[,,2][is.na(y[,,2])] <- -9
      stan_data <- list(y = y,
                        nI = n.ind,
                        nT = n.years,
                        f = f,
                        l = l,
                        Zero = rep(0,2))
      

      
      # dim(y)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # run Stan in parallel: multivariate normal (m1) and SEM (m2)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

      # Sys.time()
      print('MVN')
      print(Sys.time())
      s1 <- Sys.time()
      m1 <- sampling(cjs_MVN, 
                     data = stan_data,
                     chains = nc, iter = ni, warmup = nb, thin = nt,
                     pars = c('rho','sigma','mu','p','cov'),
                     open_progress = FALSE, refresh = 1000, cores = nc)
      e1 <- Sys.time()
      mvn.run[ii] <- e1 - s1
      # Sys.time()
      
      
      # Sys.time()
      print('SEM')
      print(Sys.time())
      s2 <- Sys.time()
      m2 <- sampling(cjs_SEM, 
                     data = stan_data,
                     chains = nc, iter = ni, warmup = nb, thin = nt,
                     pars = c('lambda','sigma','mu','p','cov'),
                     open_progress = FALSE, refresh = 1000, cores = nc)
      e2 <- Sys.time()
      sem.run[ii] <- e2 - s2
      # Sys.time()
      

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # extract Stan samples
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      # out1 <- extract(m1)
      # out2 <- extract(m2)  
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # save relevant MVN results
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      results.mvn[ii,1,] <- summary(m1)$summary[1,] # correlation (rho)
      results.mvn[ii,2:3,] <- summary(m1)$summary[2:3,] # sigmas (heterogeneity)
      results.mvn[ii,4:5,] <- summary(m1)$summary[4:5,] # means
      results.mvn[ii,6,] <- summary(m1)$summary[6,] # detection (p) 
      results.mvn[ii,7,] <- summary(m1)$summary[7,] # covariance 
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # save relevant MVN results
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      results.sem[ii,1,] <- summary(m2)$summary[1,] # loading (lambda)
      results.sem[ii,2,] <- summary(m2)$summary[2,] # sigma (heterogeneity)
      results.sem[ii,3:4,] <- summary(m2)$summary[3:4,] # means
      results.sem[ii,5,] <- summary(m2)$summary[5,] # detection (p)     
      results.sem[ii,6,] <- summary(m2)$summary[6,] # lambda * sigma^2
      # print(Sys.time())
}

#save.image('E:/O_Final/Methods/latent_fitness/simulation_Stan/Stan_latent_code/cjs_SEM_gamma_100sims_5k.RData')
save.image("C:/Users/thomas.riecke/Desktop/Simulation_runs/SEM_hazard/res1k.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
o1 <- extract(m1)
o2 <- extract(m2)

plot(o1$rho, type = 'l')
plot(o2$lambda, type = 'l')
