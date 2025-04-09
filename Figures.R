# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(vioplot)
library(latex2exp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index 1 [number of simulations]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index 2 [parameter]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   #    MVN      SEM      truth     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   1    rho      lambda   mu1       1
#   2    sigma1   sigma    mu2       2
#   3    sigma2   mu1      sigma1    3
#   4    mu1      mu2      sigma2    4
#   5    mu2      p        rho       5
#   6    p        'cov'    p         6
#   7    cov               cov       7
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index 3 [posterior summary stat order in 3rd dimension of results storage]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) mean
# 2) se_mean
# 3) sd
# 4) lower 95%
# 5) lower 50%
# 6) median
# 7) upper 50%
# 8) upper 95%
# 9) effective sample size
# 10) Rhat
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/Stan_simulation/1k_sim/res1k.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assess convergence
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conv <- rep(1, n.sims)
for (ii in 1:n.sims){
  if(any(results.sem[ii,,10] > 1.015) | 
     any(results.mvn[ii,,10] > 1.015)){
    conv[ii] <- 0
  }
}
table(conv)
used <- which(conv == 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assess convergence
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qrho <- quantile(results.mvn[used,1,1], c(0.025,0.5,0.975), na.rm = T)
qsphi <- quantile(results.mvn[used,2,1], c(0.025,0.5,0.975), na.rm = T)
qsgamma <- quantile(results.mvn[used,3,1], c(0.025,0.5,0.975), na.rm = T)
qlam <- quantile(results.sem[used,1,1], c(0.025,0.5,0.975), na.rm = T)
qseta <- quantile(results.sem[used,2,1], c(0.025,0.5,0.975), na.rm = T)
qcsem <- quantile(results.sem[used,6,1], c(0.025,0.5,0.975), na.rm = T)
qcmvn <- quantile(results.mvn[used,7,1], c(0.025, 0.5, 0.975), na.rm = T)
  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# raincloud plots figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/manuscript/Ecology_revision2/Figure4.pdf",
    width = 10, height = 4)
par(mfrow = c(1,3), family = 'serif', mar = c(2.1,2.1,2.1,2.1),
    oma = c(3,3,0,0))
# upper left
vioplot(results.mvn[used,2,1], results.mvn[used,3,1],results.sem[used,2,1], 
        wex = 0.5, side = 'right',
        col = c(adjustcolor('red',alpha = 0.5),
                adjustcolor('red',alpha = 0.5),
                adjustcolor('blue',alpha = 0.25)),
        border = c(adjustcolor('red',alpha = 0.5),
                   adjustcolor('red',alpha = 0.5),
                   adjustcolor('blue',alpha = 0.25)),
        drawRect = F,
        las = 1, names = c('','',''),
        ylim = c(0,0.8))
axis(side = 1, at = 1:3, cex.axis = 1.5, padj = 1,
     labels = c(TeX("$\\sigma_{phi}$"), TeX("$\\sigma_{\\gamma}$"), TeX("$\\sigma_{\\eta}$")))
boxplot(results.mvn[used,2,1], results.mvn[used,3,1],results.sem[used,2,1],
        boxwex = 0.1, add = T, yaxt = 'n', outline = F,
        col = c(adjustcolor('red',alpha = 0.5),
                adjustcolor('red',alpha = 0.5),
                adjustcolor('blue',alpha = 0.25)))     
points(results.mvn[used,2,1] ~ jitter(rep(0.9, length(used)), amount = 0.025), cex = 0.125, col = 'maroon')
points(results.mvn[used,3,1] ~ jitter(rep(1.9, length(used)), amount = 0.025), cex = 0.125, col = 'maroon')
points(results.sem[used,2,1] ~ jitter(rep(2.9, length(used)), amount = 0.025), cex = 0.125, col = 'navy')
mtext(TeX("Parameter estimates"), side = 2, line = 1, cex = 1.5, outer = T)

# center
vioplot(results.mvn[used,1,1], results.sem[used,1,1],
        wex = 0.5, side = 'right',
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25)), 
        border = c(adjustcolor('red',alpha = 0.5), 
                   adjustcolor('blue',alpha = 0.25)),
        drawRect = F, ylim = c(-3.5,1.75),
        las = 1, names = c('',''))
boxplot(results.mvn[used,1,1], results.sem[used,1,1],
        boxwex = 0.1, add = T, yaxt = 'n', outline = F,
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25))) 
axis(side = 1, at = 1:2, cex.axis = 1.5, padj = 1,
     labels = c(TeX("$\\rho$"),TeX("$\\lambda$")))
points(results.mvn[used,1,1] ~ jitter(rep(0.9, length(used)), amount = 0.025), cex = 0.125, col = 'maroon')
points(results.sem[used,1,1] ~ jitter(rep(1.9, length(used)), amount = 0.025), cex = 0.125, col = 'navy')



# right
vioplot(results.mvn[used,7,1], results.sem[used,6,1],
        wex = 0.5, side = 'right',
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25)), 
        border = c(adjustcolor('red',alpha = 0.5), 
                   adjustcolor('blue',alpha = 0.25)),
        drawRect = F,
        las = 1,
        ylim = c(-0.3,0.1), names = c('',''))
boxplot(results.mvn[used,7,1], results.sem[used,6,1],
        boxwex = 0.1, add = T, yaxt = 'n', outline = F,
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25))) 
points(results.mvn[used,7,1] ~ jitter(rep(0.9, length(used)), amount = 0.025), 
       cex = 0.125, col = 'maroon')
points(results.sem[used,6,1] ~ jitter(rep(1.9, length(used)), amount = 0.025), 
       cex = 0.125, col = 'navy')
axis(side = 1, at = 1:2, cex.axis = 1.5, padj = 1,
     labels = c(TeX("$\\sigma_{phi}sigma_{gamma}rho$"),
                TeX("$\\lambda\\sigma_{eta}^2$")))       
mtext('Parameter', side = 1, line = 2, cex = 1.5, outer = T)

dev.off()






















  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# figure (basic)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf("E:/O_Final/PIFL_WTDI_senescence_heterogeneity/Stan_simulation/1k_sim/Fsim.pdf",
    width = 12, height = 5)
par(mfrow = c(1,3), family = 'sans', mar = c(2.1,2.1,2.1,2.1),
    oma = c(3,3,0,0))
# upper left
vioplot(results.mvn[used,2,1], results.mvn[used,3,1],results.sem[used,2,1], 
        wex = 0.5, side = 'right',
        col = c(adjustcolor('red',alpha = 0.5),
                adjustcolor('red',alpha = 0.5),
                adjustcolor('blue',alpha = 0.25)),
        border = c(adjustcolor('red',alpha = 0.5),
                   adjustcolor('red',alpha = 0.5),
                   adjustcolor('blue',alpha = 0.25)),
        drawRect = F,
        las = 1, names = c('','',''),
        ylim = c(0,0.8))
axis(side = 1, at = 1:3, cex.axis = 1.5, padj = 1,
     labels = c(TeX("$\\sigma_{phi}$"), TeX("$\\sigma_{\\gamma}$"), TeX("$\\sigma_{\\eta}$")))
boxplot(results.mvn[used,2,1], results.mvn[used,3,1],results.sem[used,2,1],
        boxwex = 0.1, add = T, yaxt = 'n', outline = F,
        col = c(adjustcolor('red',alpha = 0.5),
                adjustcolor('red',alpha = 0.5),
                adjustcolor('blue',alpha = 0.25)))     
points(results.mvn[used,2,1] ~ jitter(rep(0.9, length(used))), cex = 0.125, col = 'maroon')
points(results.mvn[used,3,1] ~ jitter(rep(1.9, length(used))), cex = 0.125, col = 'maroon')
points(results.sem[used,2,1] ~ jitter(rep(2.9, length(used))), cex = 0.125, col = 'navy')
mtext(TeX("Estimates"), side = 2, line = 1, cex = 1.5, outer = T)

# center
vioplot(results.mvn[used,1,1], results.sem[used,1,1],
        wex = 0.5, side = 'right',
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25)), 
        border = c(adjustcolor('red',alpha = 0.5), 
                   adjustcolor('blue',alpha = 0.25)),
        drawRect = F,
        las = 1, names = c('',''),
        ylim = c(-3.5,1.5))
boxplot(results.mvn[used,1,1], results.sem[used,1,1],
        boxwex = 0.1, add = T, yaxt = 'n', outline = F,
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25))) 
axis(side = 1, at = 1:2, cex.axis = 1.5, padj = 1,
     labels = c(TeX("$\\rho$"),TeX("$\\lambda$")))
points(results.mvn[used,1,1] ~ jitter(rep(0.9, length(used))), cex = 0.125, col = 'maroon')
points(results.sem[used,1,1] ~ jitter(rep(1.9, length(used))), cex = 0.125, col = 'navy')



# right
vioplot(results.mvn[used,7,1], results.sem[used,6,1],
        wex = 0.5, side = 'right',
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25)), 
        border = c(adjustcolor('red',alpha = 0.5), 
                   adjustcolor('blue',alpha = 0.25)),
        drawRect = F,
        las = 1,
        ylim = c(-0.35,0.05), names = c('',''))
boxplot(results.mvn[used,7,1], results.sem[used,6,1],
        boxwex = 0.1, add = T, yaxt = 'n', outline = F,
        col = c(adjustcolor('red',alpha = 0.5), 
                adjustcolor('blue',alpha = 0.25))) 
points(results.mvn[used,7,1] ~ jitter(rep(0.9, length(used))), 
       cex = 0.125, col = 'maroon')
points(results.sem[used,6,1] ~ jitter(rep(1.9, length(used))), 
       cex = 0.125, col = 'navy')
axis(side = 1, at = 1:2, cex.axis = 1.5, padj = 1,
     labels = c(TeX("$\\sigma_{phi}sigma_{gamma}rho$"),
               TeX("$\\lambda\\sigma_{eta}^2$")))       
mtext('Parameter', side = 1, line = 2, cex = 1.5, outer = T)

dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# constancy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplot(results.mvn[,1,1] - truth[,5])
mean(results.mvn[,1,1] - truth[,5])
mean(results.mvn[,7,1] - truth[,7])
mean(results.sem[,6,1] - truth[,7])

length(which(results.sem[,1,8] < 0))
length(which(results.mvn[,1,8] < 0))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# convergence
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




