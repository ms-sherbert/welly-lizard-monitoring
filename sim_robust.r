#--- Update of "Ricker DM model simulation.r" to simulate robust design ---#
#--- In development, doesn't work right now ---#

library(unmarked)

#--- List of simulation inputs (to be varied by scenario) ---#
M <- 20 #Number of sites with at least one individual lizard detected
T <- 10 #Number of primary survey periods (here, this is an annual survey for 10 years)
J <- 5 #Number of secondary survey periods
lambda <- 11 #Assumed initial average abundance per 100m2 site based on Wildland Consultants (2022).
r <- -0.02 #Annual population change rate of ca. -1% over a ten-year period (r = 0 means stable population)
K <- 100 #Carrying capacity at site i - currently set high enough to not constrain population growth rate
p <- 0.23 #Individual detection probability for primary period

#--- Simulate a single scenario under Ricker population growth model ---#

# Simulate Poisson-distributed local abundance across sites
y <- N <- matrix(NA, M, T)
N[,1] <- rpois(M, lambda)
    for(t in 2:T) {
        mu <- N[,t-1]*exp(r*(1-N[,t-1]/K))
        N[,t] <- rpois(M, mu)
    }
	
	
y[] <- rbinom(M*T, N, p) # Create simulated counts based on undelying model

#--- Prepare simulated data set for unmarked analysis ---#
umf <- unmarkedFramePCO(y = y, numPrimary=T)
summary(umf)
plot(umf) #optional if you want a visual representation of the count data for each year and site

#--- Fit model and extract parameter estimates ---#
m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=30,dynamics="ricker") # K should be higher than maximum N(it) and sufficiently high to not constrain estimates of N(it).
summary(m1) #Provides parameter estimates and the associated standard errors and p-values.

#--- Power analysis ---#
# Set desired effect sizes to pass to coefs
effect_sizes <- list(gamma)
# Run power analysis and look at summary
pa <- powerAnalysis(m1, coefs=effect_sizes, alpha=0.05)