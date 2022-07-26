#--- Update of "Ricker DM model simulation.r" to simulate robust design ---#

library(unmarked)

#--- List of simulation inputs (to be varied by scenario) ---#
M <- 20 #Number of sites with at least one individual lizard detected
T <- 10 #Number of primary survey periods (here, this is an annual survey for 10 years)
J <- 5 #Number of secondary survey periods
lam <- 11 #Assumed initial average abundance per 100m2 site based on Wildland Consultants (2022).
r <- -0.02 #Annual population change rate of ca. -1% over a ten-year period (r = 0 means stable population)
K <- 100 #Carrying capacity at site i - currently set high enough to not constrain population growth rate
p <- 0.23 #Individual detection probability for primary period

#--- Simulate a single scenario under Ricker population growth model ---#

set.seed(3)

# Simulate Poisson-distributed local abundance across sites for each primary period
N <- matrix(NA, M, T)
N[,1] <- rpois(M, lam)

    for(t in 2:T) {
        mu <- N[,t-1]*exp(r*(1-N[,t-1]/K))
        N[,t] <- rpois(M, mu)
    }

# Expand matrix to add in secondary periods
y <- n <- matrix(rep(as.numeric(t(N)), each = J), nrow = nrow(N), byrow = TRUE)

# Create simulated counts based on underlying abundance matrix and detection probabilty
y[] <- rbinom(M*J*T, n, p) 

#--- Optional if you want to see what your simulated data looks like ---#

par(mfrow=c(1,2))
TotAbund <- colSums(N) 
plot(TotAbund ~ c(1:T), xlab = "Year", ylab="Total abundance")
TotCounts <- colSums(y) 
Year <- c(1:T)
Year <- rep(Year,each=J)
plot(TotCounts ~ Year, ylab="N individuals per secondary period")

par(mfrow=(c(1,1)) #Re-sets the graphics window back to one pane

#--- Prepare simulated data set for unmarked analysis ---#

umf <- unmarkedFramePCO(y = y, numPrimary=T)
summary(umf)
#plot(umf) #if you want a visual representation of the count data for each year and site


#--- Fit model and extract parameter estimates ---#
m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=40,dynamics="trend",method="Nelder-Mead") 
# K should be higher than maximum N(it) and sufficiently high to not constrain estimates of N(it).
# Used trend dynamic due to the relatively short time frame
summary(m1) #Provides parameter estimates and the associated standard errors and p-values.

# For back-transformed estimates
lamb <- exp(coef(m1, type="lambda"))
gam <- exp(coef(m1, type="gamma"))
p <- plogis(coef(m1, type="det"))

#--- Power analysis ---#
# Set desired effect sizes to pass to coefs
effect_sizes <- list(lambda=c(intercept=0),gamma=c(intercept=-0.02),det=c(intercept=0.23)) 

# Run power analysis and look at summary
pa <- powerAnalysis(m1, coefs=effect_sizes, alpha=0.05)
pa