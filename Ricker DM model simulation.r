#---- Ricker DM data simulation ----#

library(unmarked)

# Simulate one data set under Ricker population growth model

# List of inputs
M <- 30 #Number of sites with at least one detection
T <- 10 #Number of observations (=surveys) per site (here, this is an annual survey for 10 years)
lambda <- 11 #Assumed initial average abundance per 100m2 site based on Wildland Consultants (2022).
r <- 0.90 #Annual population change rate
K <- 30 #Carrying capacity at site i - currently set high enough to not constrain population growth rate
p <- 0.20 #Individual detection probability for primary period

#Simulate a single scenario under Ricker population growth model
y <- N <- matrix(NA, M, T)
N[,1] <- rpois(M, lambda)
    for(t in 2:T) {
        mu <- N[,t-1]*exp(r*(1-N[,t-1]/K))
        N[,t] <- rpois(M, mu)
    }
y[] <- rbinom(M*T, N, p)

#--- Prepare simulated data set for unmarked analysis ---#
umf <- unmarkedFramePCO(y = y, numPrimary=T)
plot(umf)
summary(umf)

#--- Fit model and extract parameter estimates ---#

m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=30,dynamics="ricker") # K should be higher than maximum N(it) and sufficently high to not constrain estimates of N(it).
summary(m1) #Provides parameter estimates and the associated standard errors and p-values. 

#Extract back-transformed parameter estimates
lam <- exp(coef(m1, type="lambda"))
gam <- exp(coef(m1, type="gamma")) #In Ricker model this is the population change rate (r)
om <- plogis(coef(m1, type="omega")) #In Ricker model this is the carrying capacity
p <- plogis(coef(m1, type="det"))


#--- Power analysis ---#
# Set desired effect sizes to pass to coefs
effect_sizes <- list(gamma,lambda)
# Run power analysis and look at summary
pa <- powerAnalysis(m1, coefs=effect_sizes, alpha=0.05)



#----- Additional stuff -----#

# Simulate and model data under Ricker model nsim times

nsim <- 10 # To be 1000, set to 10 while in development
simout.mle <- simout.cover <- matrix(NA, nsim, 4)  # (3)
colnames(simout.mle) <- colnames(simout.cover) <- c("Lambda", "r", "K", "p") # (3)

# Simulate
set.seed(345489)
lambda <- 11 # (4)
r <- -0.01     # (4)
K <- 20       # (4) #Currently set quite high so that growth not constrained by carrying capacity
p <- 0.20    # (4)
for(i in 1:nsim) {
    cat("\ndoing", i, "at", format(Sys.time()), "\n")
    data.name <- paste("data_", stub, "_", i, sep="")
    data.file <- paste("data/", data.name, ".gzip", sep="")
    sim.i <- sim.ricker(lambda=lambda, r=r, K=K, p=p, nSites=20, nYears=10) # (1)
    cat("   max(N)=", max(sim.i$N), "\n", sep="")
    assign(data.name, sim.i)
    save(list=data.name, file=data.file) # save simulated data
    umf.i <- unmarkedFramePCO(y=sim.i$y, numPrimary=ncol(sim.i$y))
    fm.i <- pcountOpen(~1, ~1, ~1, ~1, umf.i, K=200, dynamics="ricker",  # (5)
                       starts=c(log(lambda), log(r), log(K), qlogis(p))) # (5)
    model.name <- paste("fm_", stub, "_", i, sep="")
    model.file <- paste("fm/", model.name, ".gzip", sep="")
    assign(model.name, fm.i)
    save(list=model.name, file=model.file) # save fitted models
    mle.i <- coef(fm.i)
    simout.mle[i,] <- c(exp(mle.i[1]), exp(mle.i[2]), exp(mle.i[3]), plogis(mle.i[4])) # (6)
    cat("   mle=", round(simout.mle[i,], 4), "\n")
    write.csv(simout.mle, paste(stub, "_mle.csv", sep=""),
              row.names=FALSE)
    CI.lam <- confint(fm.i, type="lambda") # (6)
    CI.r <- confint(fm.i, type="gamma")    # (6)
    CI.K <- confint(fm.i, type="omega")    # (6)
    CI.p <- confint(fm.i, type="det")      # (6)
    simout.cover[i,1] <- log(lambda) >= CI.lam[1] & log(lambda) <= CI.lam[2] # (6)
    simout.cover[i,2] <- log(r) >= CI.r[1] & log(r) <= CI.r[2]               # (6)
    simout.cover[i,3] <- log(K) >= CI.K[1] & log(K) <= CI.K[2]               # (6)
    simout.cover[i,4] <- qlogis(p) >= CI.p[1] & qlogis(p) <= CI.p[2]         # (6)
    write.csv(simout.cover, paste(stub, "_cover.csv", sep=""),
              row.names=FALSE)
    rm(list=c(data.name, model.name))
}
