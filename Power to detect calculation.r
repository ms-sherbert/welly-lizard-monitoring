#--- Power calculations for repeated observations on a closed population ---#

# Here, power is defined as the probability of detection of an individual or species if present at a site

# To modify this code:
# [1] Change the values in the n.reps vector to calculate power for different numbers of replicates
# [2] Change the individual detection probabilities for your system

n.reps <- c(1:11) #[1]
# For individual detection probability n.reps is the number of surveys each a site
# For species detection probability n. reps is the number of individuals at a site (i.e. lambda in Poisson N-mixture models)

p.trial <-c(0.01,0.02,0.03,0.05,0.06,0.23) #[2]
M <- length(n.reps)
T <-length(p.trial)
p.out <- matrix(NA, M, T) #Sets up a matrix of NA values to write the for loop outputs into
colnames(p.out) <- p.trial
rownames(p.out) <- n.reps

for(i in 1:T) {
power <- 1 - (1-p.trial[i])^n.reps
p.out[,i] <- power
}

p.out #To view populated matrix of power cacluations for each of your scenarios

# Export matrix of power calculatons as a .csv file to your working directory
library(MASS)
write.matrix(p.out,file="Power_calcs.csv")

# Note that you may need to change your working directory depnding on defaults
# To check location of your working directory, use getwd()
# To change location of working directory, use something like setwd("D:/R_outputs")
