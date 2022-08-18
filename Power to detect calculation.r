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

for(i in 1:T) {
power <- 1 - (1-p.trial[i])^n.reps
p.out[,i] <- power
}

p.out #To view populated matrix of power cacluations for each of your scenarios