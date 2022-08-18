#--- Determinisitc Ricker function plots ---#

par(mfrow=c(2,3)) #Optional for making panel plots

K <- 200 #Carrying capacity across all sites
lambda <- 10
n.sites <- 10
Ninit <- lambda * n.sites #Initial N across all monitoring sites
T <- 31


Abundance <- rep(NA, times = T) #Sets up a vector of NA values to write the for loop outputs into
Abundance[1] <- Ninit

#--- Make Ricker function plots with varying annual change rates ---#

r <- 0.5 #For changing the annual population change rate, r = -0.02 produces a ca. 10% decrease after 10 years

for(i in 2:T) {
Ni <- Abundance[i-1] * exp(r*(1 - Abundance[i-1]/K))
Abundance[i] <- Ni
}

Year<-c(1:T)
plot(Abundance~Year,pch="",main=r,xlab="Year",ylab="Abundance")
lines(Abundance~Year)

# Add indication of 10-year monitoring period
y <-c(0,Abundance[11])
x <-c(11,11) 
lines(x,y,col="red")
