# Repository description
A repository for analysis and power analysis of lizard survey data collected using a repeated count 
protocol from unmarked individuals.

**Note**: This repository is still in development as at 29/09/2022. Some extensions can be made to the power analysis script to simulate more scenarios. 

## Repository contents:

`Power to detect calculation.r` R script for calculating the power of different numbers of secondary survey periods to detect individuals. 
Uses equation from Royle and Dorazio (2008).  

`Ricker DM model simulation.r` R script for simulating repeated count data from a Poisson-distributed population with imperfect individual detection probability over a number of primary periods. 
Each primary period consists of a single survey. Simulated population change between primary periods follows a Ricker (1954) function.
Individual detection probability is modeled as a binomial process with *p* = individual detection probability = [0,1].
Uses Dail-Madsen (2011) models to model data and as basis for power analysis in unmarked package.  Model uses `method = "Nelder-Mead"` because original model where `dynamics = "ricker"` returned a singular Hessian matrix.
`dynamics = "trend"` used because of the relatively short time frame (10 years) of the simulated data and the primary interest in detection declining trends. 
Modeled on code provided by Hostetler and Chandler (2015).

`sim_robust.r` R script that extends `Ricker DM model simulation.r` to simulate and run a power analysis on data from a robust design programme with 
T primary periods and J secondary periods. 

`Ricker function examples.r` R script for plotting Ricker functions with different population growth rates. Not fully automated, 
i.e. you'll need to run each iteration yourself to make each plot in the panel with different parameters.
Useful for selecting gamma values for simulations because 
the plots help you to see how gamma relates to net change in abundance over a defined time period. 



## References: 

Dail, D. and Madsen, L. (2011) Models for estimating abundance from repeated counts of an open metapopulation. *Biometrics 67(2)*: 577-587.

Fiske, I. and Chandler, R. (2011) unmarked: an R package for fitting hierarchical models of wildlife occurrence and abundance. *Journal of Statistical Software 43*: 1-23. 

Hostetler, J. A. and Chandler, R. B. (2015) Improved state-space models for inference about spatial and temporal variation in abundance from count data. *Ecology 96*: 1713-1723.

Ricker, W. E. (1954) Stock and recruitment. *Journal of the Fisheries Board of Canada 11*: 559-623.   

Royle, J. A. and Dorazio, R. M. (2008) *Hierarchical Modeling and Inference in Ecology: The Analysis of Data from Populations, Metapopulations and Communities.* Academic Press, London, UK. 