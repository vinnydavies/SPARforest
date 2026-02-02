##########################################################################
#### This file illustrates the fitting of the SPAR-Forest-ERF model with a 
#### linear ERF and no measurement error on simulated data. The commonly 
#### used GLMM with linear confounder-response associations is also fitted
##########################################################################

##############################################
#### Load the libraries and functions required
##############################################
library(sf)
library(spdep)
library(tidyverse)
library(ranger)
library(INLA)

source("SPARForestERFlinear.R")



#####################
#### Read in the data
#####################
load(file="Illustrative data for SPAR-Forest-ERF.Rdata")
Y <- dat$Y
e <- dat$e
x <- dat$x
z1 <- dat$z1
z2 <- dat$z2
z3 <- dat$z3
z4 <- dat$z4
z5 <- dat$z5
z6 <- dat$z6
z7 <- dat$z7
z8 <- dat$z8
W <- dat$W
rm(dat)



##############################
#### Format the data arguments
##############################
#### Number of data points
K <- length(Y)


#### Create an INLA graph object
W.list <- mat2listw(W, style="B")
graph <- nb2INLA("DZgraph", W.list$neighbours)
graph <- inla.read.graph(filename = "DZgraph")


#### Create the indicator variable for INLA
region <- 1:K



############################################
#### Specify the parameters of the algorithm
############################################
#### Fixed parameters
ntrees <- 500 # Number of trees to use when tuning the random forest
epsilon <- 0.0005 # Threshold value for stopping the algorithm
max.iter <- 30 # Maximum number of iterations allowed for the algorithm
n.sample <- 100 # Number of samples to draw from the spatial model to propagate uncertainty into the random forest
ntrees.sample <- 10 # Number of trees for each sampled response in the random forest


#### Specify the possible tuning parameter values for the random forest
mtry <-  c(3, 5, 8)
minnode <- c(5, 10, 20)



##########################################################################
#### Specify the half-normal prior on the standard deviation scale in INLA
##########################################################################
HN.prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);  
"



####################################################################
#### Fit the GLMM model with linear confounder-response associations
####################################################################
#### Create a final data set as a data.frame
dataset <- data.frame(Y=Y, e=e, region=region, x=x, z1=z1, z2=z2, 
                      z3=z3, z4=z4, z5=z5, z6=z6, z7=z7, z8=z8)


#### Specify the model fitting formula 
formula.glmm <- Y ~ offset(log(e)) + x + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 +
      f(region, model='bym2', graph = graph, hyper=list(theta1=list(prior = HN.prior)))   


#### Fit the model
mod.glmm <- inla(formula=formula.glmm, data=dataset, family="poisson",
                     control.compute=list(dic=TRUE, mlik=TRUE, waic=TRUE, cpo=TRUE, config=TRUE),
                     control.fixed=list(mean=0, mean.intercept=0, prec=0.00001, prec.intercept=0.00001),
                     control.predictor=list(compute=TRUE, link=1), control.inla=list(strategy="auto"))   
summary(mod.glmm)


#### Extract the estimated exposure effect and its 95% credible interval
mod.glmm$summary.fixed[2, c(1,3,5)]
## The true value used to generate the data is 0.2



########################################################
#### Fit the SPAR-Forest-ERF algorithm with a linear ERF
########################################################
#### Specify the confounder data frame
Z <- data.frame(z1, z2, z3, z4, z5, z6, z7, z8)


#### Run the model
mod.SPARForest.ERF <- SPARForestERFlinear(Y=Y, e=e, x=x, Z=Z, W=W, epsilon=epsilon, 
                                          max.iter=max.iter, ntrees=ntrees, mtry=mtry, 
                                          minnode=minnode, n.sample=n.sample, 
                                          ntrees.sample=ntrees.sample)
 
  
#### Visualise the results
summary(mod.SPARForest.ERF)
mod.SPARForest.ERF$mod.spatial$summary.fixed[2, c(1,3,5)]
## The true value used to generate the data is 0.2
