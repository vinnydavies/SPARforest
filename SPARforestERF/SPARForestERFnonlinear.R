######################################################
#### This file runs the SPARForestERFlinear() function
######################################################
SPARForestERFnonlinear <- function(Y, e, x, Z, W, epsilon, max.iter, ntrees, mtry, minnode, n.sample, ntrees.sample)
{
#########################
#### Format the arguments
#########################
#### Variables
logSMR <- log(Y / e)   


#### Create the grouped exposure variable
n.group <- 50
x.group <- inla.group(x, n=n.group, method="quantile")

   
#### Number of observations
K <- nrow(W)
region <- 1:K   


#### Create an INLA graph object
W.list <- mat2listw(W, style="B")
graph <- nb2INLA("DZgraph", W.list$neighbours)
graph <- inla.read.graph(filename = "DZgraph")



################################################
#### Set up quantities required by the algorithm
################################################
#### Set up the matrix to monitor convergence
converge.results <- data.frame(iter=0, difference=100)
diff <- converge.results$difference[nrow(converge.results)]
iter <- converge.results$iter[nrow(converge.results)]


#### Set up quantities need by the algorithm
ERF.current <- rep(0, K)
ERFandRE.current <- rep(0, K)
confounder.current <- rep(0,K)


#### Specify the data matrices for the different parts of the model
dat.rf <- data.frame(response=logSMR, Z)
dat.spatial <- data.frame(Y, e, x, x.group=x.group, region, confounder.current=rep(NA, K))


#### Specify the half-normal prior on the standard deviation scale in INLA
HN.prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);  
"



##################
#### Run the model
##################
   while(diff > epsilon & iter < max.iter)
   {
   ####################################
   #### Estimate the random forest part   
   ####################################
   #### Specify the formula
   dat.rf$response <- logSMR - ERFandRE.current
   formula.rf <- response ~ .
 
       
   #### Tune the random forest
   tune.combinations <- expand.grid(mtry, minnode)
   n.tune <- nrow(tune.combinations)
   tune.results <- data.frame(mtry=tune.combinations[ ,1], minnode=tune.combinations[ ,2], RMSE=rep(NA, n.tune))
      for(j in 1:n.tune)
      {
      mod.rftune <- ranger(formula=formula.rf, data=dat.rf, num.trees=ntrees, 
                       mtry=tune.results$mtry[j], min.node.size=tune.results$minnode[j])
      tune.results$RMSE[j] <- sqrt(mean((dat.rf$response - mod.rftune$predictions)^2))
      }
   row.optimal <- which(tune.results$RMSE==min(tune.results$RMSE))
   mtry.optimal <- tune.results$mtry[row.optimal]
   minnode.optimal <- tune.results$minnode[row.optimal]
 
       
   #### Fit the optimised RF model
      if(iter==0)
      {
      mod.rffinal <- ranger(formula=formula.rf, data=dat.rf, num.trees=ntrees, 
                             mtry=mtry.optimal, min.node.size=minnode.optimal)
      confounder.current <- mod.rffinal$predictions
      sigma2.oob <- mod.rffinal$prediction.error
      logprec.oob <- log(1 / sigma2.oob)
      }else
      {
      pred.samples <- array(NA, c(n.sample, K))
      response.samples <- array(NA, c(n.sample, K))
      n.oob.samples <- array(NA, c(n.sample, K))          
         for(j in 1:n.sample)
         {
         response.samples[j, ] <- logSMR - ERFandRE.samples[j, ]
         dat.rf$response <- response.samples[j, ]
         formula.rf <- response ~ .
         mod.rffinal <- ranger(formula=formula.rf, data=dat.rf, num.trees=ntrees.sample, 
                           mtry=mtry.optimal, min.node.size=minnode.optimal, keep.inbag=TRUE)
         pred.samples[j, ] <- mod.rffinal$predictions
         temp.mat <- do.call(rbind, mod.rffinal$inbag.counts)
         n.oob.samples[j, ] <- apply(temp.mat==0, 2, sum)    
         }
      confounder.current <- apply(pred.samples * n.oob.samples, 2, sum, na.rm=TRUE)   / apply(n.oob.samples, 2, sum)
      response.mean <- apply(response.samples * n.oob.samples, 2, sum, na.rm=TRUE)   / apply(n.oob.samples, 2, sum)
      sigma2.oob <- mean((response.mean - confounder.current)^2)
      logprec.oob <- log(1 / sigma2.oob)
      }

     
   ####################################  
   #### Estimate the spatial model part
   ####################################
   #### Add in the current confounder values
   dat.spatial$confounder.current <- confounder.current
   
   
   #### Specify the model fitting formula
   formula.spatial <- Y ~ offset(log(e)) +
         f(x.group, model="rw2", hyper=list(theta=list(prior=HN.prior))) +
         f(region, model='bym2', graph = graph, hyper=list(theta1=list(prior = HN.prior))) + 
         f(confounder.current, model="meb", values=confounder.current, scale=rep(1,K),
            hyper=list(theta1=list(initial=1, fixed=TRUE, prior = "gaussian", param = c(0, 0.00001)), 
                       theta2=list(initial = logprec.oob, fixed = TRUE, prior = "loggamma", param = c(1, 1e-5))))

   
   #### Fit the spatial model   
   mod.spatial <- inla(formula=formula.spatial, data=dat.spatial, family="poisson",
                         control.compute=list(dic=TRUE, mlik=TRUE, waic=TRUE, cpo=TRUE, config=TRUE),
                         control.fixed=list(mean=0, mean.intercept=0, prec=0.00001, prec.intercept=0.00001),
                         control.predictor=list(compute=TRUE, link=1), control.inla=list(strategy="auto")) 
     
   
   #### Extract the estimated ERF
   dat.gx <- mod.spatial$summary.random$x.group[ ,1:2]
   ERF.new <- rep(NA, K)
      for(j in 1:K)
      {
      ERF.new[j] <- dat.gx$mean[dat.gx$ID==dat.spatial$x.group[j]]  
      }

   
   #### Construct a matrix of samples for the exposure and ERF parts
   samples <- inla.posterior.sample(n=n.sample, result=mod.spatial)
   ERFandRE.samples <- array(NA, c(n.sample, K))
      for(j in 1:n.sample)
      {
      m <- length(samples[[j]]$latent)
      sample.lp <- samples[[j]]$latent[1:K]
      sample.confounder <- samples[[j]]$latent[(m-K):(m-1)]
      sample.intercept <- samples[[j]]$latent[m]
      ERFandRE.samples[j, ] <-  sample.lp - log(dat.spatial$e) - sample.confounder - sample.intercept
      }
   ERFandRE.current <- apply(ERFandRE.samples, 2, mean)      

      
   
   #####################################      
   #### Compute the convergence criteria
   #####################################
   iter <- iter + 1
   diff <- mean(abs((ERF.current - ERF.new)))
   converge.results <- rbind(converge.results, c(iter, diff))
   ERF.current <- ERF.new
   }
   

   
#######################
#### Return the results
#######################
#### Format the convergence summary
convergence.summary <- converge.results[-1, ]


#### Format the ERF summary
samples <- inla.posterior.sample(n=1000, result=mod.spatial)
mod.spatial$misc$configs$contents
ERF.samples <- array(NA, c(1000, n.group))      
   for(j in 1:1000)
   {
   sample.ERF <- samples[[j]]$latent[(K+1):(K+50)]
   ERF.samples[j, ] <-sample.ERF - sample.ERF[1]
   }

ERF.estimate <- data.frame(exposure=mod.spatial$summary.random$x.group$ID, ERF=apply(ERF.samples, 2, mean), 
                             LCI=apply(ERF.samples, 2, quantile, 0.025), 
                              UCI=apply(ERF.samples, 2, quantile, 0.975))


#### Format the confounder summary
confounder.estimate <- dat.spatial$confounder.current


#### Construct the final results object
results <- list(convergence.summary=convergence.summary, ERF.estimate=ERF.estimate, 
                confounder.estimate=confounder.estimate, mod.spatial=mod.spatial)   
return(results)   
}

