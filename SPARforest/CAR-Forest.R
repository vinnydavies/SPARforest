CARForest <- function(X.training, X.test, Y.training, coords.training, coords.test, n_tr, m_try, min_node, D, R, logscale, inla_samples)
{
####################
#### Specify lengths
####################
#### Lengths
K.training <- nrow(X.training)
K.test <- nrow(X.test)
K.total <- K.training + K.test
  
  
  
##############################
#### Initialise key quantities
##############################
#### Initialise phi
phi.training <- rep(0, K.training)
phi.test <- rep(0, K.test)
  
  
#### Neighbourhood matrix
nb.test <- knn2nb(knearneigh(coords.training, k=D))
W.nb <- make.sym.nb(nb.test)
W.training <- nb2mat(W.nb, style = "B")
W.testtraining <- array(0, dim=c(K.test, K.training))
  for(j in 1:K.test)
  {
  dist <- sqrt((coords.test[j, 1] - coords.training[ ,1])^2 + (coords.test[j, 2] - coords.training[ ,2])^2)
  indices <- order(dist)[1:D]
  W.testtraining[j, indices] <- 1
  }
W <- array(NA, dim=c(K.total, K.total))
W[1:K.training, 1:K.training] <- W.training
W[(K.training+1):K.total, 1:K.training] <- W.testtraining
W[1:K.training, (K.training+1):K.total] <- t(W.testtraining)
W[(K.training+1):K.total, (K.training+1):K.total] <- array(0, dim=c(K.test, K.test))
W.list <- mat2listw(W, style="B")
graphname <- paste("graph", ".adj", sep="")
graph <- nb2INLA(graphname, W.list$neighbours)
graph <- inla.read.graph(filename = graphname)


#### Set up the results list
results.predictions <- as.list(rep(NA, R))
names(results.predictions) <- 1:R
results.components <- as.list(rep(NA, R))
names(results.components) <- 1:R  
  

  
##################################################
#### Fit the iterative algorithm with R iterations
##################################################
   for(r in 1:R)
   {
   #### Compute the adjusted data 
   Z <- Y.training - phi.training
   dat.training <- data.frame(Z=Z, X.training)
   dat.test <- data.frame(Z=rep(0, K.test), X.test)
    
    
   #### Fit the RF model and make prediction
   mod.rf <- ranger(formula = Z~., data = dat.training, importance= 'impurity', num.trees=n_tr, mtry=m_try, min.node.size=min_node, oob.error=TRUE)
   mhat.training.oob <- mod.rf$predictions
   mhat.test <- predict(mod.rf, data=dat.test)$predictions
    
    
   #### Fit the SAR model
   dat.all <- rbind(dat.training, dat.test)
   dat.all$region <- 1:K.total
   dat.all$offset <- c(mhat.training.oob, mhat.test)
   dat.all$response <- c(Y.training, rep(NA, K.test))
   formula.car <-  response ~ offset(offset) +  f(region, model='besagproper2', graph = graph, constr=TRUE,
                                             hyper=list(lambda=list(prior="gaussian", param=c(0, 0.01)), 
                                                        prec=list(prior = "loggamma", param = c(1, 0.01))))
   mod.car <- inla(formula=formula.car, data=dat.all, family="gaussian",
                control.compute=list(dic=TRUE, mlik=TRUE, waic=TRUE, config = TRUE),
                control.fixed=list(mean=0, mean.intercept=0, prec=0.00001, prec.intercept=0.00001),
                control.family=list(hyper=list(prec=list(prior="loggamma", param=c(1, 0.01)))),
                control.predictor=list(compute=TRUE, link=1), scale=1, control.inla=list(strategy="auto"))   
   
   
   
    #### Update the estimate of phi
    phi.training <- mod.car$summary.random$region$mean[1:K.training]
    phi.test <- mod.car$summary.random$region$mean[(K.training+1):K.total]
    
    
    #### Predict from the SAR model
    samples.parameters <- inla.posterior.sample(n=inla_samples, result=mod.car, add.names=TRUE, use.improved.mean = TRUE, skew.corr=TRUE)
    samples.test <- array(NA, c(inla_samples, K.test))
      for(s in 1:inla_samples)
      {
      sigma.current <- sqrt(1 / samples.parameters[[s]]$hyperpar[1])
      fitted.current <- samples.parameters[[s]]$latent[(K.training+1):K.total]
      samples.test[s, ] <- rnorm(n=K.test, mean=fitted.current, sd=rep(sigma.current , K.test))
      }
    
    
    #### Generate the predictions
    results.test <- data.frame(prediction=rep(NA, K.test), LCI=rep(NA, K.test), UCI=rep(NA, K.test))
    results.test$prediction <- apply(samples.test, 2, mean)
    results.test[ ,2:3] <- t(apply(samples.test, 2, quantile, c(0.025, 0.975)))
    
    
    #### Backtransform them to the exponentiated scale if required
      if(logscale)
      {
      sigma2.mean <- 1 / mod.car$summary.hyperpar[1,1]
      results.test$prediction <- exp(results.test$prediction + sigma2.mean / 2)
      results.test$LCI <- exp(results.test$LCI)
      results.test$UCI <- exp(results.test$UCI)   
      }else
      {}
    
    
    #### Save the results
    components <- data.frame(mhat=c(mhat.training.oob, mhat.test), phi=c(phi.training, phi.test))
    results.predictions[[r]] <- results.test
    results.components[[r]] <- components
   }
  
  
#### Return the results
return(list(predictions=results.predictions, components=results.components))
}






