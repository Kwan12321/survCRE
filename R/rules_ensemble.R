#'Rule Ensemble procedure
#'
#'@importFrom grpreg cv.grpsurv
#'@importFrom stats coef
#'@param pre_dat a object created from sub function "pre_dat0"
#'@param rules.main set of main effect rules 
#'@param rules.treat set of treatment effect rules
#'@param nfolds number of cross validation folds
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item mod: Fitted model
#'    \item lambda: Optimal lambda value
#'    \item x_scale: Scale of features
#'    \item rules.treat: Treatment rules
#'    \item rules.main: Main effect rules
#'    \item beta.hat: Coefficients
#'    \item basehaz: Baseline hazard estimates
#'  }

rule_cut <- function(pre_dat, rules.main, rules.treat, nfolds = 5) {
  
  # Extract Parameters from pre_dat list
  y <- pre_dat$Y        # Survival time
  d <- pre_dat$D        # Treatment indicator
  z <- pre_dat$W        # Control covariates
  x <- pre_dat$X        # Feature variables
  x_scale <- pre_dat$x_scale  # Scale of feature variables
  g_index <- pre_dat$g_index  # Group index for features
  p <- length(x_scale)  # Number of scaled features
  
  # Perform Group Lasso
  
  fit.lasso <- glasso(x, Surv(y, d), group = g_index, nfolds = 5)
  
  # Extract best lambda and coefficients from the fitted model
  lambda <- fit.lasso$best.lambda
  beta.hat <- fit.lasso$beta
  beta.hat <- beta.hat[beta.hat != 0, , drop = FALSE]  # Keep only non-zero coefficients
  
  # Compute Baseline Hazard using Breslow estimator
  # Identify unique event times for individuals who experienced the event (d==1)
  time.interest <- y[d == 1]
  time.interest <- sort(unique(time.interest))
  
  # Predict the baseline hazard
  pre_basehaz <- predict(fit.lasso$model, as.matrix(x), lambda = lambda, type = "response")
  
  # Compute baseline hazard estimates
  basehaz <- basehaz.est(fx = pre_basehaz, time = y, d, t = time.interest)$H0
  
  # Create a data frame to store the baseline hazard and time
  basehaz <- data.frame(basehaze = basehaz, time = time.interest)
  
  # Summarize and store the results into a list
  res <- list(
    mod = fit.lasso$model,         # Fitted model
    lambda = lambda,               # Optimal lambda value
    x_scale = x_scale,             # Scale of features
    rules.treat = rules.treat,     # Treatment rules
    rules.main = rules.main,       # Main effect rules
    beta.hat = beta.hat,           # Coefficients
    basehaz = basehaz              # Baseline hazard estimates
  )
  
  # Return the result
  return(res)
}

#'Perform (adaptive) group lasso
#'
#'@param X matrix of covariates
#'@param y vector of outcome 
#'@param group vector of group membership
#'@param nfolds number of cross validation folds
#'@param adaptive logical: should adaptive group lasso be performed
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item beta: coefficients
#'    \item model: fitted model
#'    \item time: computation time
#'    \item best.lambda: Optimized tuning parameter lambda
#'    \item adpen: adaptive weight
#'  }

glasso <- function(
    X, y, group,
    nfolds = 5,
    adaptive = FALSE
){
  
  fit <- cv.grpsurv(X = X,y=y,group = group,nfolds = nfolds)
  bhat <- as.matrix(coef(fit))
  bhat.res <- bhat
  model.name <- "glasso"
  
  best.lambda <- fit$lambda.min
  
  adpen <- NA
  
  if(adaptive){
    
    bhat[bhat == 0] <- .Machine$double.eps
    adpen <- sapply(unique(group),function(x,beta) 1 / sqrt(sum(beta[x,]^2)),bhat)
    
    fit <- grpreg::cv.grpsurv(
      X = X, y = y,group = group,
      group.multiplier = adpen
    )
    
    bhat.res <- as.matrix(coef(fit))
    
    model.name <- "aglasso"
    
    best.lambda <- fit$lambda.min
  }
  
  aglasso.model <- list(
    "beta" = bhat.res,
    "model" = fit,
    "best.lambda" = best.lambda,
    "adpen" = adpen
  )
  
  return(aglasso.model)
}




