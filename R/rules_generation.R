#'Rule generation procedure
#'
#'@import mboost
#'@import rpart
#'@importFrom survival Surv
#'@importFrom stats as.formula predict rexp
#'@param dat Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)
#'@param depth The mean depth of each tree-base function (Default is 2)
#'@param learnrate Shrinkage rate for each boosting step (Default is 0.01)
#'@param ntrees Number of trees (Default is 333)
#'@param sampfrac The fraction size of training sample for rule generation  (Default is NULL, means min(n / 2, 100 + 6 * sqrt(n), where n is the sample size)
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item main: main effect rules
#'    \item treat: treatment effect rules
#'  }

rule_gen <- function(dat,              # Input data frame
                     depth = 2,        # The mean depth of each tree-base function
                     learnrate = 0.01, # Shrinkage rate for each boosting step
                     ntrees = 333,     # Number of trees
                     sampfrac = NULL   # The fraction size of training sample for rule generation
) {
  
  # Extract the variable names from the dataframe 'dat'
  y_names <- colnames(dat)[1]  # The response variable is in the first column
  d_names <- colnames(dat)[2]  # The censoring/event indicator is in the second column
  z_names <- colnames(dat)[3]  # The indicator variable is in the third column
  x_names <- colnames(dat)[-c(1:3)]  # The covariates are in all other columns
  
  # Recreate a new dataset for model learning
  data.learn <- data.frame(y = dat[, y_names], z = dat[, z_names], dat[, x_names])
  colnames(data.learn)[1] <- y_names
  colnames(data.learn)[2] <- z_names
  colnames(data.learn)[-c(1:2)] <- x_names
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(dat)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Generalized Boosting Machine (GBM) for rule generation (reference the 'pre' package, Fokkema et al. 2020) ===============================================
  
  # Create a formula for the Cox Proportional-Hazards model
  formula <- as.formula(paste0("pseudo_y ~ ", paste0(c(z_names, x_names), collapse = " + ")))
  
  # Create a survival object using survival package
  y <- Surv(dat[, y_names], dat[, d_names])
  
  # Initialize estimates
  eta <- rep(0, times = nrow(dat))
  
  # Get the negative gradient function for the Cox Proportional-Hazards model
  ngradient_CoxPH <- mboost::CoxPH()@ngradient
  
  # Create a new dataset including the survival object and without the original response variable
  data_with_y_learn <- cbind(dat[, -which(names(dat) == y_names), drop = FALSE], y)
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn$pseudo_y <- ngradient_CoxPH(y = y, f = eta, w = rep(1, nrow(dat)))
  
  # Initialize an empty vector to store the rules
  rules <- c()
  
  # Iterate over ntrees to build individual trees and extract rules
  for(i in 1:ntrees) {
    
    # Set up controls for rpart tree
    tree.control <- rpart::rpart.control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- rpart::rpart(formula, control = tree.control,
                         data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose the tree into rules
    paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
    paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
    
    # Remove the root rule
    paths <- paths[-1]
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn$pseudo_y <- ngradient_CoxPH(y = y, f = eta, w = rep(1, nrow(dat)))
    
  }
  
  return(rules)
}  
