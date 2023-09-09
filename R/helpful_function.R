#'Reconstruct the data for rule ensemble procedure/prediction 
#'
#'@importFrom stats quantile sd
#'@param dat.test Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)
#'@param rules.main set of main effect rules
#'@param rules.treat set of treatment effect rules
#'@param win winsorized rate (Default is 0.025)
#'@param x_scale normalized scale (Default is NULL, unspecified)
#'@param treat.index vector of specified treatment indicator (Default is NULL, similar to the input data)
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item Y: vector of outcome
#'    \item D: vector of censoring indicator
#'    \item W: vector of treatment indicator
#'    \item X: matrix of covariates
#'    \item g_index: vector of group memebership
#'    \item x_scale: vector of normalized scale for linear terms
#'    \item treat: set of treatment effect rules
#'    \item main: set of main effect rules  
#'  }

pre_rules0 <- function(dat.test, rules.main, rules.treat, win = 0.025, x_scale = NULL,treat.index = NULL){
  
  # Extracting the column names for different types of variables
  y_names <- colnames(dat.test)[1]  # Name of the survival time column
  d_names <- colnames(dat.test)[2]  # Name of the censoring/status indicator column
  z_names <- colnames(dat.test)[3]  # Name of the treatment indicator column
  x_names <- colnames(dat.test)[-c(1:3)]  # Names of covariate columns
  
  # 1: Prepare the linear terms
  
  # Extracting the linear terms (covariates)
  linear <- dat.test[, x_names]
  
  # Apply winsorization to each column in the linear terms data frame
  # (Note: Winsorizing is the process of altering "extreme" values in the data set)
  # win = 0.025 as per Friedman and Popescu (2008)
  linear.win <- apply(linear, 2, function(x) {
    
    # Skip winsorization for binary variables
    if (length(unique(x)) > 3) {
      upper <- quantile(x, 1 - win)  # Calculate the upper quantile
      bottom <- quantile(x, win)     # Calculate the lower quantile
      
      # Replace values above the upper quantile with the upper quantile value
      x[x >= upper] <- upper
      
      # Replace values below the lower quantile with the lower quantile value
      x[x <= bottom] <- bottom
    }
    
    # Return winsorized variable
    return(x)
  })
  
  # Scaling the linear terms
  # If the scale parameter (x_scale) is provided, use it to scale the data
  if (!is.null(x_scale)) {
    linear.gen <- scale(linear.win, center = FALSE, scale = x_scale)
  } else {
    # Calculate the standard deviation for each column
    linear.sd <- apply(linear.win, 2, sd)
    
    # Replace any standard deviation of 0 with 1 to avoid division by zero
    linear.sd[linear.sd == 0] <- 1
    
    # Set the scale parameter to the standard deviation divided by 0.4
    x_scale <- linear.sd / 0.4
    
    # Scale the data
    linear.gen <- scale(linear.win, center = FALSE, scale = x_scale)
  }
  
  # 2: Combine the rule terms and linear terms. Reconstruct the data for group lasso
  
  # If the treatment index is not specified, use the treatment index in the input data
  
  if(is.null(treat.index)){
    
    # Generate rule variables for treatment and main effects
    rulevars.treat.gen <- tran_rules(dat.test, rules.treat)
    rulevars.main.gen <- tran_rules(dat.test, rules.main)
    
    # Construct new features for the control group
    rules.treat.ensemble0 <- paste0(c(rules.treat, x_names), "_c")
    rulevars.ensemble0 <- cbind(rulevars.treat.gen, linear.gen) * (1 - dat.test[, z_names])
    colnames(rulevars.ensemble0) <- rules.treat.ensemble0
    g_index0 <- 1:length(rules.treat.ensemble0)
    
    # Construct new features for the treatment group
    rules.treat.ensemble1 <- paste0(c(rules.treat, x_names), "_t")
    rulevars.ensemble1 <- cbind(rulevars.treat.gen, linear.gen) * (dat.test[, z_names])
    colnames(rulevars.ensemble1) <- rules.treat.ensemble1
    g_index1 <- 1:length(rules.treat.ensemble1)
    
    # Construct new features for the main effect
    rules.main.ensemble2 <- c(rules.main, x_names)
    rulevars.ensemble2 <- cbind(rulevars.main.gen, linear.gen)
    colnames(rulevars.ensemble2) <- rules.main.ensemble2
    g_index2 <- 1:length(rules.main.ensemble2)
    
    # Combine all new features
    X <- cbind(rulevars.ensemble0, rulevars.ensemble1, rulevars.ensemble2)
    g_index <- c(g_index0, g_index1, (g_index2 + max(g_index1)))
    
  }else{
    
    # Generate rule variables for treatment and main effects
    rulevars.treat.gen <- tran_rules(dat.test,rules.treat)
    rulevars.main.gen <- tran_rules(dat.test,rules.main)
    
    # Construct new features for the control group
    rules.treat.ensemble0 <- paste0(c(rules.treat,x_names),"_c")
    rulevars.ensemble0 <- cbind(rulevars.treat.gen,linear.gen)*(1 - treat.index)
    colnames(rulevars.ensemble0) <- rules.treat.ensemble0
    g_index0 <- 1:length(rules.treat.ensemble0)
    
    # Construct new features for the treatment group
    rules.treat.ensemble1 <- paste0(c(rules.treat,x_names),"_t")
    rulevars.ensemble1 <- cbind(rulevars.treat.gen,linear.gen)*treat.index
    colnames(rulevars.ensemble1) <- rules.treat.ensemble1
    g_index1 <- 1:length(rules.treat.ensemble1)
    
    # Construct new features for the main effect
    rules.main.ensemble2 <- c(rules.main,x_names)
    rulevars.ensemble2 <- cbind(rulevars.main.gen, linear.gen)
    colnames(rulevars.ensemble2) <- rules.main.ensemble2
    g_index2 <- 1:length(rules.main.ensemble2)
    
    # Combine all new features
    X <- cbind(rulevars.ensemble0,rulevars.ensemble1,rulevars.ensemble2)
    g_index <- c(g_index0,g_index1, (g_index2 + max(g_index1)))
    
  }  
  
  return(list(Y = dat.test[,y_names],D = dat.test[,d_names], W = dat.test[,z_names],X = X,g_index = g_index,x_scale = x_scale,treat = rules.treat, main = rules.main))
}  

#'Transform rules in the dataset into binary variables (Reference from "pre" package, Fokkema et al., 2020)
#'
#'@param dat Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)
#'@param rules set of rules
#'
#'@return binary matrix for rules 

tran_rules <- function(dat, rules) {
  
  # Constructing an Expression
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  
  # Evaluating the Expression
  rulevars <- eval(expr, dat) + 0
  
  # Setting Column Names
  colnames(rulevars) <- rules
  
  # Return the Transformed Data
  return(rulevars)
}


#'Remove complement rules from the rules candidate (Reference from "pre" package, Fokkema et al., 2020)
#'
#'@param dat Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)
#'@param rules set of rules
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item rules: set of rules
#'    \item rulevars: binary matrix for rules
#'  }

remove_complement <- function(dat, rules) {
  
  # Step 1: Transform Rules to Binary Values
  # Call the tran_rules function to convert the rules into a binary matrix
  rulevars <- tran_rules(dat, rules)
  
  # Step 2: Calculate Column Variances
  # Compute the standard deviation (used as a proxy for variance here) for each column in the binary matrix
  vars <- apply(rulevars, 2, sd)
  
  # Step 3: Group Columns with Similar Variances
  # Identify columns that have almost the same standard deviation to reduce computational effort
  vars_distinct <- lapply(unique(vars), function(x) {
    idx <- which(is_almost_eq(x, vars))  # Function to check near equality of standard deviations
    list(var = x, n = length(idx), idx = idx)
  })
  
  # Step 4: Initialize Complement Indicator
  # Create a logical vector to keep track of which columns are complements
  complements <- logical(ncol(rulevars))
  
  # Step 5: Identify Complementary Columns
  # Loop through each unique standard deviation group to find complementary columns
  for (va in vars_distinct) {
    
    # Skip if there's only one column in this variance group
    if (va$n < 2L) next
    
    idx <- va$idx
    idx <- setdiff(idx, which(complements))
    if (length(idx) < 2) next
    
    n_idx <- length(idx)
    
    # Nested loop to actually compare the columns and identify complements
    for (j in 1:(n_idx - 1)) {
      if (complements[idx[j]]) next  # Skip already identified complements
      
      this_val <- rulevars[, idx[j]]
      is_compl <- which(apply(rulevars[, idx[(j + 1):n_idx], drop = FALSE], 2, function(x) all(x != this_val))) + j
      
      # Mark identified complements
      if (length(is_compl) > 0) 
        complements[idx[is_compl]] <- TRUE
    }
  }
  
  # Step 6: Filter Out Complements
  # Update rules and rulevars to only include non-complementary rules/columns
  rules <- rules[!complements]
  rulevars <- rulevars[, !complements, drop = FALSE]
  
  # Step 7: Return Results
  # Return the updated list of rules and binary matrix
  return(list(rules = rules, rulevars = rulevars))
}

#'Check near equality (Reference from 'pre' package, Fokkema et al. 2020)
#'
#'@param x value 1 
#'@param y value 2
#'@param tolerance The tolerated difference between value 1 and value 2
#'
#'@return logical: If value 1 and value 2 are similar. 

is_almost_eq <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  
  # Check for proper input: x should be numeric and of length 1.
  stopifnot(is.numeric(x), length(x) == 1L)
  
  # Calculate the absolute value of x.
  x_abs <- abs(x)
  
  # Compute the relative or absolute difference between x and y.
  # If the absolute value of x is greater than the tolerance, use relative difference; otherwise, use absolute difference.
  xy <- if (x_abs > tolerance) {
    abs(x - y) / x_abs
  } else {
    abs(x - y)
  }
  
  # Check if the computed difference is within the tolerance level.
  return(xy <= tolerance)
}

#' Breslow estimator for estimating the baseline hazard
#'
#'@param fx vector of estimated hazard ratio
#'@param time vector of outcome
#'@param delta vector of censoring indicator
#'@param t time points of interest (default is NULL, using all observed survival time point)
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item H0: vector of estimates of baseline hazard function
#'    \item time.interest: vector of time point
#'  }

basehaz.est <- function(fx,                  # Estimated hazard ratio
                        time,                # Observed survival time
                        delta,               # Censoring indicator
                        t = NULL             # Time points of interest; default is NULL
){
  
  # Step 1: Determine the Times of Interest
  # If no specific times are given (t = NULL), use the distinct event times from 'time'
  if(is.null(t)){
    time.interest <- time[delta == 1]           # Filter only for uncensored (event) times
    time.interest <- sort(time.interest)        # Sort these times in ascending order
  } else {
    time.interest <- t                          # Use the provided times
  }
  
  # Step 2: Calculate Baseline Hazard
  # Compute the baseline hazard using Breslow's estimator for each time in 'time.interest'
  H0 <- sapply(time.interest, function(x, time, delta, fx){
    y0 <- cbind(time, delta)                      # Combine time and delta into a 2-column matrix
    h0 <- apply(y0, 1, function(y, x, time, fx) { # Apply function to each row of y0
      I(x >= y[1]) * y[2] / sum(fx[time >= y[1]]) # Compute component of the Breslow estimator
    }, x = x, time = time, fx = fx)
    
    return(sum(h0))                               # Sum up the components to get the estimate at time x
  }, time = time, delta = delta, fx = fx)
  
  # Step 3: Return Results
  # Return a list containing the baseline hazard estimates and the times of interest
  return(list(H0 = unlist(H0), time.interest = time.interest))
}


