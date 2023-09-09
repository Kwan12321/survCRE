#'Survival Causal Rule Ensemble
#'
#'@importFrom survival Surv
#'@importFrom stats as.formula predict rexp
#'@param dat Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)
#'@param meandepth The mean depth of each tree-base function (Default is 2)
#'@param learnrate Shrinkage rate for each boosting step (Default is 0.01)
#'@param ntrees Number of trees (Default is 333)
#'@param sampfrac The fraction size of training sample for rule generation  (Default is NULL, means min(n / 2, 100 + 6 * sqrt(n), where n is the sample size)
#'@param win winsorized rate (Default is 0.025)
#'@param nfolds number of cross validation folds
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item main: main effect rules
#'    \item treat: treatment effect rules
#'  }
#'  
#'@examples ## Fit proposed method to ACTG175 dataset
#' library(speff2trial)
#' 
#' ## Input dataset
#' data(ACTG175)
#' data <- ACTG175
#' data <- data[data$cd40>= 200 & data$cd40<= 500,]
#' data <- data[data$arms == 0|data$arms == 1,]
#' 
#' ## Reconstruct dataset for proposed method
#' y <- data$days 
#' d <- data$cens 
#' z <- data$arms 
#' x.con <- data[,c("cd40","cd80","age","wtkg","karnof")] 
#' x.con <- data.frame(apply(x.con,2,as.numeric))
#' x.bin <- data[,c("hemo","homo","drugs","race","gender","str2","symptom")] 
#' x <- cbind(x.con,x.bin)
#' dat <- data.frame(Y = y, D = d, Z = z, x) # Construct the dataframe for analysis
#' 
#' ## Fit the proposed method to dataset
#' fit <- survCRE(dat)
#' 
#' ## Predict the HTE at the 0.8 quantile of observed surival time
#' pred <- predict_survCRE(fit, dat.test = dat)
#' 
#'@export
survCRE <- function(dat,              # Input data frame
                    meandepth = 2,        # The mean depth of each tree-base function
                    learnrate = 0.01, # Shrinkage rate for each boosting step
                    ntrees = 333,     # Number of trees
                    sampfrac = NULL,  # The fraction size of training sample for rule generation
                    win = 0.025,      # winsorized rate (Default is 0.025)
                    nfolds = 5
                    ){
                    
                    # Rules generation -----------------------------------------
                    rules.gen <- rule_gen(dat = dat,depth = meandepth, learnrate = learnrate, ntrees = ntrees, sampfrac = sampfrac)
                    print(paste0("Rule generation Complete"))  # Indicate completion of rule generation
                    print(paste0("all = ", length(rules.gen))) # Display count of generated rules
                    
                    # Rules division -----------------------------------------
                    rules.can <- rules_divide(dat = dat, rules = rules.gen)
                    rules.treat <- rules.can$treat  # Rules related to treatment
                    rules.main <- rules.can$main  # Rules related to main effect
                    print(paste0("Rule division Complete"))             # Indicate completion of rule division
                    print(paste0("treat.ini = ", length(rules.treat)))  # Display initial count of treatment-related rules
                    print(paste0("main.ini = ", length(rules.main)))    # Display initial count of main-effect-related rules
                    
                    # Remove the dumplicates and complements to reduce the computation ----------------------------------

                    # Remove duplicates
                    rules.treat <- unique(rules.treat)
                    rules.main <- unique(rules.main)
                    
                    # Remove complements
                    rules.treat <- remove_complement(dat, rules.treat)[[1]]
                    rules.main <- remove_complement(dat, rules.main)[[1]]
                    print(paste0("Remove dumplicates & complements Complete"))
                    print(paste0("treat.fin = ", length(rules.treat)))  # Display final count of treatment-related rules
                    print(paste0("main.fin = ", length(rules.main)))  # Display final count of main-effect-related rules
                    
                    # Reconstruct the dataset for rule ensemble procedure (winsorized, normalized, grouped) --------------
                    pre_dat <- pre_rules0(dat,rules.main = rules.main,rules.treat = rules.treat,win = win)
                    
                    # Rules ensemble -------------------------------------------------------------------------------------
                    rule_ensemble <- rule_cut(pre_dat,rules.main,rules.treat, nfolds = nfolds)
                    
                    # Summarize and store the results into a list
                    res <- list(
                      mod = rule_ensemble$mod,                     # Fitted model
                      lambda = rule_ensemble$lambda,               # Optimal lambda value
                      x_scale = rule_ensemble$x_scale,             # Scale of features
                      rules.treat = rule_ensemble$rules.treat,     # Treatment rules
                      rules.main = rule_ensemble$rules.main,       # Main effect rules
                      beta.hat = rule_ensemble$beta.hat,           # Coefficients
                      basehaz = rule_ensemble$basehaz              # Baseline hazard estimates
                    )
                    print(paste0("Rules ensemble Complete"))
              
                    return(res)
}