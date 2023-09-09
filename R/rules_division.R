#'Rule division procedure
#'
#'@param dat Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)
#'@param rules The set of rules created from the generalized boosting machine
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item main: main effect rules
#'    \item treat: treatment effect rules
#'  }

rules_divide <- function(dat, rules) {
  
  # Extracting the column names for different types of variables
  y_names <- colnames(dat)[1]  # Name of the survival time column
  d_names <- colnames(dat)[2]  # Name of the censoring/status indicator column
  z_names <- colnames(dat)[3]  # Name of the treatment indicator column
  x_names <- colnames(dat)[-c(1:3)]  # Names of covariate columns
  
  
  # Identifying rules that are related to treatment
  rules.treat.ini <- rules[grep(z_names, rules)]
  
  # Remove the portion of the rule that is about the treatment indicator.
  # This is done to derive the subgroup related to HTE.
  rules.treat.gen0 <- gsub(paste0(z_names, ">=0.5 & "), "", rules.treat.ini)
  rules.treat.gen1 <- gsub(paste0(" & ", z_names, ">=0.5"), "", rules.treat.gen0)
  rules.treat0 <- gsub(paste0(z_names, "< 0.5 & "), "", rules.treat.gen1)
  rules.treat.t <- gsub(paste0(" & ", z_names, "< 0.5"), "", rules.treat0)
  
  # Remove rules that are just the treatment indicator (i.e., the root rules)
  rules.treat <- rules.treat.t[(rules.treat.t != paste0(z_names, "< 0.5") & rules.treat.t != paste0(z_names, ">=0.5"))]
  
  # Remove rules which only consist on treatment indicator
  rules.remove <- rules.treat.t[!((rules.treat.t != paste0(z_names, "< 0.5") & rules.treat.t != paste0(z_names, ">=0.5")))]
  print(paste0("remove_rules = ", length(rules.remove)))
  
  # Identifying rules that are for the main effect
  # These are the rules that are not part of the treatment (HTE) rules
  rules.main <- rules[!is.element(rules, rules.treat.ini)]
  
  # Return the two sets of rules
  return(list(treat = rules.treat, main = rules.main))
}