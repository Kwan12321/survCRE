% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rules_division.R
\name{rules_divide}
\alias{rules_divide}
\title{Rule division procedure}
\usage{
rules_divide(dat, rules)
}
\arguments{
\item{dat}{Input data frame (First column is observed survival time; Second column is censoring indicator; Third column is treatment indicator; The other columns are covariates)}

\item{rules}{The set of rules created from the generalized boosting machine}
}
\value{
an object with attributes:
\itemize{
\item main: main effect rules
\item treat: treatment effect rules
}
}
\description{
Rule division procedure
}
