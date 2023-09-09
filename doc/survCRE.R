## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = TRUE, message=FALSE, warning=FALSE, results = "hide"-------------
library(speff2trial)

## Input dataset 
data(ACTG175)
data <- ACTG175
data <- data[data$cd40>= 200 & data$cd40<= 500,] # Select subjects whose baseline hazard is between 200 and 500 (Hammer et al.1996)
data <- data[data$arms == 0|data$arms == 1,] # Select the ZDL monotherapy as control,ZDL+ DID combination therapy as treatment 

## Reconstruct the dataset for proposed method
y <- data$days # Survival times
d <- data$cens # Censoring indicators
z <- data$arms # Treatment indicators
x.con <- data[,c("cd40","cd80","age","wtkg","karnof")] # Five continuous variables
x.con <- data.frame(apply(x.con,2,as.numeric))
x.bin <- data[,c("hemo","homo","drugs","race","gender","str2","symptom")] # Seven binary variables
x <- cbind(x.con,x.bin)
dat <- data.frame(Y = y, D = d, Z = z, x) # Construct the dataframe for analysis


## ---- echo = TRUE, message = FALSE, warning = FALSE---------------------------
library(survCRE)
library(ggplot2)

# Hyperparameters setting
meandepth <- 2
learnrate <- 0.01
ntrees <- 333

## Build the model
set.seed(1)
fit <- survCRE(dat = dat,meandepth = meandepth, learnrate = learnrate, ntrees = ntrees)

## ---- echo = TRUE, message = FALSE, warning = FALSE, results = "asis", out.width = "95%"----
## Base function importance
## support > 0.1 & variable importance > average variable importance
prop.vip <- prop_imp(fit,dat)
res.prop <- prop.vip[[1]]
res.prop.selected <- res.prop[res.prop[,1] > mean(res.prop[,1])& res.prop[,3] > 0.1,][,1:3]
pander::pandoc.table(res.prop.selected)

## ---- echo = TRUE, message = FALSE, warning = FALSE, fig.height = 5, fig.width = 10, out.width = "95%"----
## Variable importance
prop.vip <- sort(prop.vip[[2]],decreasing = TRUE)
variable.name <- factor(names(prop.vip), levels = names(prop.vip))
pic.prop.dat <- data.frame(var = variable.name, value = prop.vip)
gp <- ggplot(pic.prop.dat,aes(y = value,x = var)) + geom_bar(stat="identity")
gp <- gp + ylab("Variabel importance") + xlab("") + theme_minimal()
gp <- gp + theme(text=element_text(size=15)) + labs(title = "Proposed Method")
gp

## ---- echo = TRUE, message = FALSE, warning = FALSE---------------------------
## Random drawn 5 subjects from dataset and estimate their HTE
set.seed(1)
subject.id <- sample(1:nrow(dat),5)
pred <- predict_survCRE(fit,dat.test = dat[subject.id,])
hte.ran <- cbind(hte = round(pred,digit = 2),dat[subject.id,-c(1:3)])
pander::pandoc.table(hte.ran)

