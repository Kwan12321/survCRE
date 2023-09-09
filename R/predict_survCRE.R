#'Predict function for proposed method
#'
#'@param rule_ensemble An object from proposed method
#'@param dat.test The data used to estimate HTE
#'@param time Specify the time point to define the HTE, default is -1, the 0.8 quantile of observed survival time is used  
#'@param type Type of HTE (survival rate or RMST)
#'
#'@return a vector of HTE
#'
#'@export

# Prediction function for proposed method # ==============

predict_survCRE <- function(rule_ensemble,dat.test,time = -1,type = "survival rate"){
  
  pre_dat.test0 <- pre_rules0(dat.test,rule_ensemble$rules.main,rule_ensemble$rules.treat,treat.index = 0)
  pre_dat.test1 <- pre_rules0(dat.test,rule_ensemble$rules.main,rule_ensemble$rules.treat,treat.index = 1)
  hr0 <- predict(rule_ensemble$mod,pre_dat.test0$X,lambda = rule_ensemble$lambda,type = "response") 
  hr1 <- predict(rule_ensemble$mod,pre_dat.test1$X,lambda = rule_ensemble$lambda,type = "response") 
  
  if(type == "survival rate"){
    
    if(time < 0){
      time <- quantile(rule_ensemble$basehaz[,2],0.8)
    }
    
    # Estimated survival rate
    est0 <- exp(-hr0%*%t(rule_ensemble$basehaz[,1]))
    est1 <- exp(-hr1%*%t(rule_ensemble$basehaz[,1]))
    est <- est1[,which.min(abs(rule_ensemble$basehaz[,2] - time))] - est0[,which.min(abs(rule_ensemble$basehaz[,2] - time))]
    
  }
  
  if(type == "RMST"){
    
    if(time < 0){
      time <- quantile(rule_ensemble$basehaz[,2],0.8)
    }
    
    # Estimated RMST
    time.interest <- rule_ensemble$basehaz[,2]
    t.rmst0 <- diff(c(0,time.interest[1:which.min(abs(time.interest - time))]))
    prop.est0.t <- pre_dat.test0[,1:which.min(abs(time.interest - time))]
    est.rmst0 <- prop.est0.t%*%t.rmst0
    t.rmst1 <- diff(c(0,time.interest[1:which.min(abs(time.interest - time))]))
    prop.est1.t <- pre_dat.test1[,1:which.min(abs(time.interest - time))]
    est.rmst1 <- prop.est1.t%*%t.rmst1
    est <- est.rmst1 - est.rmst0
    
  }
  return(est)
}	