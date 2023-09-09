#'Base function importance & Variable importance
#'
#'@param mod.hat a object created from function survCRE() 
#'@param dat In put dataset
#'
#'@return an object with attributes:
#'  \itemize{
#'    \item base: base function importance
#'    \item var: variable importance
#'  }
#'  
#' @export

prop_imp <- function(mod.hat,dat){
  
  # Prepare for importance
  
  beta.hat <- mod.hat$beta.hat
  
  c_id <- grep("_c",rownames(beta.hat))
  t_id <- grep("_t",rownames(beta.hat))
  
  beta.hat0 <- beta.hat[c_id,]
  base.hat0 <- unlist(strsplit(rownames(beta.hat)[c_id],"_"))
  base.hat0 <- base.hat0[base.hat0 != "c"]
  names(beta.hat0) <- base.hat0
  
  beta.hat1 <- beta.hat[t_id,]
  base.hat1 <- unlist(strsplit(rownames(beta.hat)[t_id],"_"))
  base.hat1 <- base.hat1[base.hat1 != "t"]
  names(beta.hat1) <- base.hat1
  
  # Base function importance =================================
  
  # Rule importance : 
  
  trt.id0 <- grep("(<|>=)",base.hat0)
  trt.hat0 <- base.hat0[trt.id0]
  trt.mat0 <- tran_rules(dat,trt.hat0)
  trt.sd0 <- apply(trt.mat0,2,sd)
  
  trt.id1 <- grep("(<|>=)",base.hat1)
  trt.hat1 <- base.hat1[trt.id1]
  trt.mat1 <- tran_rules(dat,trt.hat1)
  trt.sd1 <- apply(trt.mat1,2,sd)
  
  trt.sd <- trt.sd0 <- trt.sd1
  trt.mat <- trt.mat0 <- trt.mat1
  
  rules.imp <- abs(beta.hat1[trt.id1] - beta.hat0[trt.id0])*trt.sd
  hazard.ratio.trt <- exp(beta.hat1[trt.id1] - beta.hat0[trt.id0])
  support.trt <- apply(trt.mat,2,mean)
  trt.imp.sum <- cbind(Importance = rules.imp,Harzard_Ratio = hazard.ratio.trt,Support = support.trt)
  
  # Linear importance :
  
  lins.id0 <- setdiff(1:length(base.hat0),trt.id0)
  lins.hat0 <- base.hat0[lins.id0]
  lins.sd0 <- mod.hat$x_scale[lins.hat0]*0.4
  
  lins.id1 <- setdiff(1:length(base.hat1),trt.id1)
  lins.hat1 <- base.hat1[lins.id1]
  lins.sd1 <- mod.hat$x_scale[lins.hat1]*0.4
  
  lins.sd <- lins.sd0 <- lins.sd1
  
  lins.imp <- abs(beta.hat1[lins.id1] - beta.hat0[lins.id0])*lins.sd
  hazard.ratio.lins <- exp(beta.hat1[lins.id1] - beta.hat0[lins.id0])
  support.lins <- rep(1,length(lins.imp))
  lins.imp.sum <- cbind(Importance = lins.imp,Harzard_Ratio = hazard.ratio.lins,Support = support.lins)
  
  # Summarize the base importance
  
  base.imp <- rbind(trt.imp.sum,lins.imp.sum)
  base.imp <- base.imp[order(base.imp[,1],decreasing = TRUE),]
  base.imp[,1] <- base.imp[,1]*100/max(base.imp[,1])
  base.imp <- apply(base.imp,2,round,digit = 2)
  
  # Variable importance ========================
  
  var_names <- names(mod.hat$x_scale)
  trt.hat.names <- rownames(trt.imp.sum)
  lins.hat.names <- rownames(lins.imp.sum)
  
  variable_importance <- rep(0,length(var_names))
  names(variable_importance) <- var_names
  
  for(i in 1:length(var_names)){
    
    stack_trt_importance <- c()
    stack_lins_importance <- 0
    
    if(is.element(var_names[i],lins.hat1)){
      stack_lins_importance <- stack_lins_importance + lins.imp.sum[var_names[i],1]
    }
    
    for(ii in 1:nrow(trt.imp.sum)){
      
      trt.hat <- trt.imp.sum[ii,]
      trt.hat.names <- rownames(trt.imp.sum)[ii]
      
      trt.var0 <- unlist(strsplit(trt.hat.names," & "))
      trt.var1 <- unlist(strsplit(trt.var0,"(<|>=)"))
      trt.var2 <- trt.var1[1:length(trt.var1)%%2 == 1]
      trt.num <- length(trt.var2)
      stack_trt_importance[ii] <- (trt.hat[1]/trt.num)*(is.element(var_names[i],trt.var2))
      #print(is.element(var_names[i],trt.var2))
    }
    variable_importance[i] <- stack_lins_importance + sum(stack_trt_importance)
  }
  
  #var.imp <- variable_importance[order(variable_importance,decreasing = TRUE)]
  var.imp <- variable_importance
  var.imp <- var.imp*100/max(var.imp)
  var.imp <- round(var.imp,digits = 2)
  
  return(list(base = base.imp,var = var.imp))
}
