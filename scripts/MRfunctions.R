library(data.table)
library(tidyverse)
library(MendelianRandomization) 

ToMR <- function(mydata){
  if (nrow(mydata) == 1) {
    mydata2 = dat_to_MRInput(mydata)
  } else
    mydata2 = dat_to_MRInput(mydata,get_correlation=TRUE) 
} #TOMR


df_data <- function(mydata2){
  dfdata=data.table(
    cbind(mydata2[[1]]@snps,
          mydata2[[1]]@betaX, 
          mydata2[[1]]@betaY,
          mydata2[[1]]@betaXse, 
          mydata2[[1]]@betaYse)
  ) 
  
  colnames(dfdata)=c("SNP","beta.x","beta.y","se.x","se.y")
  
  dfdata$beta.x=as.numeric(dfdata$beta.x)
  dfdata$beta.y=as.numeric(dfdata$beta.y)
  dfdata$se.x=as.numeric(dfdata$se.x)
  dfdata$se.y=as.numeric(dfdata$se.y)
  
  return(dfdata)
} #×ª»¯MR·ÖÎöÎÄ¼þ


# Blank results
error_df <- function(df_mr){
  n_ins <- nrow(df_mr)
  data.table(N_ins = n_ins,
             beta = NaN, SE = NaN,
             LCI = NaN, UCI = NaN,
             P_value = NaN)
}

blank_df <- function(){
  data.table(N_ins = 0,
             beta = NA, SE = NA,
             LCI = NA, UCI = NA,
             P_value = NA, Method = NA)
}


# Single instrument MR (ratio method)
single.ins.MR <- function (df_mr){
  if (nrow(df_mr) != 1) stop ("N instrument is not 1")
  df_mr[, `:=`(beta.mr = beta.y / beta.x,
               se.mr = sqrt((se.y^2/beta.x^2) + (beta.y^2 *se.x^2/beta.x^4))
  )][,`:=`(lci.mr = beta.mr - qnorm(0.975)*se.mr,
           uci.mr = beta.mr + qnorm(0.975)*se.mr,
           P.mr = 2*pnorm(-abs(beta.mr / se.mr))
  )]
  res2 <- data.table(N_ins = 1,
                     beta = df_mr[,beta.mr],
                     SE = df_mr[,se.mr],
                     LCI = df_mr[,lci.mr],
                     UCI = df_mr[,uci.mr],
                     P_value = df_mr[,P.mr])
  return(res2)
}


# MR IVW 
MR_IVW <- function (df_mr, ldrho=matrix(),  ...) {
  res2 <- tryCatch({
    MR.input <- mr_input(bx = df_mr[,beta.x],
                         bxse = df_mr[,se.x],
                         by = df_mr[,beta.y],
                         byse = df_mr[,se.y],
                         corr = ldrho,
                         snps = df_mr[,SNP])
    MR.output <- MendelianRandomization::mr_ivw(MR.input, correl = TRUE, ...)
    
    data.table(N_ins = MR.output@SNPs,
               beta = MR.output@Estimate,
               SE = MR.output@StdError,
               LCI = MR.output@CILower,
               UCI = MR.output@CIUpper,
               P_value = MR.output@Pvalue)
  },
  error = function(e) error_df(df_mr)
  )
  return(res2)
}


# MR Egger
MR_Egger <- function (df_mr, ldrho=matrix(), ...) {
  res2 <- tryCatch({
    MR.input <- mr_input(bx = df_mr[,beta.x],
                         bxse = df_mr[,se.x],
                         by = df_mr[,beta.y],
                         byse = df_mr[,se.y],
                         corr = ldrho,
                         snps = df_mr[,SNP])
    
    MR.output <- MendelianRandomization::mr_egger(MR.input, correl = TRUE, ...)
    
    data.table(N_ins = MR.output@SNPs,
               beta = c(MR.output@Estimate, MR.output@Intercept),
               SE = c(MR.output@StdError.Est, MR.output@StdError.Int),
               LCI = c(MR.output@CILower.Est, MR.output@CILower.Int),
               UCI = c(MR.output@CIUpper.Est, MR.output@CIUpper.Int),
               P_value = c(MR.output@Pvalue.Est, MR.output@Pvalue.Int)
    )
  },
  error = function(e) rbind(error_df(df_mr), error_df(df_mr))
  )
  
  res2[, Method := c("MR_Egger", "EggerInt")]
  return(res2)
}

# MR IVW-PCA method
MR_PCA <- function(df_mr, ldrho, var_exp=0.99){
  res2 <- tryCatch({
    attach(df_mr)
    Phi = (beta.x / se.y) %o% (beta.x / se.y) * ldrho
    
    K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2 /
                       sum((prcomp(Phi, scale=FALSE)$sdev^2)))
              > var_exp)[1]
    
    betaXG0 = as.numeric(beta.x%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
    betaYG0 = as.numeric(beta.y%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
    
    Omega = se.y %o% se.y * ldrho
    
    pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]
    

    beta_IVWcorrel.pc <- solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
    beta_IVWcorrel.pc <- as.numeric(beta_IVWcorrel.pc)
    se_IVWcorrel.fixed.pc <- sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)) %>%
      as.numeric

    Z_score = beta_IVWcorrel.pc / se_IVWcorrel.fixed.pc
    P_value = 2*pnorm(-abs(Z_score))
    LCI <- beta_IVWcorrel.pc - qnorm(0.975)*se_IVWcorrel.fixed.pc
    UCI <- beta_IVWcorrel.pc + qnorm(0.975)*se_IVWcorrel.fixed.pc
    
    detach(df_mr)
    
    data.table(N_ins = nrow(df_mr),
               beta = beta_IVWcorrel.pc,
               SE = se_IVWcorrel.fixed.pc,
               LCI = LCI, UCI = UCI,
               P_value = P_value)
  },
  error = function(e) error_df(df_mr)
  )
  
  return (res2)
}


#run MR using above models 
run_MR_all <- function(df_mr, ldrho){
  setDT(df_mr)
  if (nrow(df_mr) == 0) {
    res2 <- blank_df()
  } else if (nrow(df_mr) == 1) {
    res2 <- single.ins.MR(df_mr)
    res2[, `:=`(Method = "Wald_1Ins")]
  } else {
    
    insSNPs <- df_mr$SNP
    
    ldrho_ins <- ldrho[insSNPs, insSNPs]
    
    # Exclude SNP if LDcorr results in NaN
    NaN_ins <- which(is.nan(ldrho_ins), T) %>% rownames
    if (!is.null(NaN_ins)){
      ins <- rownames(ldrho_ins) %>% .[which(!. %in% NaN_ins)]
      ldrho_ins <- ldrho_ins[ins, ins]
      
      df_mr <- df_mr[SNP %in% ins]
    }
    
    res_IVW <- MR_IVW(df_mr, ldrho=ldrho_ins)
    res_IVW[, Method := "IVW"]
    
    res_PCA_0.99 <- MR_PCA(df_mr, ldrho_ins, var_exp=0.99) %>%
      .[, `:=`(Method = "PCA_0.99")]
    res_PCA_0.90 <- MR_PCA(df_mr, ldrho_ins, var_exp=0.90) %>%
      .[, `:=`(Method = "PCA_0.90")]
    
    res_Egger <- MR_Egger(df_mr, ldrho=ldrho_ins)
    
    res2 <- rbind(res_IVW, res_PCA_0.99, res_PCA_0.90, res_Egger)
  }
  return(res2)
}

