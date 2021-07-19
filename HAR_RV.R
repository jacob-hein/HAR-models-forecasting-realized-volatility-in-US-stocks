# install.packages("rstudioapi")
# install.packages("rJava")
# install.packages("readxl")
# install.packages("xts")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("xlsx")
# install.packages("e1071")
# install.packages("stringr")
# install.packages("rlang")
# install.packages(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directive to current folder
options(scipen=999) # No scientific floating points notation

library(xts) # Time-indexed data frames, ideal for plots 
library(readxl) # Reading data
library(xlsx) # Exporting results
library(e1071) # Skew & Kurtosis for Summary Statistics Table
library(stringr) # String formatting
library(rlang)
library(ggplot2) # Plots
theme_set(theme_minimal()) # Plot Window
###############
# For Robust Standard Errors
library(lmtest)
library(sandwich)
###############

#Forecast, Errors, Betas & Plots
HAR <- function(data, out_sample = 96, plot_scalar = 1, extra_plots = FALSE) {
  
  RV = data$RV
  RQ = data$RQ
  BPV = data$BPV
  RV_p = data$RV_plus
  RV_m = data$RV_minus
  
  nobs = length(RV)
  in_sample = nobs - out_sample
  
  outRV = RV[(in_sample+1):(length(RV))] # We +1 to get equal length as out_sample size
  lag = 22 # 22 days lag is equivalent to one month of trading days lag
  
  all_predsA = rep(0, times = out_sample)
  all_preds = rep(0, times = out_sample)
  all_predsQ = rep(0, times = out_sample)
  all_predsF = rep(0, times = out_sample)
  all_predsC = rep(0, times = out_sample)
  all_predsS = rep(0, times = out_sample)
  all_predsJ = rep(0, times = out_sample)
  
  all_betasA = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 4)
  all_betas = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 4)
  all_betasQ = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 5)
  all_betasF = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 7)
  all_betasC = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 4)
  all_betasS = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 5)
  all_betasJ = matrix(rep(0, times = out_sample), nrow = out_sample, ncol = 5)
  
  for (t in 1:(out_sample)) {
    # Estimation
    y = RV[(lag + t + 1):(in_sample + t)]
    XA = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 3)
    X = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 3)
    XQ = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 4)
    XF = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 6)
    XC = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 3)
    XS = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 4)
    XJ = matrix(rep(0, times = in_sample-lag), nrow = in_sample-lag, ncol = 4)
    
    for (i in 0:(in_sample - lag - 1)) {
      
      # AR(3)
      XA[i+1,1] = RV[(-1+i+1+lag + t)]
      XA[i+1,2] = RV[(-2+i+1+lag + t)]
      XA[i+1,3] = RV[(-3+i+1+lag + t)]
      
      # HAR
      X[i+1,1] = RV[(-1+i+1+lag + t)]
      X[i+1,2] = (1/5)*sum(RV[(-5+i+1+lag + t):(i+lag + t)])
      X[i+1,3] = (1/22)*sum(RV[(-22+i+1+lag + t):(i+lag + t)])
      
      # HARQ
      XQ[i+1,1] = RV[(-1+i+1+lag + t)]
      XQ[i+1,2] = (1/5)*sum(RV[(-5+i+1+lag + t):(i+lag + t)])
      XQ[i+1,3] = (1/22)*sum(RV[(-22+i+1+lag + t):(i+lag + t)])
      XQ[i+1,4] = (RQ[-1+i+1+lag + t]^(1/2) * RV[(-1+i+1+lag + t)])
      
      # HARQ-F
      XF[i+1,1] = RV[(-1+i+1+lag + t)]
      XF[i+1,2] = (1/5)*sum(RV[(-5+i+1+lag + t):(i+lag + t)])
      XF[i+1,3] = (1/22)*sum(RV[(-22+i+1+lag + t):(i+lag + t)])
      XF[i+1,4] = (RQ[-1+i+1+lag + t]^(1/2) * RV[(-1+i+1+lag + t)])
      XF[i+1,5] = (((1/5)*sum(RQ[(-5+i+1+lag + t):(i+lag + t)]))^(1/2) * ((1/5)*sum(RV[(-5+i+1+lag + t):(i+lag + t)])))
      XF[i+1,6] = (((1/22)*sum(RQ[(-22+i+1+lag + t):(i+lag + t)]))^(1/2) * ((1/22)*sum(RV[(-22+i+1+lag + t):(i+lag + t)])))
      
      # CHAR
      XC[i+1,1] = BPV[(-1+i+1+lag + t)]
      XC[i+1,2] = (1/5)*sum(BPV[(-5+i+1+lag + t):(i+lag + t)])
      XC[i+1,3] = (1/22)*sum(BPV[(-22+i+1+lag + t):(i+lag + t)])
      
      # SHAR
      XS[i+1,1] = (1/5)*sum(RV[(-5+i+1+lag + t):(i+lag + t)])
      XS[i+1,2] = (1/22)*sum(RV[(-22+i+1+lag + t):(i+lag + t)])
      XS[i+1,3] = RV_p[(-1+i+1+lag + t)]
      XS[i+1,4] = RV_m[(-1+i+1+lag + t)]
      
      # HAR-J
      XJ[i+1,1] = RV[(-1+i+1+lag + t)]
      XJ[i+1,2] = (1/5)*sum(RV[(-5+i+1+lag + t):(i+lag + t)])
      XJ[i+1,3] = (1/22)*sum(RV[(-22+i+1+lag + t):(i+lag + t)])
      XJ[i+1,4] = max((RV[(-1+i+1+lag + t)] - BPV[(-1+i+1+lag + t)]), 0)
      
    }
    
    
    # R Regression at t=1 for Standard Errors before performing any out-of-sample forecasts
    if (t==out_sample) {
      modelA = lm(y ~ XA)
      model = lm(y ~ X)
      modelQ = lm(y ~ XQ)
      modelF = lm(y ~ XF)
      modelC = lm(y ~ XC)
      modelS = lm(y ~ XS)
      # HAR-J if statement:
      if (sum(XJ[,4]) ==0) {
        XJ[1,4] = 0.1 # Ensure invertibility if XJ singular with zero column
      }
      modelJ = lm(y ~ XJ)
      models_at_t_1 = list("modelA" = modelA, "model" = model, 
                           "modelQ" = modelQ, "modelF" = modelF, 
                           "modelC" = modelC, "modelS" = modelS, 
                           "modelJ" = modelJ)
      
      # Below we retrieve R^2 & Adjusted R^2, prior to out-of-sample forecasts
      num_of_models = 7
      r_squareds = matrix(0, nrow = 2, ncol = num_of_models)
      for (val in 1:length(models_at_t_1)) {
        r_squareds[1,val] = summary(models_at_t_1[[val]])$r.squared
        r_squareds[2,val] = summary(models_at_t_1[[val]])$adj.r.squared
      }
      rownames(r_squareds) = c("R-squared", "Adj.R-squared")
      colnames(r_squareds) = c("AR(3)", "HAR", "HARQ", "HARQ-F", "CHAR", "SHAR", "HAR-J")
      
    }
    
    XA = cbind(rep(1, times = nrow(XA)), XA)
    X = cbind(rep(1, times = nrow(X)), X)
    XQ = cbind(rep(1, times = nrow(XQ)), XQ)
    XF = cbind(rep(1, times = nrow(XF)), XF)
    XC = cbind(rep(1, times = nrow(XC)), XC)
    XS = cbind(rep(1, times = nrow(XS)), XS)
    XJ = cbind(rep(1, times = nrow(XJ)), XJ)
    
    # HAR-J if statement:
    if (sum(XJ[,5]) ==0) {
      XJ[1,5] = 0.1 # Ensure invertibility if XJ singular with zero column
    }
    
    # OLS Regression
    betasA = solve(t(XA) %*% XA) %*% t(XA) %*% y
    betas = solve(t(X) %*% X) %*% t(X) %*% y
    betasQ = solve(t(XQ) %*% XQ) %*% t(XQ) %*% y
    betasF = solve(t(XF) %*% XF) %*% t(XF) %*% y
    betasC = solve(t(XC) %*% XC) %*% t(XC) %*% y
    betasS = solve(t(XS) %*% XS) %*% t(XS) %*% y
    betasJ = solve(t(XJ) %*% XJ) %*% t(XJ) %*% y
    
    b0A = betasA[1]
    b1A = betasA[2]
    b2A = betasA[3]
    b3A = betasA[4]
    
    b0 = betas[1]
    b1 = betas[2]
    b2 = betas[3]
    b3 = betas[4]
    
    b0Q = betasQ[1]
    b1Q = betasQ[2]
    b2Q = betasQ[3]
    b3Q = betasQ[4]
    b1Q_Q = betasQ[5]
    
    b0F = betasF[1]
    b1F = betasF[2]
    b2F = betasF[3]
    b3F = betasF[4]
    b1F_Q = betasF[5]
    b2F_Q = betasF[6]
    b3F_Q = betasF[7]
    
    b0C = betasC[1]
    b1C = betasC[2]
    b2C = betasC[3]
    b3C = betasC[4]
    
    b0S = betasS[1]
    b1S_P = betasS[2]
    b2S = betasS[3]
    b3S = betasS[4]
    b1S_M = betasS[5]
    
    b0J = betasJ[1]
    b1J = betasJ[2]
    b2J = betasJ[3]
    b3J = betasJ[4]
    bJ = betasJ[5]
    
    all_betasA[t,1] = b0A
    all_betasA[t,2] = b1A
    all_betasA[t,3] = b2A
    all_betasA[t,4] = b3A
    
    all_betas[t,1] = b0
    all_betas[t,2] = b1
    all_betas[t,3] = b2
    all_betas[t,4] = b3
    
    all_betasQ[t,1] = b0Q
    all_betasQ[t,2] = b1Q
    all_betasQ[t,3] = b2Q
    all_betasQ[t,4] = b3Q
    all_betasQ[t,5] = b1Q_Q
    
    all_betasF[t,1] = b0F
    all_betasF[t,2] = b1F
    all_betasF[t,3] = b2F
    all_betasF[t,4] = b3F
    all_betasF[t,5] = b1F_Q
    all_betasF[t,6] = b2F_Q
    all_betasF[t,7] = b3F_Q
    
    all_betasC[t,1] = b0C
    all_betasC[t,2] = b1C
    all_betasC[t,3] = b2C
    all_betasC[t,4] = b3C
    
    all_betasS[t,1] = b0S
    all_betasS[t,2] = b1S_P
    all_betasS[t,3] = b2S
    all_betasS[t,4] = b3S
    all_betasS[t,5] = b1S_M
    
    all_betasJ[t,1] = b0J
    all_betasJ[t,2] = b1J
    all_betasJ[t,3] = b2J
    all_betasJ[t,4] = b3J
    all_betasJ[t,5] = bJ
    
    
    
    # Prediction at time-step t
    predA = b0A + b1A*XA[nrow(XA),2] + b2A*XA[nrow(XA),3] + b3A*XA[nrow(XA),4]
    pred = b0 + b1*X[nrow(X),2] + b2*X[nrow(X),3] + b3*X[nrow(X),4]
    predQ = b0Q + b1Q*XQ[nrow(XQ),2] + b2Q*XQ[nrow(XQ),3] + b3Q*XQ[nrow(XQ),4] + b1Q_Q*XQ[nrow(XQ),5]
    predF = b0F + b1F*XF[nrow(XF),2] + b2F*XF[nrow(XF),3] + b3F*XF[nrow(XF),4] + b1F_Q*XF[nrow(XF),5] + b2F_Q*XF[nrow(XF),6] + b3F_Q*XF[nrow(XF),7]
    predC = b0C + b1C*XC[nrow(XC),2] + b2C*XC[nrow(XC),3] + b3C*XC[nrow(XC),4]
    predS = b0S + b1S_P*XS[nrow(XS),2] + b2S*XS[nrow(XS),3] + b3S*XS[nrow(XS),4] + b1S_M*XS[nrow(XS),5]
    predJ = b0J + b1J*XJ[nrow(XJ),2] + b2J*XJ[nrow(XJ),3] + b3J*XJ[nrow(XJ),4] + bJ*XJ[nrow(XJ),5]
    
    # Saving time-step t prediction within list for error computations
    all_predsA[t] = predA
    all_preds[t] = pred
    all_predsQ[t] = predQ
    all_predsF[t] = predF
    all_predsC[t] = predC
    all_predsS[t] = predS
    all_predsJ[t] = predJ
    
  }
  
  # Error Computations: Mean Squared Error & Mean Absolute Error
  AR_mse = mean((outRV - all_predsA)^2)
  AR_mae = mean(abs(outRV - all_predsA))
  
  HAR_mse = mean((outRV - all_preds)^2)
  HAR_mae = mean(abs(outRV - all_preds))
  
  HARQ_mse = mean((outRV - all_predsQ)^2)
  HARQ_mae = mean(abs(outRV - all_predsQ))
  
  HARQF_mse = mean((outRV - all_predsF)^2)
  HARQF_mae = mean(abs(outRV - all_predsF))
  
  CHAR_mse = mean((outRV - all_predsC)^2)
  CHAR_mae = mean(abs(outRV - all_predsC))
  
  SHAR_mse = mean((outRV - all_predsS)^2)
  SHAR_mae = mean(abs(outRV - all_predsS))
  
  HARJ_mse = mean((outRV - all_predsJ)^2)
  HARJ_mae = mean(abs(outRV - all_predsJ))
  
  
  errors_list = list("AR_mse" = AR_mse, "AR_mae" = AR_mae, "HAR_mse" = HAR_mse, "HAR_mae" = HAR_mae,
                "HARQ_mse" = HARQ_mse, "HARQ_mae" = HARQ_mae, "HARQF_mse" = HARQF_mse, "HARQF_mae" = HARQF_mae,
                "CHAR_mse" = CHAR_mse, "CHAR_mae" = CHAR_mae, "SHAR_mse" = SHAR_mse, "SHAR_mae" = SHAR_mae,
                "HARJ_mse" = HARJ_mse, "HARJ_mae" = HARJ_mae)
  
  # Plots
  #plot_scalar = 100
  if (plot_scalar == 1) {
    ylab_str = "Realized Variance"
  } else {
    ylab_str = paste("Realized Var. (Axis scaled by a factor of ", plot_scalar, ")", sep="")
  }
  
  plotDates = DataSet$Dates[(length(DataSet$Dates)-out_sample+1):length(DataSet$Dates)]
  
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsA ~ plotDates, col="darkred")
  legend("topright", legend=c("Actual RV", "AR(3)"), col=c("blue", "darkred"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/AR(3) R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  #plot_path_png = paste(substr(plot_path, 1, nchar(plot_path)-4), ".png", sep="")
  dev.off()
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_preds ~ plotDates, col="red")
  legend("topright", legend=c("Actual RV", "HAR"), col=c("blue", "red"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/HAR R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsQ ~ plotDates, col="green")
  legend("topright", legend=c("Actual RV", "HARQ"), col=c("blue", "green"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/HARQ R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsF ~ plotDates, col="darkgreen")
  legend("topright", legend=c("Actual RV", "HARQ-F"), col=c("blue", "darkgreen"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/HARQ-F R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsC ~ plotDates, col="black")
  legend("topright", legend=c("Actual RV", "CHAR"), col=c("blue", "black"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/CHAR R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsS ~ plotDates, col="maroon")
  legend("topright", legend=c("Actual RV", "SHAR"), col=c("blue", "maroon"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/SHAR R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsJ ~ plotDates, col="orange")
  legend("topright", legend=c("Actual RV", "HAR-J"), col=c("blue", "orange"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1)
  plot_path = paste(getwd(), "/OutSampPlots/HAR-J R Plot.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  
  if (extra_plots == TRUE) {
    # Plotting All models in one 
    plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="All Models ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_predsA ~ plotDates, col="darkred")
    #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HAR ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_preds ~ plotDates, col="red")
    #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HARQ ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_predsQ ~ plotDates, col="green")
    #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HARQ-F ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_predsF ~ plotDates, col="darkgreen")
    #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="CHAR ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_predsC ~ plotDates, col="orange")
    #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="SHAR ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_predsS ~ plotDates, col="maroon")
    #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HAR-J ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
    lines(plot_scalar*all_predsJ ~ plotDates, col="black")
    plot_path = paste(getwd(), "/OutSampPlots/AllModels R Plot.pdf", sep="")
    dev.copy(pdf, plot_path)
    dev.off()
    
  }
  
  par(mfrow=c(2,2))
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  #lines(plot_scalar*all_predsA ~ plotDates, col="darkred")
  #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HAR ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_preds ~ plotDates, col="red")
  legend("topright", legend=c("Actual RV", "HAR"), col=c("blue", "red"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsJ ~ plotDates, col="orange")
  legend("topright", legend=c("Actual RV", "HAR-J"), col=c("blue", "orange"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HARQ ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsQ ~ plotDates, col="green")
  legend("topright", legend=c("Actual RV", "HARQ"), col=c("blue", "green"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
 
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsF ~ plotDates, col="darkgreen")
  legend("topright", legend=c("Actual RV", "HARQ-F"), col=c("blue", "darkgreen"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="CHAR ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  #lines(plot_scalar*all_predsC ~ plotDates, col="orange")
  #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="SHAR ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  #lines(plot_scalar*all_predsS ~ plotDates, col="maroon")
  #plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", main="HAR-J ", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  plot_path = paste(getwd(), "/OutSampPlots/FOUR_R_Plots.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off()
  
  par(mfrow=c(2,2))
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsC ~ plotDates, col="black")
  legend("topright", legend=c("Actual RV", "CHAR"), col=c("blue", "black"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsS ~ plotDates, col="maroon")
  legend("topright", legend=c("Actual RV", "SHAR"), col=c("blue", "maroon"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsC ~ plotDates, col="black")
  legend("topright", legend=c("Actual RV", "CHAR"), col=c("blue", "black"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
  plot(plot_scalar*outRV ~ plotDates, type="l", col="blue", xlab = "Out-of-Sample Trading Days", ylab = ylab_str)
  lines(plot_scalar*all_predsS ~ plotDates, col="maroon")
  legend("topright", legend=c("Actual RV", "SHAR"), col=c("blue", "maroon"), lty=c(1,1), lwd=c(2.5,2.5), box.lty=1, bty = "n")
  
  plot_path = paste(getwd(), "/OutSampPlots/TWO_R_Plots.pdf", sep="")
  dev.copy(pdf, plot_path)
  dev.off() 
  
  
  # Output formatting: 
  output = matrix(outRV)
  out_sample_dates = data$Date[(in_sample+1):length(data$Date)] 
  output = cbind(out_sample_dates, output, all_predsA, all_preds, all_predsQ, all_predsF, all_predsC, all_predsS, all_predsJ)
  output_df = as.data.frame(output)
  colnames(output_df) = c("Date", "outRV", "all_predsA", "all_preds", "all_predsQ", "all_predsF", "all_predsC", "all_predsS", "all_predsJ")
  betas_list = list("all_betasA" = all_betasA, "all_betas" = all_betas, "all_betasQ" = all_betasQ, 
                    "all_betasF" = all_betasF, "all_betasC" = all_betasC, "all_betasS" = all_betasS, 
                    "all_betasJ" = all_betasJ)
  
  output_df_errors_betas = list(output_df, errors_list, betas_list, models_at_t_1, r_squareds)
  
  return(output_df_errors_betas)
}

#RV, RQ, BPV Estimator
estimator <- function(data) {
  RV_t_estimates = c()
  RQ_t_estimates = c()
  BPV_t_estimates = c()
  RV_t_plus_estimates = c()
  RV_t_minus_estimates = c()
  RV_t_dates = c()
  M = 0 # Intraday obs used in estimation of RV_t
  
  for (t in 1:length(data$Dates)) {
    # t accounts for the final number of daily RV_t estimates
    RV_t_i_estimates = c() # M number of r_t,i to be summed up
    M_Q = M # Counter
    
    while (substring(data$Dates[t+M], first = 1, last = 5) == substring(data$Dates[t+M+1], first = 1, last = 5)
           # The below AND condition breaks while-loop when no more intraday obs available
           & !is.na(data$Dates[t+M+1])) {
      
      # Intraday returns 
      RV_t_i = (data$Open[t+M+1] - data$Open[t+M])
      RV_t_i_estimates = c(RV_t_i_estimates, RV_t_i)
      
      M = M + 1
    }
    
    if (is.na(data$Dates[t+M])) {
      break # This if-clause breaks for-loop when the eventual NA intraday obs is reached 
    }
    
    RV_t = sum(RV_t_i_estimates^2) #Realized Variance
    #RV_t = sqrt(RV_t) #Realized Volatilty
    
    RQ_t = ((M-M_Q)/3) * sum(RV_t_i_estimates^4)
    #RQ_t = sqrt(RQ_t)
    
    # Bi-Power Variance: |r_t,i||r_t,i+1|
    
    ############### Faulty one line of code below:
    # Notice [-length(RV_....)] ensures final element is left out.
    # BPV_t = (sqrt(2/pi))^(-2) * (sum(abs(RV_t_i_estimates[-length(RV_t_i_estimates)]))*sum(abs(RV_t_i_estimates)))
    ############### Line above discarded as wrong implementation
    
    BPV_t_i_estimates = c()
    # i in 1:len(...)-1 corresponds to summing up to M-1 as in Bollerslev (2016)
    for (i in 1:(length(RV_t_i_estimates)-1)) {
      BPV_t_i_estimates = c(BPV_t_i_estimates, abs(RV_t_i_estimates[i] * abs(RV_t_i_estimates[i+1])))
    }
    BPV_t = (sqrt(2/pi))^(-2) * sum(BPV_t_i_estimates)
    #BPV = sqrt(BPV_t)
    
    # RV Plus and RV Minus for SHAR model spec
    RV_t_plus = sum(RV_t_i_estimates[RV_t_i_estimates > 0]^2)
    #RV_t_plus = sqrt(RV_t_plus)
    RV_t_minus = sum(RV_t_i_estimates[RV_t_i_estimates < 0]^2)
    #RV_t_minus = (-1)*sqrt(abs(RV_t_minus)) # Abs value and multiply by -1 to avoid sqrt'ing negative numbers
    
    RV_t_estimates = c(RV_t_estimates, RV_t)
    RQ_t_estimates = c(RQ_t_estimates, RQ_t)
    BPV_t_estimates = c(BPV_t_estimates, BPV_t)
    RV_t_plus_estimates = c(RV_t_plus_estimates, RV_t_plus)
    RV_t_minus_estimates = c(RV_t_minus_estimates, RV_t_minus)
    
    # Dates
    RV_t_date = as.numeric(substring(data$Dates[t+M], first = 1, last = 5))
    RV_t_dates = c(RV_t_dates, RV_t_date)
  }
  
  RV_df = as.data.frame(RV_t_dates)
  RV_df = cbind(RV_df, RV_t_estimates, RQ_t_estimates, BPV_t_estimates, RV_t_plus_estimates, RV_t_minus_estimates)
  colnames(RV_df) = c("Dates", "RV", "RQ", "BPV", "RV_plus", "RV_minus")
  RV_df$Dates = as.Date(RV_df$Dates, origin = "1899-12-30")
  return(RV_df)
}

# Scaling of RV, RQ, BPV, RVplus, RVminus estimates for numerical stability
DataSet_Scalar <- function(DataSet, scalar = 100000){
  for (col in 2:ncol(DataSet)) {
    DataSet[,col] = DataSet[,col]*scalar
  }
  return(DataSet)
}

# Constructing Beta Table
betaTable <- function(forecast) {
  
  robustStdErrs_A = coeftest(forecast[[4]]$modelA, vcov = vcovHC(forecast[[4]]$modelA, type="HC1"))
  robustStdErrs_ = coeftest(forecast[[4]]$model, vcov = vcovHC(forecast[[4]]$model, type="HC1"))
  robustStdErrs_Q = coeftest(forecast[[4]]$modelQ, vcov = vcovHC(forecast[[4]]$modelQ, type="HC1"))
  robustStdErrs_F = coeftest(forecast[[4]]$modelF, vcov = vcovHC(forecast[[4]]$modelF, type="HC1"))
  robustStdErrs_C = coeftest(forecast[[4]]$modelC, vcov = vcovHC(forecast[[4]]$modelC, type="HC1"))
  robustStdErrs_S = coeftest(forecast[[4]]$modelS, vcov = vcovHC(forecast[[4]]$modelS, type="HC1"))
  robustStdErrs_J = coeftest(forecast[[4]]$modelJ, vcov = vcovHC(forecast[[4]]$modelJ, type="HC1"))
  
  
  
  Robust_T_test = list("modelA" = robustStdErrs_A, "model" = robustStdErrs_,
                       "modelQ" = robustStdErrs_Q, "modelF" = robustStdErrs_F,
                       "modelC" = robustStdErrs_C, "modelS" = robustStdErrs_S,
                       "modelJ" = robustStdErrs_J)
  
  Robust_T_test_matrix = do.call("rbind", Robust_T_test)[, 1:2]
  Robust_T_test_table = matrix(NA, nrow = 10, ncol = 7)
  colnames(Robust_T_test_table) = c("AR(3)", "HAR", "HARQ", "HARQ-F", "CHAR", "SHAR", "HAR-J")
  rownames(Robust_T_test_table) = c("b0", "b1", "b2", "b3", "b1Q", "b2Q", "b3Q","b1+", "b1-", "bJ")
  
  # c = Counter that counts through the rows of #Robust_T_test_matrix to extract the needed betas
  # l = length of beta parameters of given model
  # dig = num of digits to round up estimates
  dig = 4
  c_AR = 1:length(forecast[[4]]$modelA$coefficients)
  Robust_T_test_table[c_AR , 1] = paste(round(Robust_T_test_matrix[c_AR, 1], dig), " (",round(Robust_T_test_matrix[c_AR, 2], dig),")", sep="")
  
  l_HAR = 1:length(forecast[[4]]$model$coefficients)
  c_HAR = length(c_AR) + l_HAR
  Robust_T_test_table[l_HAR , 2] = paste(round(Robust_T_test_matrix[c_HAR, 1], dig), " (",round(Robust_T_test_matrix[c_HAR, 2], dig),")", sep="")
  
  l_HARQ = 1:length(forecast[[4]]$modelQ$coefficients)
  c_HARQ = length(c_HAR) + l_HARQ + 4 # Ad-Hoc added values found by checking Robust_T_test_matrix[c_HAR, ] throughout
  Robust_T_test_table[l_HARQ , 3] = paste(round(Robust_T_test_matrix[c_HARQ, 1], dig), " (",round(Robust_T_test_matrix[c_HARQ, 2], dig),")", sep="")
  
  l_HARQF = 1:length(forecast[[4]]$modelF$coefficients)
  c_HARQF = length(c_HARQ) + l_HARQF + 8
  Robust_T_test_table[l_HARQF , 4] = paste(round(Robust_T_test_matrix[c_HARQF, 1], dig), " (",round(Robust_T_test_matrix[c_HARQF, 2], dig),")", sep="")
  
  l_CHAR = 1:length(forecast[[4]]$modelC$coefficients)
  c_CHAR = length(l_HARQF) + l_CHAR + 13
  Robust_T_test_table[l_CHAR , 5] = paste(round(Robust_T_test_matrix[c_CHAR, 1], dig), " (",round(Robust_T_test_matrix[c_CHAR, 2], dig),")", sep="")
  
  l_SHAR = 1:length(forecast[[4]]$modelS$coefficients)
  c_SHAR = length(l_CHAR) + l_SHAR + 20
  Robust_T_test_table[c(1,3,4,8,9) , 6] = paste(round(Robust_T_test_matrix[c_SHAR, 1], dig), " (",round(Robust_T_test_matrix[c_SHAR, 2], dig),")", sep="")
  
  l_HARJ = 1:length(forecast[[4]]$modelJ$coefficients)
  c_HARJ = length(l_SHAR) + l_HARJ + 24
  Robust_T_test_table[c(1:4,10) , 7] = paste(round(Robust_T_test_matrix[c_HARJ, 1], dig), " (",round(Robust_T_test_matrix[c_HARJ, 2], dig),")", sep="")
  
  #print(Robust_T_test_matrix)
  #print(Robust_T_test_table)
  
  # Appending R^squareds to table
  Robust_T_test_table_w_rsquareds = rbind(Robust_T_test_table, NA*c(1:7), round(forecast[[5]][1,], dig), round(forecast[[5]][2,], dig))
  rownames(Robust_T_test_table_w_rsquareds)[c(12,13)] = c("R^2", "Adj.R^2")
  # Exporting Results
  write.xlsx(Robust_T_test_table_w_rsquareds, paste(getwd(), "/Results/Betas_w_all_rsquareds.xlsx", sep=""), sheetName="Sheet1", 
             col.names=TRUE, row.names=TRUE, append=FALSE)
  
  return(Robust_T_test_table_w_rsquareds)
}

# Summary Stats
summaryStats <- function(stocks, scalar = 1, freq = "5min_extended") {
  sumStats = matrix(NA, nrow = length(stocks), ncol = 7)
  colnames(sumStats) = c("Min", "Mean", "Median", "Max", "Std. Dev.", "Skewness", "Kurtosis")
  Symbol = c()
  for (stockname in 1:length(stocks)) {
    Symbol = c(Symbol,substr(stocks[stockname],1, (nchar(stocks[stockname])-10)))
  }
  rownames(sumStats) = Symbol
  
  # Realized Measures of all in assets in 'stocks' variable 
  DataSets = list(rep(NA, times = length(stocks)))
  for (stock in 1:length(stocks)) {
    excel_file = paste("Data/",freq,".xlsx" , sep="")
    data_name = stocks[stock]
    data = as.data.frame(read_excel(excel_file, sheet = data_name))
    data = data[c("BarTp", "Trade")]
    colnames(data) = c("Dates", "Open")
    data = data[-c(1:4),]
    data$Open = as.numeric(data$Open)
    data_log = cbind(data$Dates, as.data.frame(log(data$Open)))
    colnames(data_log) = c("Dates", "Open")
    
    #DataSet = estimator(data)
    DataSet = estimator(data_log)
    DataSet = DataSet_Scalar(DataSet) # DataSet Scalar Function!
    
    # Appending all data.frames of RV measures to a list of data.frames
    DataSets[stock] = list(DataSet)
  }
  names(DataSets) = stocks
  
  # Summary Stats to summStats table. Scaled by factor of variable 'scalar' 
  for (stk in 1:length(stocks)) {
    curRV = DataSets[[stk]]$RV # Realized Vol of current stock in iteration
    sumStats[stk, 1] = round(min(curRV)*scalar,3) 
    sumStats[stk, 2] = round(mean(curRV)*scalar,3)
    sumStats[stk, 3] = round(median(curRV)*scalar,3)
    sumStats[stk, 4] = round(max(curRV)*scalar ,3)
    sumStats[stk, 5] = round(sd(curRV)*scalar,3)
    sumStats[stk, 6] = round(skewness(curRV)*scalar,1)
    sumStats[stk, 7] = round(kurtosis(curRV)*scalar,1) 
  }
  sumStats = cbind(Symbol, sumStats)
  Symbol_str = str_c(Symbol,collapse='_') #stringr function dependency
  sumStats_file_path = paste(getwd(), "/Results/SummaryStats_", Symbol_str, "_", freq, ".xlsx", sep="")
  
  # Exporting Summary Stats
  write.xlsx(sumStats, sumStats_file_path, sheetName="Sheet1",
             col.names=TRUE, row.names=TRUE, append=FALSE)
  print(paste("All values in the Summary Statistics are scaled by a factor of: ", scalar, sep=""))
  return(sumStats)
}

# Error Computation & Table Construction
errorsTable <- function(stocks, dig = 5, errorsScalar = 1) {
  errorsMatrix = matrix(NA, nrow = length(stocks)+1, ncol = 14)
  Symbol = c()
  for (stockname in 1:length(stocks)) {
    Symbol = c(Symbol,substr(stocks[stockname],1, (nchar(stocks[stockname])-10)))
  }
  rownames(errorsMatrix) = c("", Symbol)
  errorsMatrix[1, ] = rep(c("MSE", "MAE"), times = 7)
  
  model_names = c("AR(3)", "HAR", "HARQ", "HARQ-F", "CHAR", "SHAR", "HAR-J")
  model_names2 = rep(NA, times = ncol(errorsMatrix))
  
  #odds = odds[lapply(odds, "%%", 2) != 0]
  mod = 0
  for (m in 1:length(model_names)) {
    model_names2[(mod+m)] = model_names[m]
    mod = mod + 1
  }
  for (m in 1:length(model_names2)) {
    if(is.na(model_names2[m])) {
      model_names2[m] = model_names2[(m-1)]
    }
  }
  colnames(errorsMatrix) = model_names2
  
  
  # Realized Measures of all in assets in 'stocks' variable 
  forecasts = list(rep(NA, times = length(stocks)))
  for (stock in 1:length(stocks)) {
    excel_file = paste("Data/",freq,".xlsx" , sep="")
    data_name = stocks[stock]
    data = as.data.frame(read_excel(excel_file, sheet = data_name))
    data = data[c("BarTp", "Trade")]
    colnames(data) = c("Dates", "Open")
    data = data[-c(1:4),]
    data$Open = as.numeric(data$Open)
    data_log = cbind(data$Dates, as.data.frame(log(data$Open)))
    colnames(data_log) = c("Dates", "Open")
    
    #DataSet = estimator(data)
    DataSet = estimator(data_log)
    DataSet = DataSet_Scalar(DataSet)
    forecast = HAR(DataSet)
    
    # Appending all data.frames of RV measures to a list of data.frames
    forecasts[stock] = list(forecast)
  }
  names(forecasts) = stocks
  
  #errorsScalar = 100
  #dig = 5
  for (s in 1:length(stocks)) {
    # Loop over errors and models
    all_errors = forecasts[[s]][[2]]
    for (e in 1:length(all_errors)) {
      errorsMatrix[(s+1), e] = round(all_errors[[e]]*errorsScalar, dig)
    }
  }
  
  errorsMatrix_BM = duplicate(errorsMatrix)
  HAR_MSE_idx = grep("^HAR$", colnames(errorsMatrix_BM))[1]
  HAR_MAE_idx = grep("^HAR$", colnames(errorsMatrix_BM))[2]
  
  #as.numeric(errorsMatrix_BM[1:length(stocks)+1, HAR_MSE_idx])
  
  MSE_cols = (1:ncol(errorsMatrix_BM))[1:ncol(errorsMatrix_BM) %% 2 != 0]
  MAE_cols = (1:ncol(errorsMatrix_BM))[1:ncol(errorsMatrix_BM) %% 2 == 0]
  
  for (e in 1:ncol(errorsMatrix_BM)) {
    if (e %in% MSE_cols) {
      errorsMatrix_BM[1:length(stocks)+1,e] = as.numeric(errorsMatrix[1:length(stocks)+1,e]) / as.numeric(errorsMatrix[1:length(stocks)+1, HAR_MSE_idx]) # MSE_idx
    } else if (e %in% MAE_cols) {
      errorsMatrix_BM[1:length(stocks)+1,e] = as.numeric(errorsMatrix[1:length(stocks)+1,e]) / as.numeric(errorsMatrix[1:length(stocks)+1, HAR_MAE_idx]) #MAE_idx
    }
  }
  
  #Unbenchmarked Errors
  errorsMatrix_file_path = paste(getwd(), "/Results/Errors_UnBenchmarked.xlsx", sep="")
  errorsMatrix_BM_file_path = paste(getwd(), "/Results/Errors_Benchmarked_to_HAR.xlsx", sep="")
  # Exporting Errors Table
  write.xlsx(errorsMatrix, errorsMatrix_file_path, sheetName="Sheet1",
             col.names=TRUE, row.names=TRUE, append=FALSE)
  write.xlsx(errorsMatrix_BM, errorsMatrix_BM_file_path, sheetName="Sheet1",
             col.names=TRUE, row.names=TRUE, append=FALSE)
  
  errorsTables = c(errorsMatrix, errorsMatrix_BM)
  return(errorsTables)
}


######## SETTINGS ########
stocks = c("SPY US Equity", "MSFT US Equity", "MCD US Equity", "JPM US Equity", "DIS US Equity")
freq = "5min_extended"
##########################

# Data Prep
log_returns_scalar = 1
excel_file = paste("Data/",freq,".xlsx" , sep="")
data_name = stocks[1]
data = as.data.frame(read_excel(excel_file, sheet = data_name))
data = data[c("BarTp", "Trade")]
colnames(data) = c("Dates", "Open")
data = data[-c(1:4),]
data$Open = as.numeric(data$Open)
data_log = cbind(data$Dates, as.data.frame(log(data$Open)*log_returns_scalar))
colnames(data_log) = c("Dates", "Open")

#DataSet = estimator(data)
DataSet = estimator(data_log)
DataSet = DataSet_Scalar(DataSet)
forecast = HAR(DataSet, out_sample = 96, extra_plots = TRUE)


#### FINAL SECTION: Generating tables, retrieving out-of-sample trading dates for plot, etc ######

###### This section generates the data exports for the tables in the paper
# beta_table = betaTable(forecast)
#summStats = summaryStats(stocks)
#errorsTables = errorsTable(stocks)

# Function for plotDates (out-of-sample forecasting)
# plotDates = substring(data$Dates, first = 1, last = 5)
# plotDates = unique(plotDates)[which(!is.na(unique(plotDates)))]
# plotDates = as.Date(as.numeric(plotDates), origin = "1899-12-30")
# length(plotDates)

# Realized Quarticity Plots of the 5 stocks
# for (i in 1:5) {
#   data_name = stocks[i]
#   data = as.data.frame(read_excel(excel_file, sheet = data_name))
#   data = data[c("BarTp", "Trade")]
#   colnames(data) = c("Dates", "Open")
#   data = data[-c(1:4),]
#   data$Open = as.numeric(data$Open)
#   data_log = cbind(data$Dates, as.data.frame(log(data$Open)*log_returns_scalar))
#   colnames(data_log) = c("Dates", "Open")
#   
#   #DataSet = estimator(data)
#   DataSet = estimator(data_log)
#   DataSet = DataSet_Scalar(DataSet)
#   plot(DataSet$RQ, type ="l", main = data_name)
#   print(data_name)
#   print(summary(DataSet$RQ))
# }

