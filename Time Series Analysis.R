library(TSA)
library(tseries)
library(readr)
library(lmtest)
library(dplyr)
library(lubridate)
library(forecast)
# install.packages("fBasics")
# library(fBasics)
setwd("F:/Subbu/RMIT/sem 3/Time Series Analysis/assign 2")
assign_ds <- read_csv("sunspots.csv",col_names = TRUE)

# str(assign_ds$Date)

# class(assign_ds)


#FILTER
assign_ds <- tail(assign_ds,841)
assign_ds

# RESIDUAL ANALYSIS FUNCTION

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH", "fGARCH")[1]){
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else if (class == "fGARCH"){
    res.model = model@residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}

# Monthly Mean Total Sunspot Number

sunspotsts <- ts(assign_ds$`Monthly Mean Total Sunspot Number` ,start=c(1951, 1), end=c(2021, 1), frequency = 12)

plot(sunspotsts, ylab = "Sunspots",  type = "o", main = "Time series plot for the Mean number of Sunspots along the years")
#time <-  time(ozonethicknessTS)

# AUTO-CORRELATION

plot(y=sunspotsts,x=zlag(sunspotsts),xlab = "Previous Year Number of sunspots",
     ylab= "Average number of sunspots", main = "Scatter plot of Average number of sunspots vs the first lag ")

# ACF AND PACF

par(mfrow = c(1,2))
acf(sunspotsts, lag.max = 100, main = "ACF of the sunspot time series")
pacf(sunspotsts, lag.max = 100, main = "PACF of the sunspot time series")
par(mfrow = c(1,1))

# UNIT ROOT TEST OF NON-STATIONARITY

adf.test(sunspotsts)

# BOX-COX TRANSFORMATION

sunspots_BC = sunspotsts+abs(min(sunspotsts))+1
BC_sunspots = BoxCox.ar(y = sunspots_BC)
BC_sunspots
lambda <- BC_sunspots$lambda[which(max(BC_sunspots$loglike) == BC_sunspots$loglike)]
lambda

sunspotsTSBC <- ((sunspots_BC^lambda) - 1) / lambda

par(mfrow = c(1,1))

# TIME SERIES PLOT OF BOX-COX TRANSFORMED SERIES
plot(sunspotsTSBC)

# CHECKING NON-STATIONARITY

adf.test(sunspotsTSBC)
acf(sunspotsTSBC)
pacf(sunspotsTSBC)

# Differencing

par(mfrow = c(1,1))

# First Difference
sunspotsTS_diff = diff(sunspotsTSBC,differences = 1)
plot(sunspotsTS_diff)

adf.test(sunspotsTS_diff)
pp.test(sunspotsTS_diff)

# series is stationary

par(mfrow = c(1,2))

acf(sunspotsTS_diff)
pacf(sunspotsTS_diff)

par(mfrow = c(1,1))

# Model Specification
# ACF AND PACF

par(mfrow = c(1,2))
acf(sunspotsTS_diff)
pacf(sunspotsTS_diff)
par(mfrow = c(1,1))
# ARIMA(3,1,2)

#EACF
eacf(sunspotsTS_diff)

# ARIMA(0,1,2) , ARIMA(0,1,3), ARIMA(1,1,3)

#BIC TABLE
sunspots_bic = armasubsets(y=sunspotsTS_diff , nar=5 , nma=5, y.name='Average number of sunspots', ar.method='ols')
plot(sunspots_bic)

par(mfrow = c(1,1))

# ARIMA(1,1,1),ARIMA(2,1,1)

#POSSIBLE MODELS

# ARIMA(3,1,2), ARIMA(0,1,2) , ARIMA(0,1,3) , ARIMA(1,1,3), ARIMA(2,1,1), ARIMA(1,1,1)

#FITTING THE MODELS

model.312 = arima(sunspotsTS_diff,order=c(3,1,2),method='CSS')
coeftest(model.312)

model.312_ML = arima(sunspotsTS_diff,order=c(3,1,2),method='ML')
coeftest(model.312_ML)

model.012 = arima(sunspotsTS_diff,order=c(0,1,2),method='CSS')
coeftest(model.012)

model.012_ML = arima(sunspotsTS_diff,order=c(0,1,2),method='ML')
coeftest(model.012_ML)

model.013 = arima(sunspotsTS_diff,order=c(0,1,3),method='CSS')
coeftest(model.013)

model.013_ML = arima(sunspotsTS_diff,order=c(0,1,3),method='ML')
coeftest(model.013_ML)

model.113 = arima(sunspotsTS_diff,order=c(1,1,3),method='CSS')
coeftest(model.113)

model.113_ML = arima(sunspotsTS_diff,order=c(1,1,3),method='ML')
coeftest(model.113_ML)

model.211 = arima(sunspotsTS_diff,order=c(2,1,1),method='CSS')
coeftest(model.211)

model.211_ML = arima(sunspotsTS_diff,order=c(2,1,1),method='ML')
coeftest(model.211_ML)


model.111 = arima(sunspotsTS_diff,order=c(1,1,1),method='CSS')
coeftest(model.111)

model.111_ML = arima(sunspotsTS_diff,order=c(1,1,1),method='ML')
coeftest(model.111_ML)

# FUNCTION FOR SORTING AIC AND BIC

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

sunspot_AIC = AIC(model.111_ML,model.211_ML,model.012_ML, model.013_ML, model.113_ML, model.312_ML)
sunspot_BIC = BIC(model.111_ML,model.211_ML,model.012_ML, model.013_ML, model.113_ML, model.312_ML)

sort.score(sunspot_AIC, score = "aic")
sort.score(sunspot_BIC, score = "bic")

# ARIMA(0,1,3 - Best Model)

# Overfitting - ARIMA(1,1,3) AND ARIMA (0,1,4)

model.213_OF = arima(sunspotsts,order=c(2,1,3),method='ML')
coeftest(model.213_OF)

model.014_OF = arima(sunspotsts,order=c(0,1,4),method='ML')
coeftest(model.014_OF)


# Residual analysis

par(mar=c(1,1,1,1))

residual.analysis(model = model.211_ML)
residual.analysis(model = model.012_ML)
residual.analysis(model = model.111_ML)
residual.analysis(model = model.113_ML)
residual.analysis(model = model.312_ML)
residual.analysis(model = model.013_ML) # BEST MODEL


par(mar = c(5.1, 4.1, 4.1, 2.1))

# Best Model - ARIMA(0,1,3)


# FORECASTING

fit = Arima(sunspotsts,c(0,1,3))
fitFrc = forecast(fit,h=60)
fitFrc

plot(fitFrc)
