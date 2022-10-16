
# might not need all these libraries but *shrugs*
library(DescTools)
library(lmtest)
library(ggplot2)
library(fBasics)
library(fpp3)
library(forecast)
library(urca)
library(NTS)
library(fUnitRoots)
library(TSA)



# - the functions below are just simple wrappers for various time-series data 
#   tests to help me understand the results in a more layman's perspective
# - (most) tests below are tested against their respective null hypothesis (H0)
#     TRUE: FAIL to reject H0
#     FALSE: reject H0
#   (check each specific function for details)



# ** DISCLAIMER **
# - posted as-is without warranty
# - use at your discretion
# - may contain bugs/errors/changes
# - far from perfect; to be used as growing learning tool/reference
# - collaborative changes/submissions to help improve this script
#     is more than welcome, thank your for your support! =)



# version history:
# v0.01:
#   - made parameter names more explanatory
#   - document inputs and returns for each test
#   - mymcleodlitest - fixed:
#     result <- length(which(ml_m$p.values <= 0.05)) == 0
# v0.02:
#   - emphasized NO/NOT/FAIL in documentation for visibility
#   - rewrote documentation for each test
#   - removed TRUE/FALSE comments for each test:
#       refer to print statement for each test based on (result)
#   - myadftest - updated with drift explanation 
# v0.03:
#   - emphasized *NO*/*NOT*/*FAIL* in documentation for visibility
#   - updated KPSS/ADF information



# ----------------------------------



# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
myttest <- function(data, print_output=TRUE) {
  # Mean 0: T-test
  # t.test() test: check if 95% confidence interval ($conf.int[1] lower, $conf.int[2] upper) 
  #   contains 0 and p-value ($p.value) >= 0.05 (H0)
  ttest_data <- t.test(data)
  result <- ttest_data$conf.int[1] <=0 && ttest_data$conf.int[2] >=0 && ttest_data$p.value >= 0.05
  if (print_output) {
    print(ttest_data)
    if(result)
      print("T-Test: mean is statistically zero, linear trend *REMOVED* -> *FAIL* to reject H0")
    else
      print("T-Test: mean *NOT* statistically zero, linear trend present -> reject H0")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result (i.e. layman terms) - default TRUE
# returns:
#  TRUE: NO skewness
#  FALSE: has skewness
myskewtest <- function(data, print_output=TRUE) {
  # Skew: test for NO skew
  # Skew() test for 0: check if 95% confidence interval ([2] lower, [3] upper) contains zero
  skew_data <- Skew(data, method = 3, conf.level = 0.05, ci.type = "norm", R = 1000)
  result <- skew_data[2] <=0 && skew_data[3] >=0
  if (print_output) {
    print(skew_data)
    if(result)
      print("Skew: *NO* skewness, property conforms to normality and Gaussian PDF")
    else
    {
      skewType <- "LEFT"
      if (skew_data[1] > 0) skewType <- "RIGHT"
      print(sprintf("Skew: has *%s* skewness, property does *NOT* conform to normality and Gaussian PDF", skewType))
    }
  }
  # TRUE: NO skewness
  # FALSE: has skewness
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result (i.e. layman terms) - default TRUE
# returns:
#  TRUE: NO (excess) Kurtosis
#  FALSE: has (excess) Kurtosis
mykurttest <- function(data, print_output=TRUE) {
  # Kurtosis: test for NO (excess) Kurtosis
  # Kurt() test for 0: check if 95% confidence interval ([2] lower, [3] upper) contains zero
  kurt_data <- Kurt(data, method = 3, conf.level = 0.05, ci.type = "norm", R = 1000)
  result <- kurt_data[2] <=0 && kurt_data[3] >=0
  if (print_output) {
    print(kurt_data)
    if(result)
      print("Kurt: *NO* (excess) kurtosis, property conforms to normality and Gaussian PDF")
    else {
      kurtType <- "FLAT thin-tailed"
      if (kurt_data[1] > 0) kurtType <- "TALL thick-tailed"
      print(sprintf("Kurt: has *%s* (excess) kurtosis, property does *NOT* conform to normality and Gaussian PDF", kurtType))
    }
  }
  # TRUE: NO kurtosis
  # FALSE: has kurtosis
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
#   lags: number of lags - default 30
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
myboxljungtest <- function(data,print_output=TRUE,lags=30) {
  # Independence: Box-Ljung test
  # Box.test(type='Ljung') test: check to see if p-value ($p.value) >= 0.05 (H0)
  bl_data <- Box.test(data,lag=lags,type='Ljung') 
  result <- bl_data$p.value >= 0.05
  if (print_output) {
    print(bl_data)
    if(result)
      print(sprintf("Box-Ljung: implies independence over %i lags, *NO* autocorrelation -> *FAIL* to reject H0", lags))
    else
      print(sprintf("Box-Ljung: implies dependency present over %i lags, autocorrelation present -> reject H0", lags))
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
mykpsstest <- function(data, print_output=TRUE) {
  # Stationarity: KPSS Test
  # ur.kpss() test: check if t-stat (@teststat) < critical value for 5% (@cval[2]) (H0)
  kpss_data <- ur.kpss(data,type="tau",lags="short")
  result <- kpss_data@teststat < kpss_data@cval[2]
  if (print_output) {
    print(summary(kpss_data))
    if(result)
      print("KPSS: *NO* unit roots, *NO* linear trend, slope zero, series is trend stationary -> *FAIL* to reject H0")
    else
      print("KPSS: contains unit roots, linear trend present, slope *NOT* zero, series *NOT* trend stationary -> reject H0")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
#   lags: number of lags - default 30
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
myadftest <- function(data,print_output=TRUE,lags=30) {
  # temporarily suppress warning message: 
  # In adfTest(rm, lags = lags, type = "nc"): p-value smaller than printed p-value
  defaultW <- getOption("warn") 
  options(warn = -1) 
  # Stationarity: ADF Test
  # adfTest() test: check to see if p-value (@test$p.value) >= 0.05 (H0)
  adf_data <- adfTest(data,lags=lags,type="nc")
  result <- adf_data@test$p.value >= 0.05
  if (print_output) {
    print(adf_data)
    if(result)
      print(sprintf("ADF: presence of unit roots over %i lags, indicates mean drift, business cycles present, series is *NOT* stationary -> *FAIL* to reject H0", lags))
    else
      print(sprintf("ADF: contains *NO* unit roots over %i lags, indicates *NO* mean drift, business cycles *NOT* present, series is stationary -> reject H0", lags))
  }
  # turn warnings back on
  options(warn = defaultW)
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   model: the Arima() model
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
mymcleodlitest <- function(model, print_output=TRUE) {
  # Constant Variance: McLeod-Li Test
  # McLeod.Li.test() test: check to see if all lags' p-values ($p.values) are greater than 0.05 (H0)
  ml_m <- McLeod.Li.test(model) 
  result <- length(which(ml_m$p.values < 0.05)) == 0
  if (print_output) {
    ml_m
    if(result)
      print("McLeod-Li: constant variance, homoscedastic -> *FAIL* to reject H0")
    else
      print("McLeod-Li: *NON*-constant variance, heteroscedastic -> reject H0")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
mybptest <- function(data, print_output=TRUE) {
  # Constant Variance: Breusch-Pagan Test
  # bptest() test: check to see if p-value ($p.value) >= 0.05 (H0)
  bp_data <- bptest(lm(data ~ seq(1,length(data))))
  result <- bp_data$p.value >= 0.05
  if (print_output) {
    print(bp_data)
    if(result)
      print("Breusch-Pagan: constant variance, homoscedastic -> *FAIL* to reject H0")
    else
      print("Breusch-Pagan: *NON*-constant variance, possible clustering, heteroscedastic -> reject H0")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}
