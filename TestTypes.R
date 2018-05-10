# Different tests for H0: beta = 0, HA: beta != 0, where beta is coefficient of the explanatory variable of interest
# obsSum - summary of fitted linear model of original observation
# permSum - summary of fitted linear model of permuted observations
# permYZRes - permuted residuals of the reduced model (Freedman-Lane)
# redYZFit - fitted linear model of the reduced model (Freedman-Lane)
# xObs - observations of explanatory variable of interest
# zObs - observations of nuicance variables

# Simple R-squared test
# returns p-value
simpleRSquaredTest <- function(obsSum, permSum){
  rObs <- obsSum$r.squared
  rVals <- c(rObs, sapply(1:length(permSum), function(i) permSum[[i]]$r.squared))
  return(sum(rObs <= rVals)/length(rVals))
}

# Multivariate R-squared
# returns p-value
multRSquaredTest <- function(xObs, zObs, redYZFit, permYZRes){
  redXZFit <- lm(xObs ~ zObs)
  rObs <- sum(redYZFit$residuals*redXZFit$residuals)^2/(sum(redYZFit$residuals^2)*sum(redXZFit$residuals^2))
  rVals <- c(rObs, apply(permYZRes, 2, function(i) sum(i*redXZFit$residuals)^2/(sum(i^2)*sum(redXZFit$residuals^2))))
  return(sum(rObs <= rVals)/length(rVals))
}

# Two sided t-test of permutationed dependent variables (including original observation)
# Uses 2min(low,high) to deal with case of non-symmetric reference distributions and absolute p-value for comparison
# Values for independent variable analysed is on second row of summary(...)$coefficients
# returns p-value
tTest <- function(obsSum, permSum){
  tObs <- obsSum$coefficients[2,3]
  tVals <- sapply(1:length(permSum), function(i) permSum[[i]]$coefficients[2,3])

  low <- sum(tObs <= tVals)
  high <- sum(tObs >= tVals)
  return(c((2*min(low,high)+1)/(length(tVals)+1), (sum(abs(tObs) <= abs(tVals))+1)/(length(tVals)+1)))
}

# t-test
# returns T-statistic
normTest <- function(obsSum){
  return(obsSum$coefficients[2,3])
}