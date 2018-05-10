#Performs multiple tests on a given data type.
# dataFun - function generating the data, should return the data as specified in "DataTypes.R".
# dataParams - parameters fed to dataFun
# testFuns - vector containing functions for the tests to be applied
# testParams - nested lists containing parameters for testFuns
# nPerm - number of permutations preformed

#Multivariate permutation test
multPermMultTest <- function(dataFun, dataParams, testFuns, testParams, nPerm) {
  obs <- do.call(dataFun, as.list(dataParams)) #Generates data, "odd" function call to be able to pass arguments as vector
  
  nCols <- ncol(obs)
  yObs <- obs[,nCols] #Response variable
  xObs <- obs[,(nCols-1)] #Explanatory variable of interest
  zObs <- obs[,1:(nCols-2)] #Nuisance variables
  
  #Fit of original observation
  obsFit <- lm(yObs ~ xObs + zObs, model=FALSE)
  
  #Fit of reduced model (Freedman-Lane)
  redYZFit <- lm(yObs ~ zObs, model=FALSE)
  
  #Permute reduced model and fit full model using permutations
  permYZRes <- replicate(nPerm, sample(redYZFit$residuals))
  permY <- permYZRes + zObs*redYZFit$coefficients["zObs"]
  permFit <- lm(permY ~ xObs + zObs)
  
  #Summaries of fits for extracting value of test statistics
  obsSum <- summary(obsFit)
  permSum <- summary(permFit)
  
  #Convert strings to variable names
  testParams <- lapply(testParams, function(x) lapply(x, as.name))
  
  #Perform the tests
  res <- sapply(1:length(testFuns), function(x) do.call(testFuns[[x]], testParams[[x]]))
  return(c(res, kurtosis(obsFit$residuals), skewness(obsFit$residuals)))
}

#Function used to flip signs randomly for the sign-flip permutation test
flipfun <- function(len){
  temp <- rbinom(len, 1, prob=0.5)
  temp[temp==0] <- -1
  return(temp)
}

#Multivariate sign flip and permutation test
multSignAndPermMultTest <- function(dataFun, dataParams, testFuns, testParams, nPerm) {
  obs <- do.call(dataFun, as.list(dataParams)) #Generates data, "odd" function call to be able to pass arguments as vector
  
  nCols <- ncol(obs)
  yObs <- obs[,nCols] #Response variable
  xObs <- obs[,(nCols-1)] #Explanatory variable of interest
  zObs <- obs[,1:(nCols-2)] #Nuisance variables
  
  #Fit of original observation
  obsFit <- lm(yObs ~ xObs + zObs, model=FALSE)
  
  #Fit of reduced model (Freedman-Lane)
  redYZFit <- lm(yObs ~ zObs, model=FALSE)
  
  #Permute residuals of reduced model and fit full model using permuted residuals
  permYZRes <- replicate(nPerm, sample(redYZFit$residuals))
  permY <- permYZRes + zObs*redYZFit$coefficients["zObs"]
  permFit <- lm(permY ~ xObs + zObs)
  
  #Sign flip residuals of reduced model and fit full model using new residuals
  signYZRes <- redYZFit$residuals*replicate(nPerm, sample(flipfun(length(redYZFit$residuals))))
  signY <- signYZRes + zObs*redYZFit$coefficients["zObs"]
  signFit <- lm(signY ~ xObs + zObs)
  
  #Summaries of fits for extracting value of test statistics
  obsSum <- summary(obsFit)
  permSum <- summary(permFit)
  signSum <- summary(signFit)
  
  #For sign test, has to be changed manually
  signTestParams <- list(c("xObs", "zObs", "redYZFit", "signYZRes"), c("obsSum", "signSum"))
  signTestParams <- lapply(signTestParams, function(x) lapply(x, as.name))
  signTestFuns <- c("multRSquaredTest", "tTest")
  
  #Convert strings to variable names
  testParams <- lapply(testParams, function(x) lapply(x, as.name))
  
  #Perform the tests
  permRes <- sapply(1:length(testFuns), function(x) do.call(testFuns[[x]], testParams[[x]]))
  signRes <- sapply(1:length(signTestFuns), function(x) do.call(signTestFuns[[x]], signTestParams[[x]]))
  return(c(signRes, permRes, kurtosis(obsFit$residuals), skewness(obsFit$residuals)))
}


#Simple permutation test
simplePermMultTest <- function(dataFun, dataParams, testFuns, testParams, nPerm) {
  obs <- do.call(dataFun, as.list(dataParams)) #Generates data, "odd" function call to be able to pass arguments as vector

  #Fit of original observation
  obsFit <- lm(obs[,2]~obs[,1], model=FALSE)
  
  permutations <- replicate(nPerm, sample(obs[,2]))  #Permuted y-values
  permFit <- lm(permutations~obs[,1]) #Fits of permuted dependent variables
  
  #Summaries of fits for extracting value of test statistics
  permSum <- summary(permFit)
  obsSum <- summary(obsFit)
  
  #Convert strings to variable names
  testParams <- lapply(testParams, function(x) lapply(x, as.name))
  
  #Perform the tests
  res <- sapply(1:length(testFuns), function(x) do.call(testFuns[[x]], testParams[[x]]))
  return(c(res, kurtosis(obsFit$residuals), skewness(obsFit$residuals)))
}


#Simple sign flip and permutation test
simpleSignAndPermMultTest <- function(dataFun, dataParams, testFuns, testParams, nPerm) {
  obs <- do.call(dataFun, as.list(dataParams)) #Generates data, "odd" function call to be able to pass arguments as vector
  
  #Fit of original observation
  obsFit <- lm(obs[,2]~obs[,1], model=FALSE)
  
  permutations <- replicate(nPerm, sample(obs[,2]))  #Permuted y-values
  permFit <- lm(permutations~obs[,1]) #Fits of permuted dependent variables
  
  signfliped <- obs[,2]*replicate(nPerm, flipfun(length(obs[,2])))  #Sign fliped y-values
  signFit <- lm(signfliped~obs[,1]) #Fits of permuted dependent variables
  
  #Summaries of fits for extracting value of test statistics
  permSum <- summary(permFit)
  signSum <- summary(signFit)
  obsSum <- summary(obsFit)
  
  #For sign test, has to be changed manually
  signTestParams <- list(c("obsSum", "signSum"), c("obsSum", "signSum"))
  signTestParams <- lapply(signTestParams, function(x) lapply(x, as.name))
  signTestFuns <- c("simpleRSquaredTest", "tTest")
  
  #Convert strings to variable names
  testParams <- lapply(testParams, function(x) lapply(x, as.name))
  
  #Perform the tests
  permRes <- sapply(1:length(testFuns), function(x) do.call(testFuns[[x]], testParams[[x]]))
  signRes <- sapply(1:length(signTestFuns), function(x) do.call(signTestFuns[[x]], signTestParams[[x]]))
  return(c(signRes, permRes, kurtosis(obsFit$residuals), skewness(obsFit$residuals)))
}