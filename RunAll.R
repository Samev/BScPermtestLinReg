#Script for running all simulations presented in the report. Notice that a seed is set at row 40, this was the seed used for the final simulation. 
#To replicate the results of the report the complete script should be run in one go.
#All simulations took roughly 10h to run on 4 cores on a stock clocked intel i5-6600K.

#Loads different types of data and testFuns from working directory
source("DataTypes.R")
source("TestTypes.R")
source("MultipleTests.R")

#For creating and naming directories
dirFunction <- function(dataFun){
  dirBase <- paste("./", dataFun, "_", Sys.Date(), sep="")
  dirName <- dirBase
  dirCount <- 0
  while(dir.exists(dirName)){
    dirCount <- dirCount + 1
    dirName <- paste(dirBase, dirCount, sep="_")
  }
  dir.create(dirName)
  return(dirName)
}

########
#Run from here after loading functions from above.
library('moments')
library('parallel')

reps <- 500 #Number of replicates per parameter setup
perms <- 999 #Number of permutations
nVec <- 10^(1:3) #Sample sizes
betaVec <- c(0, 0.1, 0.5, 1) #Coefficients for variable of interest
testFuns <- c("simpleRSquaredTest", "tTest", "normTest") #Tests
testPars<- list(c("obsSum", "permSum"), c("obsSum", "permSum"), "obsSum") #Parameters of the tests
signDFName <- c("multRSquaredTestSign", "tTestSign") #Only used in cases where sign permutations occurs
rho <- 0.5 #Correlation, only used in multivariate case

no_cores <- detectCores() #Number of cores for parallel computation. detectCores()-1 is nice if u want to use computer during simulation
cl <- makeCluster(no_cores)
clusterExport(cl, ls()) #Exports current workspace to cluster
clusterEvalQ(cl, library('moments'))
clusterEvalQ(cl, library('MASS'))
clusterSetRNGStream(cl, 8491) #Set seed for the simulations

####

dataFun <- "simpleHeavyTail"
paramVec <- c(3:6, 8, 10, 25, 100) #Degrees of freedom of t-distribution used to generate heavy tails
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- c(n, beta, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){simplePermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

####
dataFun <- "simpleExpSkew_XUnif"
paramVec <- c(1, 1.5, 2, 2.5, 3, 5, 7, 10) #Skewness factor
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- c(n, beta, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){simplePermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))


####

betaVec <- c(-1, -0.5, -0.1, 0, 0.1, 0.5, 1)
clusterExport(cl, "betaVec")
dataFun <- "simpleExpSkew_XExp"
paramVec <- c(1, 1.5, 2, 2.5, 3, 5, 7, 10) #Skewness factor
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- c(n, beta, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){simplePermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

betaVec <- c(0, 0.1, 0.5, 1)
clusterExport(cl, "betaVec")

####

dataFun <- "simpleOutliers"
paramVec <- c(1:10/100) #Percentage of outliers
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- c(n, beta, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){simplePermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

####


dataFun <- "simpleNonConstVar"
paramVec <- c(1:10) #Coefficient determinening how quickly variance changes with x
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- c(n, beta, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){simpleSignAndPermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(signDFName, testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))


####

#Changing parameters for use in multivariate case
testFuns <- c("multRSquaredTest", "tTest", "normTest") #Tests
testPars<- list(c("xObs", "zObs", "redYZFit", "permYZRes"), c("obsSum", "permSum"), "obsSum") #Parameters of the tests
clusterExport(cl, c("testFuns", "testPars"))

####


dataFun <- "multHeavyTail"
paramVec <- c(3:6, 8, 10, 25, 100) #Degrees of freedom of t-distribution used to generate heavy tails
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- list(n, c(0.5, beta), rho, 1, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){multPermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

####


dataFun <- "multExpSkew"
paramVec <- c(1, 1.5, 2, 2.5, 3, 5, 7, 10) #Skewness factor
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- list(n, c(0.5, beta), rho, 1, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){multPermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

####


dataFun <- "multOutliers"
paramVec <- c(1:10/100) #Percentage of outliers
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- list(n, c(0.5, beta), rho, 1, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){multPermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

####


dataFun <- "multNonConstVar"
paramVec <- c(1:10) #Coefficient determinening how quickly variance changes with x
dirName <- dirFunction(dataFun) #Setup directory for storing results
clusterExport(cl, c("dataFun", "paramVec")) #Exports current workspace to cluster

runTime <- system.time(
  for(n in nVec){
    for(beta in betaVec){
      for(v in paramVec){
        dataParams <- list(n, c(0.5, beta), rho, 1, v)
        clusterExport(cl, "dataParams")
        pRes <- parSapply(cl, 1:reps, function(x){multSignAndPermMultTest(dataFun, dataParams, testFuns, testPars, perms)})
        filename <- paste(dirName, "/", dataFun, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n, "_reps", reps, "_perms", perms, "_", paste(testFuns, collapse="-"), ".rds", sep="")
        pResSave <- data.frame(pRes, row.names=c(signDFName, testFuns,"kurtosis","skewness"))
        saveRDS(pResSave, file=filename)
      }
    }
  }
)

save(list=c("reps", "perms", "nVec", "betaVec", "paramVec", "testFuns", "testPars", "dataFun", "no_cores", "runTime"), file=paste(dirName, "/parameters.RData", sep=""))

stopCluster(cl)