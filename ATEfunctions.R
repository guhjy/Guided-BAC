setwd(paste("/Users/josephantonelli/Documents/joey/Research/Air_pollution/",
      "control_confounding/multiple_imputation/code/Github", sep=''))

library("MASS")
library("mvtnorm")
library("truncnorm")

source("ContinuousBinary.R")
source("BinaryBinary.R")
