ContinuousBinary = function(Cobs, Cmis, covType, X, Y, savePosterior = FALSE,
                                nScans=5000, burn=1000, thin=4, omega=1000000) {
  
  ##################################
  ## Function that calculates     ##
  ## ATE when treatment is        ##
  ## binary and outcome is        ##
  ## continuous                   ##
  ##################################
  
  N = nrow(Cobs)
  n = sum(!is.na(Cmis[,1]))
  
  M = ncol(Cobs)
  P = ncol(Cobs) + ncol(Cmis)
  
  ## indicator of whether subject is in main or validation data (1 represents validation)
  StudyInd = c(rep(1, n), rep(0, N-n))
  
  ## Combine fully observed and missing data
  C = cbind(Cobs, Cmis)
  
  ## which MCMC scans to keep
  keep = (burn+1) + 0:((nScans - burn)/thin - 1)*thin
  
  ## MCMC storage
  alphaXpost = alphaYpost = array(NA, dim=c(nScans, P))
  sigmaXpost = sigmaYpost = array(NA, dim=c(nScans, 1))
  sigmaCpost = array(NA, dim=c(nScans, (P-M)))
  thetaXpost = array(NA, dim=c(nScans, P+M+1))
  thetaYpost = array(NA, dim=c(nScans, P+M+3))
  thetaCpost = array(NA, dim=c(nScans, (P-M), P))
  alphaCpost = array(NA, dim=c(nScans, (P-M), P))
  misCpost = latentCpost = matrix(NA, nrow=P-M, ncol=N)
  
  ## set starting values
  
  alphaXpost[1,] = sample(c(0,1), P, replace=TRUE)
  alphaYpost[1,] = sample(c(0,1), P, replace=TRUE)
  sigmaXpost[1] = 1
  sigmaYpost[1] = 1
  sigmaCpost[1,] = rep(.1, (P-M))
  thetaXpost[1,] = rep(0, P+M+1)
  thetaYpost[1,] = rep(0, P+M+3)
  thetaCpost[1,,] = rep(0, P*(P-M))
  alphaCpost[1,,] = rep(0, P*(P-M))
  
  for (p in 1 : (P-M)) {
    misCpost[p,] = Cmis[,p]
    latentCpost[p,] = Cmis[,p]
    misCpost[p,(n+1):N] = mean(Cmis[,p], na.rm=TRUE)
    latentCpost[p,(n+1):N] = mean(Cmis[,p], na.rm=TRUE)
  }
  
  acc.alphaX = acc.alphaY = 0
  
  ## prior variance of regression parameters
  k=1000000
  
  ## prior hyperparameters for inverse gammas
  a=0.001
  b=0.001
  
  for (i in 2 : nScans) {
    
    if (i %% 100 == 0) print(i)
    
    ## update full data matrix conditional on missing values imputed
    Cupdate = cbind(Cobs, t(misCpost), (Cobs*StudyInd))
    CupdateLatent = cbind(Cobs, t(latentCpost), (Cobs*StudyInd))
    
    ## draw from alphaX using BIC approximation
    
    ## latent continuous gibbs step
    Zx = rep(NA, N)
    
    whichColX = c(which(alphaXpost[i-1,]==1), which(alphaXpost[i-1,1:M]==1) + P)
    meanZx = cbind(rep(1,N), Cupdate[,whichColX]) %*% 
      thetaXpost[i-1,c(1,whichColX+1)]
    
    Zx[X==1] = rtruncnorm(sum(X==1), a=0, mean = meanZx[X==1], sd=1)
    Zx[X==0] = rtruncnorm(sum(X==0), b=0, mean = meanZx[X==0], sd=1)
    
    alphaXprop = alphaXpost[i-1, ]
    wX = sample(1:P, 1)
    alphaXprop[wX] = 1 - alphaXpost[i-1, wX]
    
    wXprop = c(which(alphaXprop == 1), which(alphaXprop[1:M] == 1)+P)
    wXpost = c(which(alphaXpost[i-1,] == 1), which(alphaXpost[i-1,1:M] == 1)+P)
    
    temp.omega = omega^(alphaYpost[i-1,wX]==0)
    temp.odds = temp.omega
    
    if (alphaXprop[wX]==1) temp.odds = 1/temp.omega
    
    if (length(wXprop) >= 1) {
      fitXprop = glm(X ~ Cupdate[,wXprop], family=binomial(link="probit"))
    } else fitXprop = glm(X ~ 1, family=binomial(link="probit"))
    
    if (length(wXpost) >= 1) {
      fitXpost = glm(X ~ Cupdate[,wXpost], family=binomial(link="probit"))
    } else fitXpost = glm(X ~ 1, family=binomial(link="probit"))
    
    BICpropX = (-1)*BIC(fitXprop)/2
    BICpostX = (-1)*BIC(fitXpost)/2
    logBFX = BICpropX - BICpostX
    BFX = exp(logBFX)*temp.odds
    
    if (runif(1) > BFX) {
      alphaXpost[i,] = alphaXpost[i-1,]
    } else {
      alphaXpost[i,] = alphaXprop
      acc.alphaX = acc.alphaX + 1
    }
    
    ## draw from alphaY using BIC approximation
    
    tildeY = Y - (thetaYpost[i-1, 2]*X) - thetaYpost[i-1, 3]*X*StudyInd
    
    alphaYprop = alphaYpost[i-1, ]
    wY = sample(1:P, 1)
    alphaYprop[wY] = 1 - alphaYpost[i-1, wY]
    
    wYprop = c(which(alphaYprop == 1), which(alphaYprop[1:M] == 1)+P)
    wYpost = c(which(alphaYpost[i-1,] == 1), which(alphaYpost[i-1,1:M] == 1)+P)
    
    if (length(wYprop) >= 1) {
      fitYprop = lm(tildeY ~ Cupdate[,wYprop])
    } else fitYprop = lm(tildeY ~ 1)
    
    if (length(wYpost) >= 1) {
      fitYpost = lm(tildeY ~ Cupdate[,wYpost])
    } else fitYpost = lm(tildeY ~ 1)
    
    temp.omega = omega^(alphaXpost[i,wY]==1)
    
    temp.odds = temp.omega
    if (alphaYprop[wY]==0) temp.odds = 1/temp.omega
    
    ## adding in conditional BAC priors
    
    BICpropY = (-1)*BIC(fitYprop)/2
    BICpostY = (-1)*BIC(fitYpost)/2
    logBFY = BICpropY - BICpostY
    BFY = exp(logBFY)*temp.odds
    
    if (runif(1) > BFY) {
      alphaYpost[i,] = alphaYpost[i-1,]
    } else {
      alphaYpost[i,] = alphaYprop
      acc.alphaY = acc.alphaY + 1
    }
    
    ## draw from full conditional of thetaX
    
    wXtheta = c(which(alphaXpost[i,] == 1), which(alphaXpost[i,1:M] == 1)+P)
    Wx = cbind(rep(1, N), Cupdate[, wXtheta])
    
    muX = solve((t(Wx)%*%Wx) + (diag(dim(Wx)[2])*sigmaXpost[i-1]/k)) %*% t(Wx) %*% Zx
    covX = solve(((t(Wx)%*%Wx)/sigmaXpost[i-1]) + (diag(dim(Wx)[2])/k))
    
    thetaXpost[i,] = rep(0, length(thetaXpost[i,]))
    thetaXpost[i,c(1, wXtheta+1)] = rmvnorm(1, muX, covX)
    
    ## draw from full conditional of sigmaX
    
    sigmaXpost[i] = 1
    
    ## draw from full conditional of thetaY
    
    wYtheta = c(which(alphaYpost[i,] == 1), which(alphaYpost[i,1:M] == 1)+P)
    Wy = cbind(rep(1, N), X, X*StudyInd, Cupdate[, wYtheta])
    
    muY = solve((t(Wy)%*%Wy) + (diag(dim(Wy)[2])*sigmaYpost[i-1]/k)) %*% t(Wy) %*% Y
    covY = solve(((t(Wy)%*%Wy)/sigmaYpost[i-1]) + (diag(dim(Wy)[2])/k))
    
    thetaYpost[i,] = rep(0, length(thetaYpost[i,]))
    thetaYpost[i, c(1, 2, 3, wYtheta+3)] = rmvnorm(1, muY, covY)
    
    ## draw from full conditional of sigmaY
    
    temp.thetaY = thetaYpost[i, c(1, 2, 3, wYtheta+3)]
    aY = a + N/2
    bY = as.numeric(b + (t((Y - Wy %*% temp.thetaY)) %*% (Y - Wy %*% temp.thetaY))/2)
    
    sigmaYpost[i] = 1/rgamma(1, aY, bY)
    
    for (p in (M+1) : P) {
      
      if (covType[p]=="Binary") {
        ## variable inclusions
        alphaCprop = alphaCpost[i-1, p-M, 1:(p-1)]
        wC = sample(1:(p-1), 1)
        alphaCprop[wC] = 1 - alphaCpost[i-1, p-M, wC]
        
        wCprop = which(alphaCprop == 1)
        wCpost = which(alphaCpost[i-1, p-M, 1:(p-1)] == 1)
        
        if (length(wCprop) >= 1) {
          fitCprop = lm(CupdateLatent[,p] ~ CupdateLatent[,wCprop])
        } else fitCprop = lm(CupdateLatent[,p] ~ 1)
        
        if (length(wCpost) >= 1) {
          fitCpost = lm(CupdateLatent[,p] ~ CupdateLatent[,wCpost])
        } else fitCpost = lm(CupdateLatent[,p] ~ 1)
        
        BICpropC = (-1)*BIC(fitCprop)/2
        BICpostC = (-1)*BIC(fitCpost)/2
        logBFC = BICpropC - BICpostC
        BFC = exp(logBFC)
        
        if (runif(1) > BFC) {
          alphaCpost[i,p-M,1:(p-1)] = alphaCpost[i-1,p-M,1:(p-1)]
        } else {
          alphaCpost[i,p-M,1:(p-1)] = alphaCprop
        }
        
        wCtheta = which(alphaCpost[i,p-M,1:(p-1)] == 1)
        
        Wp = cbind(rep(1, N), CupdateLatent[,wCtheta])
        
        muP = solve((t(Wp) %*% Wp) + (diag(dim(Wp)[2])*sigmaCpost[i-1,p-M]/k)) %*% t(Wp) %*% latentCpost[p-M, ]
        covP = solve(((t(Wp) %*% Wp)/sigmaCpost[i-1,p-M]) + (diag(dim(Wp)[2])/k))
        
        thetaCpost[i,p-M,(1:p)] = rep(0, p)
        thetaCpost[i,p-M,c(1, wCtheta+1)] = rmvnorm(1, muP, covP)
        
        sigmaCpost[i,p-M] = 1
        
        ## Now update latent continuous variable distribution
        
        ## get component of distribution from variable p
        
        muPcomp = rep(thetaCpost[i,p-M,1], dim(Cupdate)[1])
        for (num.p in 1:(p-1)) {
          muPcomp = muPcomp + thetaCpost[i,p-M,num.p+1]*CupdateLatent[,num.p]
        }
        muPcomp = muPcomp/sigmaCpost[i,p-M]
        varPcomp = 1/sigmaCpost[i,p-M]
        
        ## get component from other missing variables
        
        muOtherPcomp = 0
        varOtherPcomp = 0
        
        if (p == P) {
          muOtherPcomp = 0
          varOtherPcomp = 0
        } else {
          for (j in (p+1) : P) {
            varOtherPcomp = varOtherPcomp + ((thetaCpost[i-1,j-M,p+1]^2) / sigmaCpost[i-1,j-M])
            tempWp = cbind(rep(1,N), CupdateLatent)
            tempOtherPcomp = latentCpost[j-M,] - tempWp[,(1:j)[-(p+1)]]%*%thetaCpost[i-1,j-M,(1:j)[-(p+1)]]
            muOtherPcomp = muOtherPcomp + (thetaCpost[i-1,j-M,p+1] / sigmaCpost[i-1,j-M])*tempOtherPcomp
          }
        }
        
        covM = 1/(varPcomp + varOtherPcomp)
        muM = covM*(muPcomp + muOtherPcomp)
        
        ## Now we need components from Y and X models for the probability of which truncated normal to sample from
        
        ## create new alpha vector so that p is always zero
        alphaYzero = c(alphaYpost[i,], alphaYpost[i,1:M])
        alphaYzero[p] = 0
        alphaXzero = c(alphaXpost[i,], alphaXpost[i,1:M])
        alphaXzero[p] = 0
        
        if (alphaYpost[i,p] == 1) {
          muYcomp = Y - thetaYpost[i,2]*X - thetaYpost[i,3]*X*StudyInd - rep(thetaYpost[i,1], length(X))
          for (num.p in 1:(P+M)) {
            if(alphaYzero[num.p] == 1) {
              muYcomp = muYcomp - alphaYzero[num.p]*thetaYpost[i,num.p+3]*Cupdate[,num.p]
            }
          }
        } else {
          muYcomp = NA
        }
        
        if (alphaXpost[i,p] == 1) {
          muXcomp = Zx - rep(thetaXpost[i,1], length(X))
          for (num.p in 1:(P+M)) {
            if(alphaXzero[num.p] == 1) {
              muXcomp = muXcomp - alphaXzero[num.p]*thetaXpost[i,num.p+1]*Cupdate[,num.p]
            }
          }
        } else {
          muXcomp = NA
        }
        
        ## TODO do I need this or can I have one all encompassing scenario
        ## I think it can be put into one if I turn my thetaY and thetaX NA's into zeros
        
        if (is.na(muYcomp)[1]==TRUE & is.na(muXcomp)[1]==TRUE) {
          PI = pnorm(muM/sqrt(covM))
        } else if (is.na(muYcomp)[1]==TRUE & is.na(muXcomp)[1]==FALSE) {
          PI = ((pnorm(muM/sqrt(covM)))*dnorm(muXcomp - thetaXpost[i,p+1], sd=sqrt(sigmaXpost[i]))) /
            (((pnorm(muM/sqrt(covM)))*dnorm(muXcomp - thetaXpost[i,p+1], sd=sqrt(sigmaXpost[i]))) +
               ((1-pnorm(muM/sqrt(covM)))*dnorm(muXcomp, sd=sqrt(sigmaXpost[i]))))
        } else if (is.na(muYcomp)[1]==FALSE & is.na(muXcomp)[1]==TRUE) {
          PI = ((pnorm(muM/sqrt(covM)))*dnorm(muYcomp - thetaYpost[i,p+3], sd=sqrt(sigmaYpost[i]))) /
            (((pnorm(muM/sqrt(covM)))*dnorm(muYcomp - thetaYpost[i,p+3], sd=sqrt(sigmaYpost[i]))) +
               ((1 - pnorm(muM/sqrt(covM)))*dnorm(muYcomp, sd=sqrt(sigmaYpost[i]))))
        } else {
          PI = ((pnorm(muM/sqrt(covM)))*dnorm(muYcomp - thetaYpost[i,p+3], sd=sqrt(sigmaYpost[i]))*
                  dnorm(muXcomp - thetaXpost[i,p+1], sd=sqrt(sigmaXpost[i]))) /
            (((pnorm(muM/sqrt(covM)))*dnorm(muYcomp - thetaYpost[i,p+3], sd=sqrt(sigmaYpost[i]))*
                dnorm(muXcomp - thetaXpost[i,p+1], sd=sqrt(sigmaXpost[i]))) +
               ((1 - pnorm(muM/sqrt(covM)))*dnorm(muYcomp, sd=sqrt(sigmaYpost[i]))*
                  dnorm(muXcomp, sd=sqrt(sigmaXpost[i]))))
        }
        
        misCpost[p-M,] = rbinom(N, 1, p=PI)
        misCpost[p-M, 1:n] = Cmis[1:n, p-M]
        
        Cupdate[,p] = misCpost[p-M,]
        
        ## first update latent variable for all subjects
        latentCpost[p-M, ][Cupdate[,p]==1] = rtruncnorm(sum(Cupdate[,p]), a=0, 
                                                        mean=muM[Cupdate[,p]==1], 
                                                        sd=sqrt(covM))
        
        latentCpost[p-M, ][Cupdate[,p]==0] = rtruncnorm(sum(Cupdate[,p]==0), b=0, 
                                                        mean=muM[Cupdate[,p]==0], 
                                                        sd=sqrt(covM))
        
        
        CupdateLatent[,p] = latentCpost[p-M,]
        
        
        if (i %% 1000 == 0) print(table(misCpost[p-M,], Cmis[,p-M])[1:4] - 
                                    table(misCpost[p-M,1:n], Cmis[1:n,p-M])[1:4])
      } else {
        
        ## variable inclusions
        alphaCprop = alphaCpost[i-1, p-M, 1:(p-1)]
        wC = sample(1:(p-1), 1)
        alphaCprop[wC] = 1 - alphaCpost[i-1, p-M, wC]
        
        wCprop = which(alphaCprop == 1)
        wCpost = which(alphaCpost[i-1, p-M, 1:(p-1)] == 1)
        
        if (length(wCprop) >= 1) {
          fitCprop = lm(CupdateLatent[,p] ~ CupdateLatent[,wCprop])
        } else fitCprop = lm(CupdateLatent[,p] ~ 1)
        
        if (length(wCpost) >= 1) {
          fitCpost = lm(CupdateLatent[,p] ~ CupdateLatent[,wCpost])
        } else fitCpost = lm(CupdateLatent[,p] ~ 1)
        
        BICpropC = (-1)*BIC(fitCprop)/2
        BICpostC = (-1)*BIC(fitCpost)/2
        logBFC = BICpropC - BICpostC
        BFC = exp(logBFC)
        
        if (runif(1) > BFC) {
          alphaCpost[i,p-M,1:(p-1)] = alphaCpost[i-1,p-M,1:(p-1)]
        } else {
          alphaCpost[i,p-M,1:(p-1)] = alphaCprop
        }
        
        wCtheta = which(alphaCpost[i,p-M,1:(p-1)] == 1)
        
        Wp = cbind(rep(1, N), CupdateLatent[,wCtheta])
        
        muP = solve((t(Wp) %*% Wp) + (diag(dim(Wp)[2])*sigmaCpost[i-1,p-M]/k)) %*% t(Wp) %*% latentCpost[p-M, ]
        covP = solve(((t(Wp) %*% Wp)/sigmaCpost[i-1,p-M]) + (diag(dim(Wp)[2])/k))
        
        thetaCpost[i,p-M,(1:p)] = rep(0, p)
        thetaCpost[i,p-M,c(1, wCtheta+1)] = rmvnorm(1, muP, covP)
        
        ## Need to update confounder p's variance component
        
        aP = a + N/2
        bP = as.numeric(b + (t((latentCpost[p-M,] - Wp %*% thetaCpost[i,p-M,c(1,wCtheta+1)])) %*% 
                               (latentCpost[p-M,] - Wp %*% thetaCpost[i,p-M,c(1,wCtheta+1)]))/2)
        
        sigmaCpost[i,p-M] = 1/rgamma(1, aP, bP)
        
        ## Need to update missing data for confounder p
        
        ## create new alpha vector so that p is always zero
        alphaYzero = c(alphaYpost[i,], alphaYpost[i,1:M])
        alphaYzero[p] = 0
        alphaXzero = c(alphaXpost[i,], alphaXpost[i,1:M])
        alphaXzero[p] = 0
        
        if (alphaYpost[i,p] == 1) {
          varYcomp = alphaYpost[i,p]*thetaYpost[i,p+3]^2/sigmaYpost[i]
          muYcomp = Y - thetaYpost[i,2]*X - thetaYpost[i,3]*X*StudyInd - rep(thetaYpost[i,1], length(X))
          for (num.p in 1:(P+M)) {
            if(alphaYzero[num.p] == 1) {
              muYcomp = muYcomp - alphaYzero[num.p]*thetaYpost[i,num.p+3]*Cupdate[,num.p]
            }
          }
          muYcomp = (alphaYpost[i,p]*thetaYpost[i,p+3]/sigmaYpost[i])*muYcomp
        } else {
          muYcomp = 0
          varYcomp = 0
        }
        
        if (alphaXpost[i,p] == 1) {
          varXcomp = alphaXpost[i,p]*thetaXpost[i,p+1]^2/sigmaXpost[i]
          muXcomp = Zx - rep(thetaXpost[i,1], length(X))
          for (num.p in 1:(P+M)) {
            if(alphaXzero[num.p] == 1) {
              muXcomp = muXcomp - alphaXzero[num.p]*thetaXpost[i,num.p+1]*Cupdate[,num.p]
            }
          }
          muXcomp = (alphaXpost[i,p]*thetaXpost[i,p+1]/sigmaXpost[i])*muXcomp
        } else {
          varXcomp = 0
          muXcomp = 0
        }
        
        muPcomp = rep(thetaCpost[i,p-M,1], dim(Cupdate)[1])
        for (num.p in 1:(p-1)) {
          muPcomp = muPcomp + thetaCpost[i,p-M,num.p+1]*CupdateLatent[,num.p]
        }
        muPcomp = muPcomp/sigmaCpost[i,p-M]
        varPcomp = 1/sigmaCpost[i,p-M]
        
        muOtherPcomp = 0
        varOtherPcomp = 0
        
        if (p == P) {
          muOtherPcomp = 0
          varOtherPcomp = 0
        } else {
          for (j in (p+1) : P) {
            varOtherPcomp = varOtherPcomp + ((thetaCpost[i-1,j-M,p+1]^2) / sigmaCpost[i-1,j-M])
            tempWp = cbind(rep(1,N), CupdateLatent)
            tempOtherPcomp = latentCpost[j-M,] - tempWp[,(1:j)[-(p+1)]]%*%thetaCpost[i-1,j-M,(1:j)[-(p+1)]]
            muOtherPcomp = muOtherPcomp + (thetaCpost[i-1,j-M,p+1] / sigmaCpost[i-1,j-M])*tempOtherPcomp
          }
        }
        
        covM = 1/(varPcomp + varXcomp + varYcomp + varOtherPcomp)
        muM = covM*(muYcomp + muXcomp + muPcomp + muOtherPcomp)
        
        misCpost[p-M,] = rnorm(N, mean=muM, sd=sqrt(covM))
        misCpost[p-M, 1:n] = Cmis[1:n, p-M]
        latentCpost[p-M,] = misCpost[p-M,]
        
        Cupdate[,p] = misCpost[p-M,]
        CupdateLatent[,p] = misCpost[p-M,]
      }
    }
  }
  
  l = list(InclusionProbX = apply(alphaXpost[keep,], 2, mean, na.rm=TRUE),
           InclusionProbY = apply(alphaYpost[keep,], 2, mean, na.rm=TRUE),
           ATE = mean(thetaYpost[keep,2], na.rm=TRUE),
           ATE_CI = quantile(thetaYpost[keep,2], c(.025, .975), na.rm=TRUE))
  
  if (savePosterior == FALSE) {
    l = list(InclusionProbX = apply(alphaXpost[keep,], 2, mean, na.rm=TRUE),
             InclusionProbY = apply(alphaYpost[keep,], 2, mean, na.rm=TRUE),
             ATE = mean(thetaYpost[keep,2], na.rm=TRUE),
             ATE_CI = quantile(thetaYpost[keep,2], c(.025, .975), na.rm=TRUE))
  } else {
    l = list(InclusionProbX = alphaXpost[keep,],
             InclusionProbY = alphaYpost[keep,],
             ATE = thetaYpost[keep,2])
  }
  
  return(l)
}