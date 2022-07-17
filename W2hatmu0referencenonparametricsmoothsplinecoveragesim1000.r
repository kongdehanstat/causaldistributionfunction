### This code is the data adaptive cross fitting estimator for Scenario 2, 
### where both the propensity score model and the outcome regression model are nonlinear functions of X

rm(list=ls())
library(MASS)
library(GoFKernel)
library(dplyr)
library(gss)

# Sample size
n<-1000

# Number of MC
nrep<-1100

# Dimension of confounders
p<-1

# Number of samples for empirical CDF for each subject
m<-1000

# correlation between X
#rho<-0.5

#Covariance of X
#Sigma <- rho^t(sapply(1:p, function(i, j) abs(i-j), 1:p))


truegamma0<-1
truegamma<-1

truetreatmentfun<-function(u){
  truetreatmentfun<-sin(pi*u)/8
}

truebetafun<-function(u){
  truebetafun<-sin(pi*u)/8
}

gfun<-function(u){
  gfun<-(1+exp(-u))^(-1)
}

set.seed(2020)
# generate confounders
calX <- runif(100000, -1, 1)

# calculate expectation of A
expectationA<-mean(gfun(truegamma0+sin(pi*calX)*truegamma))

truebeta0fun<-function(u){
  truebeta0fun<-(-expectationA)*truetreatmentfun(u)
}

epsilonfun<-function(u){
  epsilonfun<-sin(pi*u)/8
}

lgrid<-1001
grid<-seq(0, 1, length.out=lgrid)

rawtruetreatmentvalue<-truetreatmentfun(grid)
truebetavalue<-truebetafun(grid)
truebeta0value<-truebeta0fun(grid)

truemu0inversevalue<-truebeta0value+grid
### linear interpolation 
truemu0inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u!=1){
    truemu0inversefunction<-truemu0inversevalue[tempindex]+(u-grid[tempindex])*(truemu0inversevalue[tempindex+1]-truemu0inversevalue[tempindex])*(lgrid-1)
  } else {
    truemu0inversefunction<-1
  }	
}

mu0fun<-inverse(truemu0inversefunction, lower=0, upper=1)

truemu0<-sapply(grid, mu0fun)


#truemu1inversevalue<-truebeta0value+rawtruetreatmentvalue+grid


truetreatmentvalue<-truetreatmentfun(truemu0)

increasingindicator<-rep(NA, nrep)

IPWbias<-outcomebias<-DRbias<-matrix(NA, nrep, lgrid)
IPWMSE<-outcomeMSE<-DRMSE<-rep(NA, nrep)
IPWmedianbias<-outcomemedianbias<-DRmedianbias<-matrix(NA, nrep,9)
DRcoverage<-rep(NA, nrep)


for (s in 1:nrep){
  set.seed(2020+1000*s)
  
  # generate confounders
  X <- runif(n, -1, 1)
  
  temp<-truegamma0+sin(pi*X)*truegamma
  probability<-gfun(temp)
  A<-rep(NA, n)
  Yinverse<-matrix(NA, n, lgrid)
  
  for (i in 1:n){
    A[i]<-rbinom(1,1,probability[i])
    individualgrid<-runif(1001)
    epsilonvalue<-runif(1,-0.5,0.5)*epsilonfun(individualgrid)
    #tempY<-truebeta0value+A[i]*rawtruetreatmentvalue+X[i]*truebetavalue+epsilonvalue
    tempY<-truebeta0fun(individualgrid)+A[i]*truetreatmentfun(individualgrid)+sin(pi*X[i])*truebetafun(individualgrid)++epsilonvalue
    individualhvalue<-tempY+individualgrid
    rval<-approxfun(c(0, individualgrid,1), c(0,individualhvalue, 1),
                    method = "constant")
    hvalue<-sapply(grid, rval) 
    ### G_i inverse
    #rawYinverse[i, ]<-hvalue	
    #Yinverse[i,]<-rawYinverse[i, ]
    #Yinverse[i,which(hvalue<0)]<-0
    #Yinverse[i,which(hvalue>1)]<-1
    Yinverse[i,]<-hvalue
  }
  
  ### G bar inverse
  #Yinversemean<-apply(Yinverse, 2, mean)
  
  ### linear interpolation 
  #Gbarinversefunction<-function(u){
  #	tempindex<-floor(u*(lgrid-1))+1
  #if (u!=1){
  #	Gbarinversefunction<-Yinversemean[tempindex]+(u-grid[tempindex])*(Yinversemean[tempindex+1]-Yinversemean[tempindex])*(lgrid-1)
  #} else {
  #	Gbarinversefunction<-1
  #}	
  #}
  
  ### G bar function
  #Gbarfun<-inverse(Gbarinversefunction, lower=0, upper=1)
  
  ### G bar，Estimated Frechet mean (all checked correct, see testhatmu.r)
  #hatmu<-sapply(grid, Gbarfun)
  
  
  identity<-grid
  
  # fit smoothing spline propensity score model
  binaryspline.fit <- gssanova(A~X,family="binomial")
  
  est <- binaryspline.fit$eta
  pihat <- exp(est)/(1+exp(est))
  
  
  if (min(pihat)>0 & max(pihat)<1){
    
    estimatemu0LY<-estimatemu0IPWcomponent1<-estimatemu0IPWcomponent2<-matrix(NA, n, lgrid)
    for (i in 1:n){
      estimatemu0Yinversefunction<-function(u){
        tempindex<-floor(u*(lgrid-1))+1
        if (u>0 & u<1){
          estimatemu0Yinversefunction<-Yinverse[i,tempindex]+(u-grid[tempindex])*(Yinverse[i, tempindex+1]-Yinverse[i, tempindex])*(lgrid-1)
        } else if (u<=0){
          estimatemu0Yinversefunction<-0
        } else {
          estimatemu0Yinversefunction<-1
        }	
      }
      
      estimatemu0LY[i,]<-sapply(identity, estimatemu0Yinversefunction)
      
      estimatemu0IPWcomponent1[i,]<-A[i]*estimatemu0LY[i,]/pihat[i]
      estimatemu0IPWcomponent2[i,]<-(1-A[i])*estimatemu0LY[i,]/(1-pihat[i])
    }
    
    ### IPW estimate
    estimatemu0IPWpsi1<-apply(estimatemu0IPWcomponent1, 2, mean)
    estimatemu0IPWpsi0<-apply(estimatemu0IPWcomponent2, 2, mean)
    
    #IPWmu0inverse<-estimatemu0IPWpsi0
    
    
    #plot(IPWtxeffect, type="l")
    
    ### regression estimate
    #estimatemu0lmfit<-lm(estimatemu0LY~A+X)
    #estimatemu0outcometxeffect<-estimatemu0lmfit$coefficients[2,]
    
    
    # Fit outcome adaptive propensity score model
    treatmentindex<-which(A==1)
    controlindex<-which(A==0)
    
    estimatemu0m1fit<-estimatemu0m0fit<-matrix(NA, n, lgrid)
    for (j in 1:lgrid){
      splinefit<-smooth.spline(estimatemu0LY[treatmentindex,j]~X[treatmentindex])
      if (class(splinefit)!="try-error"){
        temp1<-predict(splinefit, X)$y
        estimatemu0m1fit[,j]<-temp1
      } else {
        estimatemu0m1fit[,j]<-NA 
      }
      splinefit<-try(smooth.spline(estimatemu0LY[controlindex,j]~X[controlindex]))
      if (class(splinefit)!="try-error"){
        temp0<-predict(splinefit, X)$y
        estimatemu0m0fit[,j]<-temp0
      } else {
        estimatemu0m0fit[,j]<-NA
      }
    }
    
    if (sum(is.na(estimatemu0m1fit))+sum(is.na(estimatemu0m0fit))==0){
      
      ###### check estimatemu0m1fit and estimatemu0m0fit
      
      
      # newdesign1<-cbind(rep(1,n), rep(1, n), X)
      # newdesign0<-cbind(rep(1,n), rep(0, n), X)
      
      # estimatemu0m1fit<-newdesign1%*%estimatemu0lmfit$coefficients
      # estimatemu0m0fit<-newdesign0%*%estimatemu0lmfit$coefficients
      
      estimatemu0outcomepsi1<-apply(estimatemu0m1fit, 2, mean)
      estimatemu0outcomepsi0<-apply(estimatemu0m0fit, 2, mean)
      
      # outcomemu0inverse<-estimatemu0outcomepsi0
      
      estimatemu0DRcomponent1<-estimatemu0DRcomponent2<-matrix(NA, n, lgrid)
      for (i in 1:n){
        estimatemu0DRcomponent1[i, ]<-estimatemu0IPWcomponent1[i,]-(A[i]/pihat[i]-1)*estimatemu0m1fit[i,]
        estimatemu0DRcomponent2[i, ]<-estimatemu0IPWcomponent2[i,]-((1-A[i])/(1-pihat[i])-1)*estimatemu0m0fit[i,]
      }
      
      estimatemu0DRpsi1<-apply(estimatemu0DRcomponent1, 2, mean)
      estimatemu0DRpsi0<-apply(estimatemu0DRcomponent2, 2, mean)
      
      DRmu0inverse<-estimatemu0DRpsi0
      
      ### linear interpolation 
      DRmu0inversefunction<-function(u){
        tempindex<-floor(u*(lgrid-1))+1
        if (u>0 & u<1){
          DRmu0inversefunction<-DRmu0inverse[tempindex]+(u-grid[tempindex])*(DRmu0inverse[tempindex+1]-DRmu0inverse[tempindex])*(lgrid-1)
        } else if (u<=0) {
          DRmu0inversefunction<-0
        } else {
          DRmu0inversefunction<-1
        }	
      }
      
      ### G bar function, mu0 function
      DRmu0inversefun<-inverse(DRmu0inversefunction, lower=0, upper=1)
      
      ### G bar，Estimated Frechet mean (all checked correct, see testhatmu.r)
      DRhatmu0<-sapply(grid, DRmu0inversefun)
      
      
      DRLY<-DRIPWcomponent1<-DRIPWcomponent2<-matrix(NA, n, lgrid)
      for (i in 1:n){
        Yinversefunction<-function(u){
          tempindex<-floor(u*(lgrid-1))+1
          if (u!=1){
            Yinversefunction<-Yinverse[i,tempindex]+(u-grid[tempindex])*(Yinverse[i, tempindex+1]-Yinverse[i, tempindex])*(lgrid-1)
          } else {
            Yinversefunction<-1
          }	
        }
        
        DRLY[i,]<-sapply(DRhatmu0, Yinversefunction)
        
        DRIPWcomponent1[i,]<-A[i]*DRLY[i,]/pihat[i]
        DRIPWcomponent2[i,]<-(1-A[i])*DRLY[i,]/(1-pihat[i])
      }
      
      
      m1fit<-m0fit<-matrix(NA, n, lgrid)
      for (j in 1:lgrid){
        splinefit<-smooth.spline(DRLY[treatmentindex,j]~X[treatmentindex])
        if (class(splinefit)!="try-error"){
          temp1<-predict(splinefit, X)$y
          m1fit[,j]<-temp1
        } else {
          m1fit[,j]<-NA 
        }
        splinefit<-smooth.spline(DRLY[controlindex,j]~X[controlindex])
        if (class(splinefit)!="try-error"){
          temp0<-predict(splinefit, X)$y
          m0fit[,j]<-temp0
        } else {
          m0fit[,j]<-NA 
        }
      }
      
      if (sum(is.na(m1fit))+sum(is.na(m0fit))==0){
        
        DRcomponent1<-DRcomponent2<-matrix(NA, n, lgrid)
        for (i in 1:n){
          DRcomponent1[i, ]<-DRIPWcomponent1[i,]-(A[i]/pihat[i]-1)*m1fit[i,]
          DRcomponent2[i, ]<-DRIPWcomponent2[i,]-((1-A[i])/(1-pihat[i])-1)*m0fit[i,]
        }
        
        DRtxeffect<-apply(DRcomponent1, 2, mean)-apply(DRcomponent2, 2, mean)
        
        
        DRpsi1<-apply(DRcomponent1, 2, mean)
        DRpsi0<-apply(DRcomponent2, 2, mean)
        
        DRindividual<-DRcomponent1-DRcomponent2
        
        calsequence<-seq(5, lgrid-1, 10)
        
        calDRindividual<-DRindividual[,calsequence]
        
        covcalDRindividual<-cov(calDRindividual)
        meancalDRindividual<-apply(calDRindividual, 2, mean)
        
        generatenormalvariable<-mvrnorm(1000, rep(0, length(calsequence)), covcalDRindividual)
        
        maxgeneratenormalvariable<-apply(abs(generatenormalvariable), 1, max)
        sortmaxgeneratenormalvariable<-sort(maxgeneratenormalvariable)
        DRband<-sortmaxgeneratenormalvariable[950]/sqrt(n)
        
        DRhigherband<-DRtxeffect+DRband
        DRlowerband<-DRtxeffect-DRband
        
        
        highercoverageindicator<-sum(DRhigherband<truetreatmentvalue)
        lowercoverageindicator<-sum(DRlowerband>truetreatmentvalue)
        
        if (sum(highercoverageindicator+lowercoverageindicator)==0){
          DRcoverage[s]<-1
        } else {
          DRcoverage[s]<-0
        }
        
        
        
        DRtransportvalue<-sapply(truemu0, DRmu0inversefunction)
        
        DRtxeffectfunction<-function(u){
          tempindex<-floor(u*(lgrid-1))+1
          if (u>0 & u<1){
            DRtxeffectfunction<-DRtxeffect[tempindex]+(u-grid[tempindex])*(DRtxeffect[tempindex+1]-DRtxeffect[tempindex])*(lgrid-1)
          } else {
            DRtxeffectfunction<-0
          }
        }
        
        DRparralleleffect<-sapply(DRtransportvalue, DRtxeffectfunction)
        
        DRbias[s, ]<-DRparralleleffect-truetreatmentvalue
        
        DRMSE[s]<-sqrt(mean((DRbias[s,])^2))
        
        tempsequence<-seq(101, 901, by=100)
        
        DRmedianbias[s,]<-DRparralleleffect[tempsequence]-truetreatmentvalue[tempsequence]
      } else {
        DRcoverage[s]<-NA
        DRbias[s,]<-rep(NA, lgrid)
        DRMSE[s]<-NA
        DRmedianbias[s,]<-rep(NA, 9)
      } 
    } else {
      DRcoverage[s]<-NA
      DRbias[s,]<-rep(NA, lgrid)
      DRMSE[s]<-NA
      DRmedianbias[s,]<-rep(NA, 9)
    } 
  } else{
    DRcoverage[s]<-NA
    DRbias[s,]<-rep(NA, lgrid)
    DRMSE[s]<-NA
    DRmedianbias[s,]<-rep(NA, 9)
  }
  #plot(truetreatmentvalue)
  #lines(outcomeparralleleffect)
  #lines(DRparralleleffect)
  #lines(IPWparralleleffect)
}


save(DRMSE, file=paste("W2hatmu0referencenonparametricsmoothsplineDRMSE_n", n, "_p", p, ".dat",sep=""))

save(DRmedianbias, file=paste("W2hatmu0referencenonparametricsmoothsplineDRmedianbias_n", n, "_p", p, ".dat",sep=""))

save(DRcoverage, file=paste("W2hatmu0referencenonparametricsmoothsplineDRcoverage_n", n, "_p", p, ".dat",sep=""))



























