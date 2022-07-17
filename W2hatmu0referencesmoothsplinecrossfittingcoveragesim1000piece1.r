### This code is the data adaptive cross fitting estimator for Scenario 1, 
### where both the propensity score model and the outcome regression model are linear functions of X


rm(list=ls())
library(MASS)
library(GoFKernel)
library(gss)
library(dplyr)

# Sample size
n<-1000

# Number of MC
nrep<-50

# piece number
piece<-1

### median number
mediannumber<-100

### folds number
K<-5

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
expectationA<-mean(gfun(truegamma0+calX*truegamma))

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


implement<-function(s){
set.seed(2020+1000*s)

# generate confounders
X <- runif(n, -1, 1)

temp<-truegamma0+X*truegamma
probability<-gfun(temp)
A<-rep(NA, n)
Yinverse<-matrix(NA, n, lgrid)

for (i in 1:n){
  A[i]<-rbinom(1,1,probability[i])
  individualgrid<-runif(1001)
  epsilonvalue<-runif(1,-0.5,0.5)*epsilonfun(individualgrid)
  #tempY<-truebeta0value+A[i]*rawtruetreatmentvalue+X[i]*truebetavalue+epsilonvalue
  tempY<-truebeta0fun(individualgrid)+A[i]*truetreatmentfun(individualgrid)+X[i]*truebetafun(individualgrid)++epsilonvalue
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


# Fit outcome adaptive propensity score model
treatmentindex<-which(A==1)
controlindex<-which(A==0)

estimatemu0m1fit<-estimatemu0m0fit<-matrix(NA, n, lgrid)
for (j in 1:lgrid){
  splinefit<-smooth.spline(estimatemu0LY[treatmentindex,j]~X[treatmentindex])
  temp1<-predict(splinefit, X)$y
  estimatemu0m1fit[,j]<-temp1
  splinefit<-try(smooth.spline(estimatemu0LY[controlindex,j]~X[controlindex]))
  temp0<-predict(splinefit, X)$y
  estimatemu0m0fit[,j]<-temp0
}

estimatemu0outcomepsi1<-apply(estimatemu0m1fit, 2, mean)
estimatemu0outcomepsi0<-apply(estimatemu0m0fit, 2, mean)  


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
  temp1<-predict(splinefit, X)$y
  m1fit[,j]<-temp1
  splinefit<-smooth.spline(DRLY[controlindex,j]~X[controlindex])
  temp0<-predict(splinefit, X)$y
  m0fit[,j]<-temp0
}

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


LY<-matrix(NA, n, lgrid)
for (i in 1:n){
  Yinversefunction<-function(u){
    tempindex<-floor(u*(lgrid-1))+1
    if (u!=1){
      Yinversefunction<-Yinverse[i,tempindex]+(u-grid[tempindex])*(Yinverse[i, tempindex+1]-Yinverse[i, tempindex])*(lgrid-1)
    } else {
      Yinversefunction<-1
    }	
  }
  
  LY[i,]<-sapply(identity, Yinversefunction)
}


####### Remember to set seed
#######
#######
mu0inverseallmedianestimate<-matrix(NA, mediannumber, lgrid)
for (l in 1:mediannumber){
  
  folds <- cut(seq(1,n),breaks=K,labels=FALSE)
  randomindex<-sample(n,n)
  randomfolds<-folds[randomindex]
  
  savemu0inverse<-matrix(NA, K, lgrid)
  
  for (j in 1:K){
    portion2index<-which(randomfolds==j)
    portion1index<-setdiff(1:n, portion2index)	
    
    tempX1<-X[portion1index]
    tempX2<-X[portion2index]
    
    # fit smoothing spline propensity score model
    estimatemu0binaryspline.fit <- gssanova(A[portion1index]~tempX1,family="binomial")
    
    
    ###### problems here
    #est <-  estimatemu0binaryspline.fit$eta
    
    estimatemu0portion2est<-predict(estimatemu0binaryspline.fit, data.frame(tempX1=tempX2))
    estimatemu0pihat2 <- exp(est)/(1+exp(est))
    
    estimatemu0IPWcomponent1<-estimatemu0IPWcomponent2<-matrix(NA, n/K, lgrid)
    for (i in 1:(n/K)){
      subjectindex<-portion2index[i]
      
      estimatemu0IPWcomponent1[i,]<-A[subjectindex]*LY[subjectindex,]/estimatemu0pihat2[i]
      estimatemu0IPWcomponent2[i,]<-(1-A[subjectindex])*LY[subjectindex,]/(1-estimatemu0pihat2[i])
    }
    
    
    tempA<-A[portion1index]
    # Fit outcome adaptive propensity score model
    temptreatmentindex<-which(tempA==1)
    tempcontrolindex<-which(tempA==0)
    
    estimatemu0m1fit<-estimatemu0m0fit<-matrix(NA, n/K, lgrid)  
    for (kk in 1:lgrid){
      splinefit<-smooth.spline(LY[temptreatmentindex,kk]~X[temptreatmentindex])
      temp1<-predict(splinefit, X[portion2index])$y
      estimatemu0m1fit[,kk]<-temp1
      splinefit<-smooth.spline(LY[tempcontrolindex,kk]~X[tempcontrolindex])
      temp0<-predict(splinefit, X[portion2index])$y
      estimatemu0m0fit[,kk]<-temp0
    }
    
    estimatemu0DRcomponent1<-estimatemu0DRcomponent2<-matrix(NA, n/K, lgrid)
    for (i in 1:(n/K)){
      subjectindex<-portion2index[i]
      estimatemu0DRcomponent1[i, ]<-estimatemu0IPWcomponent1[i,]-(A[subjectindex]/estimatemu0pihat2[i]-1)*estimatemu0m1fit[i,]
      estimatemu0DRcomponent2[i, ]<-estimatemu0IPWcomponent2[i,]-((1-A[subjectindex])/(1-estimatemu0pihat2[i])-1)*estimatemu0m0fit[i,]
    }
    
    tempmu1inverse<-apply(estimatemu0DRcomponent1, 2, mean)
    tempmu0inverse<-apply(estimatemu0DRcomponent2, 2, mean)
    
    
    
    savemu0inverse[j,]<-tempmu0inverse
    
  }
  
  finalmu0inverse<-apply(savemu0inverse, 2, mean)
  
  mu0inverseallmedianestimate[l,]<-finalmu0inverse
}

naomitmu0inverseallmedianestimate<-na.omit(mu0inverseallmedianestimate)

hatmu0inverse<-apply(naomitmu0inverseallmedianestimate, 2, median)



### linear interpolation 
mu0inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    mu0inversefunction<-hatmu0inverse[tempindex]+(u-grid[tempindex])*(hatmu0inverse[tempindex+1]-hatmu0inverse[tempindex])*(lgrid-1)
  } else if (u<=0) {
    mu0inversefunction<-0
  } else {
    mu0inversefunction<-1
  }	
}

### G bar function
mu0inversefun<-inverse(mu0inversefunction, lower=0, upper=1)

### G bar，Estimated Frechet mean (all checked correct, see testhatmu.r)
hatmu0<-sapply(grid, mu0inversefun)



#increasingindicator[s]<-all(diff(Yinversemean)>=0)

#hatmufunction<-function(u){
#	tempdistance<-abs(u-grid)
#	hatmufunction<-hatmuinverse[which.min(tempdistance)]
#}

#hatmufun<-inverse(hatmuinversefunction, lower=0, upper=1)

# Estimated Frechet mean
#hatmu<-sapply(grid, hatmufun)


LY<-IPWcomponent1<-IPWcomponent2<-matrix(NA, n, lgrid)
for (i in 1:n){
  Yinversefunction<-function(u){
    tempindex<-floor(u*(lgrid-1))+1
    if (u!=1){
      Yinversefunction<-Yinverse[i,tempindex]+(u-grid[tempindex])*(Yinverse[i, tempindex+1]-Yinverse[i, tempindex])*(lgrid-1)
    } else {
      Yinversefunction<-1
    }	
  }
  
  LY[i,]<-sapply(hatmu0, Yinversefunction)
}








####### Remember to set seed
#######
#######
allmedianestimate<-matrix(NA, mediannumber, lgrid)
for (l in 1:mediannumber){
  
  folds <- cut(seq(1,n),breaks=K,labels=FALSE)
  randomindex<-sample(n,n)
  randomfolds<-folds[randomindex]
  
  saveDRtxeffect<-matrix(NA, K, lgrid)
  
  for (j in 1:K){
    portion2index<-which(randomfolds==j)
    portion1index<-setdiff(1:n, portion2index)	
    
    
    tempX1<-X[portion1index]
    tempX2<-X[portion2index]
    
    # fit smoothing spline propensity score model
    binaryspline.fit <- gssanova(A[portion1index]~tempX1,family="binomial")
    
    
    ###### problems here
    #est <-  estimatemu0binaryspline.fit$eta
    
    est<-predict(binaryspline.fit, data.frame(tempX1=tempX2))
    pihat2 <- exp(est)/(1+exp(est))
    
    IPWcomponent1<-IPWcomponent2<-matrix(NA, n/K, lgrid)
    for (i in 1:(n/K)){
      subjectindex<-portion2index[i]
      
      tempLY<-sapply(hatmu0, Yinversefunction)
      
      IPWcomponent1[i,]<-A[subjectindex]*LY[subjectindex,]/pihat2[i]
      IPWcomponent2[i,]<-(1-A[subjectindex])*LY[subjectindex,]/(1-pihat2[i])
    }
    
    
    tempA<-A[portion1index]
    # Fit outcome adaptive propensity score model
    temptreatmentindex<-which(tempA==1)
    tempcontrolindex<-which(tempA==0)
    
    estimatemu0m1fit<-estimatemu0m0fit<-matrix(NA, n/K, lgrid)  
    for (kk in 1:lgrid){
      splinefit<-smooth.spline(LY[temptreatmentindex,kk]~X[temptreatmentindex])
      temp1<-predict(splinefit, X[portion2index])$y
      estimatemu0m1fit[,kk]<-temp1
      splinefit<-try(smooth.spline(LY[tempcontrolindex,j]~X[tempcontrolindex]))
      temp0<-predict(splinefit, X[portion2index])$y
      estimatemu0m0fit[,kk]<-temp0
    }
    
    
    
    
    
    
    
    
    # Fit propensity score model
    logisticfit<-glm(A[portion1index]~X[portion1index], family=binomial(link='logit'))	
    portion2designmatrix<-cbind(rep(1,n/K), X[portion2index])
    pihat2<-exp(portion2designmatrix%*%logisticfit$coefficients)/(1+exp(portion2designmatrix%*%logisticfit$coefficients))
    
    
    
    DRcomponent1<-DRcomponent2<-matrix(NA, n/K, lgrid)
    for (i in 1:(n/K)){
      subjectindex<-portion2index[i]
      DRcomponent1[i, ]<-IPWcomponent1[i,]-(A[subjectindex]/pihat2[i]-1)*m1fit[i,]
      DRcomponent2[i, ]<-IPWcomponent2[i,]-((1-A[subjectindex])/(1-pihat2[i])-1)*m0fit[i,]
    }
    
    DRtxeffect<-apply(DRcomponent1, 2, mean)-apply(DRcomponent2, 2, mean)
    
    saveDRtxeffect[j,]<-DRtxeffect
    
  }
  
  finalDRtxeffect<-apply(saveDRtxeffect, 2, mean)
  
  allmedianestimate[l,]<-finalDRtxeffect
}

naomitallmedianestimate<-na.omit(allmedianestimate)

finalmedianDRtxeffect<-apply(naomitallmedianestimate, 2, median)


DRmedianallmatrix<-list()
DRmedianoperatornorm<-rep(NA, mediannumber)
for (j in 1:mediannumber){
  DRmedianallmatrix[[j]]<-covcalDRindividual+(allmedianestimate[j,calsequence]-finalmedianDRtxeffect[calsequence])%*%t(allmedianestimate[j,calsequence]-finalmedianDRtxeffect[calsequence])
  DRmedianoperatornorm[j]<-norm(DRmedianallmatrix[[j]], type="2")	
}
sortDRmedianoperatornorm<-sort(DRmedianoperatornorm)
DRmedianindex<-which(DRmedianoperatornorm==sortDRmedianoperatornorm[mediannumber/2])


covcalDRmedian<-DRmedianallmatrix[[DRmedianindex]]


DRmediangeneratenormalvariable<-mvrnorm(1000, rep(0, length(calsequence)), covcalDRmedian)

maxDRmediangeneratenormalvariable<-apply(abs(DRmediangeneratenormalvariable), 1, max)
sortmaxDRmediangeneratenormalvariable<-sort(maxDRmediangeneratenormalvariable)
DRmedianband<-sortmaxDRmediangeneratenormalvariable[950]/sqrt(n)

DRmedianhigherband<-finalmedianDRtxeffect+DRmedianband
DRmedianlowerband<-finalmedianDRtxeffect-DRmedianband

highercoverageindicator<-sum(DRmedianhigherband<truetreatmentvalue)
lowercoverageindicator<-sum(DRmedianlowerband>truetreatmentvalue)

if (sum(highercoverageindicator+lowercoverageindicator)==0){
  returnDRmediancoverage<-1
} else {
  returnDRmediancoverage<-0
}

finalmedianDRtxeffectfunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    finalmedianDRtxeffectfunction<-finalmedianDRtxeffect[tempindex]+(u-grid[tempindex])*(finalmedianDRtxeffect[tempindex+1]-finalmedianDRtxeffect[tempindex])*(lgrid-1)
  } else {
    finalmedianDRtxeffectfunction<-0
  }
}

transportvalue<-sapply(truemu0, mu0inversefunction)

finalmedianDRparralleleffect<-sapply(transportvalue, finalmedianDRtxeffectfunction)

returnDRbias<-finalmedianDRparralleleffect-truetreatmentvalue

returnDRMSE<-sqrt(mean((returnDRbias)^2))


tempsequence<-seq(101, 901, by=100)
returnDRmedianbias<-finalmedianDRparralleleffect[tempsequence]-truetreatmentvalue[tempsequence]


mylist=list(mediancoverage=returnDRmediancoverage, bias=returnDRbias, MSE=returnDRMSE, medianbias=returnDRmedianbias)
return(mylist)
}



increasingindicator<-rep(NA, nrep)

DRbias<-matrix(NA, nrep, lgrid)
DRMSE<-rep(NA, nrep)
DRmedianbias<-matrix(NA, nrep,9)
DRmediancoverage<-rep(NA, nrep)


for (s in 1:nrep){
  runprogram<-try(implement(s+nrep*piece-nrep), silent=F)
  if (class(runprogram)!="try-error"){
  DRmediancoverage[s]<-runprogram$mediancoverage
  DRbias[s, ]<-runprogram$bias
  DRMSE[s]<-runprogram$MSE
  DRmedianbias[s,]<-runprogram$medianbias} else {
    DRmediancoverage[s]<-NA
    DRbias[s, ]<-rep(NA, lgrid)
    DRMSE[s]<-NA
    DRmedianbias[s,]<-rep(NA, 9)
  }
}

save(DRMSE, file=paste("W2mu0referencesmoothingsplinecrossfittingsplit", K, "mediannumber", mediannumber,"DRMSE_n", n, "_p", p, "_piece", piece, ".dat",sep=""))

save(DRmedianbias, file=paste("W2mu0referencesmoothingsplinecrossfittingsplit", K, "mediannumber", mediannumber, "DRmedianbias_n", n, "_p", p, "_piece", piece, ".dat",sep=""))

save(DRmediancoverage, file=paste("W2hatmu0referencesmoothingsplinecrossfittingsplit", K, "mediannumber", mediannumber, "DRmediancoverage_n", n, "_p", p, "_piece", piece, ".dat",sep=""))
