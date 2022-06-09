## This code implements our method for the real data application

library(GoFKernel)
library(MASS)

load("/Users/dehankong/Documents/Research/causal manifold/code for submission/real data/savedemo.dat")
load(file="/Users/dehankong/Documents/Research/causal manifold/code for submission/real data/remainingidrange1to1000.dat")

load(file="/Users/dehankong/Documents/Research/causal manifold/code for submission/real data/remainingmaritaltreatmentidrange1to1000.dat")

load(file="/Users/dehankong/Documents/Research/causal manifold/code for submission/real data/allquantilefunctionrange1to1000.dat")
load(file="/Users/dehankong/Documents/Research/causal manifold/code for submission/real data/allcdffunctionrange1to1000.dat")

cutvalue<-1000

n<-length(maritaltreatmentid)
confounders<-matrix(NA, n, 2)
maritaltx<-rep(NA, n)
quantilefunction<-matrix(NA, n, ncol(allquantilefunction))
cdffunction<-matrix(NA, n, ncol(allcdffunction))
for (i in 1:n){
  tempid<-maritaltreatmentid[i]
  index<-which(savedemo[,1]==tempid)
  confounders[i, ]<-savedemo[index, c(2,3)]
  ## married 1, unmarried 0
  maritaltx[i]<-as.numeric((savedemo[index, 4]>5 | savedemo[index, 4]<2))	
  index2<-which(remainingid==tempid)
  quantilefunction[i,]<-allquantilefunction[index2,]
  cdffunction[i,]<-allcdffunction[index2,]
}


lgrid<-1001
grid<-seq(0, 1, length.out=lgrid)


identity<-grid

A<-maritaltx

### First column age, second column gender, female is 2. 
X<-confounders

# Fit propensity score model
logisticfit<-glm(A~X, family=binomial(link='logit'))
pihat<-logisticfit$fitted.values

estimatemu0LY<-estimatemu0IPWcomponent1<-estimatemu0IPWcomponent2<-matrix(NA, n, lgrid)
for (i in 1:n){
  estimatemu0Yinversefunction<-function(u){
    tempindex<-floor(u*(lgrid-1))+1
    if (u!=1 & u!=0){
      estimatemu0Yinversefunction<-quantilefunction[i,tempindex]+(u-grid[tempindex])*(quantilefunction[i, tempindex+1]-quantilefunction[i, tempindex])*(lgrid-1)
    } else if (u==0){
      estimatemu0Yinversefunction<-1
    } else {
      estimatemu0Yinversefunction<-cutvalue
    }	
  }
  
  estimatemu0LY[i,]<-sapply(identity, estimatemu0Yinversefunction)
  
  estimatemu0IPWcomponent1[i,]<-A[i]*estimatemu0LY[i,]/pihat[i]
  estimatemu0IPWcomponent2[i,]<-(1-A[i])*estimatemu0LY[i,]/(1-pihat[i])
}

### IPW estimate
estimatemu0IPWpsi1<-apply(estimatemu0IPWcomponent1, 2, mean)
estimatemu0IPWpsi0<-apply(estimatemu0IPWcomponent2, 2, mean)

IPWmu0inverse<-estimatemu0IPWpsi0


### regression estimate
estimatemu0lmfit<-lm(estimatemu0LY~A+X)
estimatemu0outcometxeffect<-estimatemu0lmfit$coefficients[2,]

newdesign1<-cbind(rep(1,n), rep(1, n), X)
newdesign0<-cbind(rep(1,n), rep(0, n), X)

estimatemu0m1fit<-newdesign1%*%estimatemu0lmfit$coefficients
estimatemu0m0fit<-newdesign0%*%estimatemu0lmfit$coefficients

estimatemu0outcomepsi1<-apply(estimatemu0m1fit, 2, mean)
estimatemu0outcomepsi0<-apply(estimatemu0m0fit, 2, mean)

outcomemu0inverse<-estimatemu0outcomepsi0
outcomemu1inverse<-estimatemu0outcomepsi1

estimatemu0DRcomponent1<-estimatemu0DRcomponent2<-matrix(NA, n, lgrid)
for (i in 1:n){
  estimatemu0DRcomponent1[i, ]<-estimatemu0IPWcomponent1[i,]-(A[i]/pihat[i]-1)*estimatemu0m1fit[i,]
  estimatemu0DRcomponent2[i, ]<-estimatemu0IPWcomponent2[i,]-((1-A[i])/(1-pihat[i])-1)*estimatemu0m0fit[i,]
}

estimatemu0DRpsi1<-apply(estimatemu0DRcomponent1, 2, mean)
estimatemu0DRpsi0<-apply(estimatemu0DRcomponent2, 2, mean)

DRmu0inverse<-estimatemu0DRpsi0
DRmu1inverse<-estimatemu0DRpsi1


### linear interpolation 
IPWmu0inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    IPWmu0inversefunction<-min(IPWmu0inverse[tempindex]+(u-grid[tempindex])*(IPWmu0inverse[tempindex+1]-IPWmu0inverse[tempindex])*(lgrid-1),cutvalue)
  } else if (u<=0) {
    IPWmu0inversefunction<-1
  } else {
    IPWmu0inversefunction<-cutvalue
  }	
}

### mu0 function
IPWmu0inversefun<-inverse(IPWmu0inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu0
IPWhatmu0<-sapply(1:cutvalue, IPWmu0inversefun)

### linear interpolation 
outcomemu0inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    outcomemu0inversefunction<-outcomemu0inverse[tempindex]+(u-grid[tempindex])*(outcomemu0inverse[tempindex+1]-outcomemu0inverse[tempindex])*(lgrid-1)
  } else if (u<=0) {
    outcomemu0inversefunction<-1
  } else {
    outcomemu0inversefunction<-cutvalue
  }	
}

### mu0 function
outcomemu0inversefun<-inverse(outcomemu0inversefunction, lower=0, upper=1)

### G bar，Estimated Frechet mean hatmu0
outcomehatmu0<-sapply(1:cutvalue, outcomemu0inversefun)


### linear interpolation 
outcomemu1inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    outcomemu1inversefunction<-outcomemu1inverse[tempindex]+(u-grid[tempindex])*(outcomemu1inverse[tempindex+1]-outcomemu1inverse[tempindex])*(lgrid-1)
  } else if (u<=0) {
    outcomemu1inversefunction<-1
  } else {
    outcomemu1inversefunction<-cutvalue
  }	
}

### mu1 function
outcomemu1inversefun<-inverse(outcomemu1inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu1
outcomehatmu1<-sapply(1:cutvalue, outcomemu1inversefun)


### linear interpolation 
DRmu0inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    DRmu0inversefunction<-DRmu0inverse[tempindex]+(u-grid[tempindex])*(DRmu0inverse[tempindex+1]-DRmu0inverse[tempindex])*(lgrid-1)
  } else if (u<=0) {
    DRmu0inversefunction<-1
  } else {
    DRmu0inversefunction<-cutvalue
  }	
}

### mu0 function
DRmu0inversefun<-inverse(DRmu0inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu0
DRhatmu0<-sapply(1:cutvalue, DRmu0inversefun)


### linear interpolation 
DRmu1inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    DRmu1inversefunction<-DRmu1inverse[tempindex]+(u-grid[tempindex])*(DRmu1inverse[tempindex+1]-DRmu1inverse[tempindex])*(lgrid-1)
  } else if (u<=0) {
    DRmu1inversefunction<-1
  } else {
    DRmu1inversefunction<-cutvalue
  }	
}

### mu1 function
DRmu1inversefun<-inverse(DRmu1inversefunction, lower=0, upper=1)

### G bar，Estimated Frechet mean hatmu1
DRhatmu1<-sapply(1:cutvalue, DRmu1inversefun)




mediannumber<-100
K<-5

mu0inverseallmedianestimate<-mu1inverseallmedianestimate<-matrix(NA, mediannumber, lgrid)
for (l in 1:mediannumber){
  set.seed(2021+1000*l)
  folds <- cut(seq(1,n),breaks=K,labels=FALSE)
  randomindex<-sample(n,n)
  randomfolds<-folds[randomindex]
  
  savemu0inverse<-matrix(NA, K, lgrid)
  savemu1inverse<-matrix(NA, K, lgrid)
  
  for (j in 1:K){
    portion2index<-which(randomfolds==j)
    portion1index<-setdiff(1:n, portion2index)	
    
    # Fit propensity score model
    estimatemu0logisticfit<-glm(A[portion1index]~X[portion1index], family=binomial(link='logit'))	
    estimatemu0portion2designmatrix<-cbind(rep(1,length(portion2index)), X[portion2index])
    estimatemu0pihat2<-exp(estimatemu0portion2designmatrix%*%estimatemu0logisticfit$coefficients)/(1+exp(estimatemu0portion2designmatrix%*%estimatemu0logisticfit$coefficients))
    
    
    
    estimatemu0IPWcomponent1<-estimatemu0IPWcomponent2<-matrix(NA, length(portion2index), lgrid)
    for (i in 1:length(portion2index)){
      subjectindex<-portion2index[i]
      estimatemu0IPWcomponent1[i,]<-A[subjectindex]*estimatemu0LY[subjectindex,]/estimatemu0pihat2[i]
      estimatemu0IPWcomponent2[i,]<-(1-A[subjectindex])*estimatemu0LY[subjectindex,]/(1-estimatemu0pihat2[i])
    }
    
    ### regression estimate
    estimatemu0lmfit<-lm(estimatemu0LY[portion1index,]~A[portion1index]+X[portion1index])
    
    estimatemu0newdesign1<-cbind(rep(1,length(portion2index)), rep(1, length(portion2index)), X[portion2index])
    estimatemu0newdesign0<-cbind(rep(1,length(portion2index)), rep(0, length(portion2index)), X[portion2index])
    
    estimatemu0m1fit<-estimatemu0newdesign1%*%estimatemu0lmfit$coefficients
    estimatemu0m0fit<-estimatemu0newdesign0%*%estimatemu0lmfit$coefficients
    
    estimatemu0DRcomponent1<-estimatemu0DRcomponent2<-matrix(NA, length(portion2index), lgrid)
    for (i in 1:length(portion2index)){
      subjectindex<-portion2index[i]
      estimatemu0DRcomponent1[i, ]<-estimatemu0IPWcomponent1[i,]-(A[subjectindex]/estimatemu0pihat2[i]-1)*estimatemu0m1fit[i,]
      estimatemu0DRcomponent2[i, ]<-estimatemu0IPWcomponent2[i,]-((1-A[subjectindex])/(1-estimatemu0pihat2[i])-1)*estimatemu0m0fit[i,]
    }
    
    tempmu1inverse<-apply(estimatemu0DRcomponent1, 2, mean)
    tempmu0inverse<-apply(estimatemu0DRcomponent2, 2, mean)
    
    
    
    savemu0inverse[j,]<-tempmu0inverse
    savemu1inverse[j,]<-tempmu1inverse
    
  }
  
  finalmu0inverse<-apply(savemu0inverse, 2, mean)
  finalmu1inverse<-apply(savemu1inverse, 2, mean)
  
  mu0inverseallmedianestimate[l,]<-finalmu0inverse
  mu1inverseallmedianestimate[l,]<-finalmu1inverse
}

DRmedianhatmu0inverse<-apply(mu0inverseallmedianestimate, 2, median)
DRmedianhatmu1inverse<-apply(mu1inverseallmedianestimate, 2, median)


### linear interpolation 
DRmedianmu0inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    DRmedianmu0inversefunction<-max(min(DRmedianhatmu0inverse[tempindex]+(u-grid[tempindex])*(DRmedianhatmu0inverse[tempindex+1]-DRmedianhatmu0inverse[tempindex])*(lgrid-1), cutvalue),1)
  } else if (u<=0) {
    DRmedianmu0inversefunction<-1
  } else {
    DRmedianmu0inversefunction<-cutvalue
  }	
}

### mu0 function
DRmedianmu0inversefun<-inverse(DRmedianmu0inversefunction, lower=0, upper=1)

### G bar，Estimated Frechet mean hatmu0
DRmedianhatmu0<-sapply(1:cutvalue, DRmedianmu0inversefun)


### linear interpolation 
DRmedianmu1inversefunction<-function(u){
  tempindex<-floor(u*(lgrid-1))+1
  if (u>0 & u<1){
    DRmedianmu1inversefunction<-max(min(DRmedianhatmu1inverse[tempindex]+(u-grid[tempindex])*(DRmedianhatmu1inverse[tempindex+1]-DRmedianhatmu1inverse[tempindex])*(lgrid-1), cutvalue),1)
  } else if (u<=0) {
    DRmedianmu1inversefunction<-1
  } else {
    DRmedianmu1inversefunction<-cutvalue
  }	
}

### mu1 function
DRmedianmu1inversefun<-inverse(DRmedianmu1inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu1
DRmedianhatmu1<-sapply(1:cutvalue, DRmedianmu1inversefun)




pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/DRhatmu1mu0range1to1000.pdf")
plot(1:cutvalue,DRhatmu1,type="l",col="black", lty=1, lwd=3, ylim=c(0,1), xlab="", ylab="", cex.axis=2)
lines(1:cutvalue,DRhatmu0, col="blue",lty=2,lwd=3)
dev.off()


pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/outcomehatmu1mu0range1to1000.pdf")
plot(1:cutvalue,outcomehatmu1,type="l",col="black", lty=1, lwd=3, ylim=c(0,1), xlab="", ylab="", cex.axis=2)
lines(1:cutvalue,outcomehatmu0, col="blue",lty=2,lwd=3)
dev.off()

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/DRmedianhatmu1mu0range1to1000.pdf")
plot(1:cutvalue,DRmedianhatmu1,type="l",col="black", lty=1, lwd=3, ylim=c(0,1), xlab="", ylab="", cex.axis=2)
lines(1:cutvalue,DRmedianhatmu0, col="blue",lty=2,lwd=3)
dev.off()


### hatmu0 reference

domaingrid<-cutvalue

A<-maritaltx

### First column age, second column gender, female is 2. 
X<-confounders

# Fit propensity score model
logisticfit<-glm(A~X, family=binomial(link='logit'))
pihat<-logisticfit$fitted.values

IPWLY<-outcomeLY<-DRLY<-DRmedianLY<-IPWcomponent1<-IPWcomponent2<-DRIPWcomponent1<-DRIPWcomponent2<-matrix(NA, n, domaingrid)
for (i in 1:n){
  Yinversefunction<-function(u){
    tempindex<-floor(u*(lgrid-1))+1
    if (u!=1 & u!=0){
      Yinversefunction<-quantilefunction[i,tempindex]+(u-grid[tempindex])*(quantilefunction[i, tempindex+1]-quantilefunction[i, tempindex])*(lgrid-1)
    } else if (u==0){
      Yinversefunction<-1
    } else {
      Yinversefunction<-cutvalue
    }	
  }
  
  IPWLY[i,]<-sapply(IPWhatmu0, Yinversefunction)
  outcomeLY[i,]<-sapply(outcomehatmu0, Yinversefunction)
  DRLY[i,]<-sapply(DRhatmu0, Yinversefunction)
  DRmedianLY[i,]<-sapply(DRmedianhatmu0, Yinversefunction)
  
  IPWcomponent1[i,]<-A[i]*IPWLY[i,]/pihat[i]
  IPWcomponent2[i,]<-(1-A[i])*IPWLY[i,]/(1-pihat[i])
  
  DRIPWcomponent1[i,]<-A[i]*DRLY[i,]/pihat[i]
  DRIPWcomponent2[i,]<-(1-A[i])*DRLY[i,]/(1-pihat[i])
}

### IPW estimate
IPWpsi1<-apply(IPWcomponent1, 2, mean)
IPWpsi0<-apply(IPWcomponent2, 2, mean)

IPWtxeffect<-IPWpsi1-IPWpsi0

### regression estimate
outcomelmfit<-lm(outcomeLY~A+X)
outcometxeffect<-outcomelmfit$coefficients[2,]

outcomeband<-rep(NA, ncol(outcomeLY))
for (k in 1:ncol(outcomeLY)){
  templmfit<-lm(outcomeLY[,k]~A+X)
  summarytemplmfit<-summary(templmfit)
  outcomeband[k]<-summarytemplmfit$coefficients[2,2]*1.96
}

DRlmfit<-lm(DRLY~A+X)

newdesign1<-cbind(rep(1,n), rep(1, n), X)
newdesign0<-cbind(rep(1,n), rep(0, n), X)

m1fit<-newdesign1%*%DRlmfit$coefficients
m0fit<-newdesign0%*%DRlmfit$coefficients

DRcomponent1<-DRcomponent2<-matrix(NA, n, domaingrid)
for (i in 1:n){
  DRcomponent1[i, ]<-DRIPWcomponent1[i,]-(A[i]/pihat[i]-1)*m1fit[i,]
  DRcomponent2[i, ]<-DRIPWcomponent2[i,]-((1-A[i])/(1-pihat[i])-1)*m0fit[i,]
}

DRtxeffect<-apply(DRcomponent1, 2, mean)-apply(DRcomponent2, 2, mean)

DRpsi1<-apply(DRcomponent1, 2, mean)
DRpsi0<-apply(DRcomponent2, 2, mean)

DRindividual<-DRcomponent1-DRcomponent2

calsequence<-seq(5, cutvalue-5, 10)

calDRindividual<-DRindividual[,calsequence]

covcalDRindividual<-cov(calDRindividual)
meancalDRindividual<-apply(calDRindividual, 2, mean)

generatenormalvariable<-mvrnorm(1000, rep(0, length(calsequence)), covcalDRindividual)

maxgeneratenormalvariable<-apply(abs(generatenormalvariable), 1, max)
sortmaxgeneratenormalvariable<-sort(maxgeneratenormalvariable)
DRband<-sortmaxgeneratenormalvariable[950]/sqrt(n)


#### calculate Wasserstein distance between DR mu1 and mu0
rawdensityhatmu0<-(DRhatmu0[2:1000]-DRhatmu0[1:999])
densityhatmu0<-c(rawdensityhatmu0,0)

W2distance<-(sum(DRtxeffect^2*densityhatmu0))^(1/2)

tempvariable<-rep(NA, 1000)
for (kk in 1:1000){
  tempvariable[kk]<-(sum(DRtxeffect[calsequence]*generatenormalvariable[kk,]*densityhatmu0[calsequence])*10)
}

W2sigma2<-var(tempvariable)

W2asymptoticvariance<-W2sigma2/W2distance/W2distance/n

W2distancehigherCI<-W2distance+1.96*sqrt(W2asymptoticvariance)
W2distancelowerCI<-W2distance-1.96*sqrt(W2asymptoticvariance)

W2distancehigherCI
W2distancelowerCI





allmedianestimate<-matrix(NA, mediannumber, domaingrid)
for (l in 1:mediannumber){
  ###### 5-fold cross-fitting 100 random splits
  set.seed(2021+1000*l)
  K<-5
  folds <- cut(seq(1,n),breaks=K,labels=FALSE)
  randomindex<-sample(n,n)
  randomfolds<-folds[randomindex]
  
  saveDRtxeffect<-matrix(NA, K, domaingrid)
  
  for (j in 1:K){
    portion2index<-which(randomfolds==j)
    portion1index<-setdiff(1:n, portion2index)	
    
    # Fit propensity score model
    logisticfit<-glm(A[portion1index]~X[portion1index], family=binomial(link='logit'))	
    portion2designmatrix<-cbind(rep(1,length(portion2index)), X[portion2index])
    pihat2<-exp(portion2designmatrix%*%logisticfit$coefficients)/(1+exp(portion2designmatrix%*%logisticfit$coefficients))
    
    
    
    
    IPWcomponent1<-IPWcomponent2<-matrix(NA, length(portion2index), domaingrid)
    for (i in 1:length(portion2index)){
      subjectindex<-portion2index[i]
      IPWcomponent1[i,]<-A[subjectindex]*DRmedianLY[subjectindex,]/pihat2[i]
      IPWcomponent2[i,]<-(1-A[subjectindex])*DRmedianLY[subjectindex, ]/(1-pihat2[i])
    }
    
    ### regression estimate
    lmfit<-lm(DRmedianLY[portion1index,]~A[portion1index]+X[portion1index])
    
    newdesign1<-cbind(rep(1,length(portion2index)), rep(1, length(portion2index)), X[portion2index])
    newdesign0<-cbind(rep(1,length(portion2index)), rep(0, length(portion2index)), X[portion2index])
    
    m1fit<-newdesign1%*%lmfit$coefficients
    m0fit<-newdesign0%*%lmfit$coefficients
    
    DRcomponent1<-DRcomponent2<-matrix(NA, length(portion2index), domaingrid)
    for (i in 1:length(portion2index)){
      subjectindex<-portion2index[i]
      DRcomponent1[i, ]<-IPWcomponent1[i,]-(A[subjectindex]/pihat2[i]-1)*m1fit[i,]
      DRcomponent2[i, ]<-IPWcomponent2[i,]-((1-A[subjectindex])/(1-pihat2[i])-1)*m0fit[i,]
    }
    
    tempDRtxeffect<-apply(DRcomponent1, 2, mean)-apply(DRcomponent2, 2, mean)
    
    saveDRtxeffect[j,]<-tempDRtxeffect
    
  }
  
  finalDRtxeffect<-apply(saveDRtxeffect, 2, mean)
  
  allmedianestimate[l,]<-finalDRtxeffect
}

DRmediantxeffect<-apply(allmedianestimate, 2, median)


DRmedianallmatrix<-list()
DRmedianoperatornorm<-rep(NA, mediannumber)
for (j in 1:mediannumber){
  DRmedianallmatrix[[j]]<-covcalDRindividual+(allmedianestimate[j,calsequence]-DRmediantxeffect[calsequence])%*%t(allmedianestimate[j,calsequence]-DRmediantxeffect[calsequence])
  DRmedianoperatornorm[j]<-norm(DRmedianallmatrix[[j]], type="2")	
}
sortDRmedianoperatornorm<-sort(DRmedianoperatornorm)
DRmedianindex<-which(DRmedianoperatornorm==sortDRmedianoperatornorm[50])


covcalDRmedian<-DRmedianallmatrix[[DRmedianindex]]


DRmediangeneratenormalvariable<-mvrnorm(1000, rep(0, length(calsequence)), covcalDRmedian)

maxDRmediangeneratenormalvariable<-apply(abs(DRmediangeneratenormalvariable), 1, max)
sortmaxDRmediangeneratenormalvariable<-sort(maxDRmediangeneratenormalvariable)
DRmedianband<-sortmaxDRmediangeneratenormalvariable[950]/sqrt(n)




save(IPWtxeffect, file="/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaltxIPWtxeffectrange1to1000.dat")
save(outcometxeffect, file="/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaltxoutcometxeffectrange1to1000.dat")
save(DRtxeffect, file="/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaltxDRtxeffectrange1to1000.dat")
save(DRmediantxeffect, file="/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaltxmedianDRtxeffectrange1to1000.dat")


pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaltxfigurerange1to1000.pdf")
plot(1:cutvalue,DRtxeffect,type="l",col="black", lty=1, lwd=3, ylim=c(0, 90), xlab="", ylab="", cex.axis=2)
lines(1:cutvalue, outcometxeffect,col="red",lty=2,lwd=3)
lines(1:cutvalue, DRmediantxeffect,col="blue",lty=3,lwd=3)
dev.off()

DRhigherband<-DRtxeffect+DRband
DRlowerband<-DRtxeffect-DRband

DRtxeffect[200]
DRhigherband[200]
DRlowerband[200]


DRtxeffect[400]
DRhigherband[400]
DRlowerband[400]


DRtxeffect[800]
DRhigherband[800]
DRlowerband[800]


tempx<-1:cutvalue

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritalDRtxCIrange1to1000.pdf")
plot(tempx, DRhigherband, col="white",lty=1,lwd=3, ylim=c(-10,90), xlab="", ylab="", cex.axis=2)
lines(tempx, DRlowerband, col="white",lty=1,lwd=3)
polygon(c(tempx, rev(tempx)), c(DRlowerband, rev(DRhigherband)), col="yellow", border= NA)
lines(tempx,DRtxeffect, col="black", lty=1,  lwd=3)
dev.off()


### DR median
DRmedianhigherband<-DRmediantxeffect+DRmedianband
DRmedianlowerband<-DRmediantxeffect-DRmedianband
pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritalDRmediantxCIrange1to1000.pdf")
plot(tempx, DRmedianhigherband, col="white",lty=4,lwd=3, ylim=c(-10,90), xlab="", ylab="", cex.axis=2)
lines(tempx, DRmedianlowerband, col="white",lty=4,lwd=3)
polygon(c(tempx, rev(tempx)), c(DRmedianlowerband, rev(DRmedianhigherband)), col="yellow", border= NA)
lines(tempx,DRmediantxeffect, col="blue", lty=4,  lwd=3)
dev.off()

outcomehigherband<-outcometxeffect+outcomeband
outcomelowerband<-outcometxeffect-outcomeband

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaloutcometxCIrange1to1000.pdf")
plot(tempx, outcomehigherband, col="white",lty=1,lwd=3, ylim=c(-10,90), xlab="", ylab="", cex.axis=2)
lines(tempx, outcomelowerband, col="white",lty=1,lwd=3)
polygon(c(tempx, rev(tempx)), c(outcomelowerband, rev(outcomehigherband)), col="yellow", border= NA)
lines(tempx,outcometxeffect, col="red",lty=2, lwd=3)
dev.off()

### plot causal transport map

DRtransportmap<-DRtxeffect+1:cutvalue
DRtransportmaphigherband<-DRhigherband+1:cutvalue
DRtransportmaplowerband<-DRlowerband+1:cutvalue


DRmediantransportmap<-DRmediantxeffect+1:cutvalue
DRmediantransportmaphigherband<-DRmedianhigherband+1:cutvalue
DRmediantransportmaplowerband<-DRmedianlowerband+1:cutvalue


outcometransportmap<-outcometxeffect+1:cutvalue
outcometransportmaphigherband<-outcomehigherband+1:cutvalue
outcometransportmaplowerband<-outcomelowerband+1:cutvalue

##### Need to modify later

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritalDRtransportmapfigurerange1to1000.pdf")
plot(tempx, DRtransportmaphigherband, type="l", col="white",lty=1,lwd=3, ylim=c(-10,1010), xlab="", ylab="", cex.axis=2)
lines(tempx, DRtransportmaplowerband, col="white",lty=1,lwd=3)
polygon(c(tempx, rev(tempx)), c(DRtransportmaplowerband, rev(DRtransportmaphigherband)), col="yellow", border= NA)
lines(tempx,DRtransportmap, col="black", lty=1,  lwd=3)
#### Add dashed lines
lines(c(200, 200), c(0, DRtransportmap[200]),lty=2, lwd=2, col="red")
lines(c(400, 400), c(0, DRtransportmap[400]),lty=2, lwd=2, col="red")
lines(c(800, 800), c(0, DRtransportmap[800]),lty=2, lwd=2, col="red")
dev.off()





pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencemaritaltransportmapfigurerange1to1000.pdf")
plot(1:cutvalue,DRtransportmap,type="l",col="black", lty=1, lwd=3, xlab="", ylab="", cex.axis=2)
lines(1:cutvalue, outcometransportmap,col="red",lty=2,lwd=3)
lines(1:cutvalue, DRmediantransportmap,col="blue",lty=3,lwd=3)
lines(1:cutvalue, 1:cutvalue, col="green", lty=1, lwd=3)
dev.off()







### plot delta*mu0^-1 function

### DRtxeffectfunction I to I, 1:1000 to 1:1000

DRtxeffectfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    DRtxeffectfunction<-DRtxeffect[tempindex]+(u-tempindex)*(DRtxeffect[tempindex+1]-DRtxeffect[tempindex])
  } else {
    DRtxeffectfunction<-0
  }	
}


DRdeltamuinverse<-sapply(DRmu0inverse, DRtxeffectfunction)





DRhigherbandfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    DRhigherbandfunction<-DRhigherband[tempindex]+(u-tempindex)*(DRhigherband[tempindex+1]-DRhigherband[tempindex])
  } else {
    DRhigherbandfunction<-0
  }	
}

DRlowerbandfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    DRlowerbandfunction<-DRlowerband[tempindex]+(u-tempindex)*(DRlowerband[tempindex+1]-DRlowerband[tempindex])
  } else {
    DRlowerbandfunction<-0
  }	
}

DRdeltamuinversehigherband<-sapply(DRmu0inverse, DRhigherbandfunction)
DRdeltamuinverselowerband<-sapply(DRmu0inverse, DRlowerbandfunction)

### evaluate at median
DRdeltamuinverse[501]
DRdeltamuinversehigherband[501]
DRdeltamuinverselowerband[501]


### DRmediantxeffectfunction I to I, 1:1000 to 1:1000

DRmediantxeffectfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    DRmediantxeffectfunction<-DRmediantxeffect[tempindex]+(u-tempindex)*(DRmediantxeffect[tempindex+1]-DRmediantxeffect[tempindex])
  } else {
    DRmediantxeffectfunction<-0
  }	
}

DRmediandeltamuinverse<-sapply(DRmedianhatmu0inverse, DRmediantxeffectfunction)


DRmedianhigherbandfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    DRmedianhigherbandfunction<-DRmedianhigherband[tempindex]+(u-tempindex)*(DRmedianhigherband[tempindex+1]-DRmedianhigherband[tempindex])
  } else {
    DRmedianhigherbandfunction<-0
  }	
}

DRmedianlowerbandfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    DRmedianlowerbandfunction<-DRmedianlowerband[tempindex]+(u-tempindex)*(DRmedianlowerband[tempindex+1]-DRmedianlowerband[tempindex])
  } else {
    DRmedianlowerbandfunction<-0
  }	
}

DRmediandeltamuinversehigherband<-sapply(DRmedianhatmu0inverse, DRmedianhigherbandfunction)
DRmediandeltamuinverselowerband<-sapply(DRmedianhatmu0inverse, DRmedianlowerbandfunction)

### outcometxeffectfunction I to I, 1:1000 to 1:1000

outcometxeffectfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    outcometxeffectfunction<-outcometxeffect[tempindex]+(u-tempindex)*(outcometxeffect[tempindex+1]-outcometxeffect[tempindex])
  } else {
    outcometxeffectfunction<-0
  }	
}

outcomedeltamuinverse<-sapply(outcomemu0inverse, outcometxeffectfunction)


outcomehigherbandfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    outcomehigherbandfunction<-outcomehigherband[tempindex]+(u-tempindex)*(outcomehigherband[tempindex+1]-outcomehigherband[tempindex])
  } else {
    outcomehigherbandfunction<-0
  }	
}

outcomelowerbandfunction<-function(u){
  tempindex<-floor(u)
  if (tempindex>0 & tempindex<cutvalue){
    outcomelowerbandfunction<-outcomelowerband[tempindex]+(u-tempindex)*(outcomelowerband[tempindex+1]-outcomelowerband[tempindex])
  } else {
    outcomelowerbandfunction<-0
  }	
}

outcomedeltamuinversehigherband<-sapply(outcomemu0inverse, outcomehigherbandfunction)
outcomedeltamuinverselowerband<-sapply(outcomemu0inverse, outcomelowerbandfunction)



pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencedeltamuinversefigurerange1to1000.pdf")
plot(grid, DRdeltamuinverse,type="l",col="black", lty=1, lwd=3, ylim=c(0, 90), xlab="", ylab="", cex.axis=2)
lines(grid,outcomedeltamuinverse,col="red",lty=2,lwd=3)
lines(grid, DRmediandeltamuinverse, col="blue",lty=3,lwd=3)
dev.off()


tempx<-grid

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencedeltamuinverseoutcomeCIrange1to1000.pdf")
plot(tempx, outcomedeltamuinversehigherband, col="white",lty=1,lwd=3, ylim=c(-10,90), xlab="", ylab="", cex.axis=2)
lines(tempx, outcomedeltamuinverselowerband, col="white",lty=1,lwd=3)
polygon(c(tempx, rev(tempx)), c(outcomedeltamuinverselowerband, rev(outcomedeltamuinversehigherband)), col="yellow", border= NA)
lines(tempx,outcomedeltamuinverse, col="red",lty=2,  lwd=3)
dev.off()

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencedeltamuinverseDRCIrange1to1000.pdf")
plot(tempx, DRdeltamuinversehigherband, col="white",lty=1,lwd=3, ylim=c(-10,90), xlab="", ylab="", cex.axis=2)
lines(tempx, DRdeltamuinverselowerband, col="white",lty=1,lwd=3)
polygon(c(tempx, rev(tempx)), c(DRdeltamuinverselowerband, rev(DRdeltamuinversehigherband)), col="yellow", border= NA)
lines(tempx,DRdeltamuinverse, col="black", lty=1,  lwd=3)
dev.off()

pdf(file = "/Users/dehankong/Documents/Research/causal manifold/revision data/hatmu0referencedeltamuinverseDRmedianCIrange1to1000.pdf")
plot(tempx, DRmediandeltamuinversehigherband, col="white",lty=1,lwd=3, ylim=c(-10,90), xlab="", ylab="", cex.axis=2)
lines(tempx, DRmediandeltamuinverselowerband, col="white",lty=1,lwd=3)
polygon(c(tempx, rev(tempx)), c(DRmediandeltamuinverselowerband, rev(DRmediandeltamuinversehigherband)), col="yellow", border= NA)
lines(tempx,DRmediandeltamuinverse, col="blue",lty=3,  lwd=3)
dev.off()
