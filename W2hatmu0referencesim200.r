### The code implements our method (IPW, outcome regression and doubly robust estimators) for the simulation study in the paper when both models are correctly specified.   

rm(list=ls())
library(MASS)
library(GoFKernel)

# Sample size
n<-200

# Number of MC
nrep<-1000

# Dimension of confounders
p<-1

# Number of samples for empirical CDF for each subject
m<-1000


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


truetreatmentvalue<-truetreatmentfun(truemu0)

increasingindicator<-rep(NA, nrep)

IPWbias<-outcomebias<-DRbias<-matrix(NA, nrep, lgrid)
IPWMSE<-outcomeMSE<-DRMSE<-rep(NA, nrep)
IPWmedianbias<-outcomemedianbias<-DRmedianbias<-matrix(NA, nrep,9)

for (s in 1:nrep){
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
	tempY<-truebeta0fun(individualgrid)+A[i]*truetreatmentfun(individualgrid)+X[i]*truebetafun(individualgrid)+epsilonvalue
	individualhvalue<-tempY+individualgrid
	rval<-approxfun(c(0, individualgrid,1), c(0,individualhvalue, 1),
        method = "constant")
    hvalue<-sapply(grid, rval) 
    Yinverse[i,]<-hvalue
}




identity<-grid


# Fit propensity score model
logisticfit<-glm(A~X, family=binomial(link='logit'))
pihat<-logisticfit$fitted.values

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

estimatemu0DRcomponent1<-estimatemu0DRcomponent2<-matrix(NA, n, lgrid)
for (i in 1:n){
	estimatemu0DRcomponent1[i, ]<-estimatemu0IPWcomponent1[i,]-(A[i]/pihat[i]-1)*estimatemu0m1fit[i,]
		estimatemu0DRcomponent2[i, ]<-estimatemu0IPWcomponent2[i,]-((1-A[i])/(1-pihat[i])-1)*estimatemu0m0fit[i,]
}

estimatemu0DRpsi1<-apply(estimatemu0DRcomponent1, 2, mean)
estimatemu0DRpsi0<-apply(estimatemu0DRcomponent2, 2, mean)

DRmu0inverse<-estimatemu0DRpsi0

### linear interpolation 
IPWmu0inversefunction<-function(u){
	tempindex<-floor(u*(lgrid-1))+1
if (u>0 & u<1){
	IPWmu0inversefunction<-IPWmu0inverse[tempindex]+(u-grid[tempindex])*(IPWmu0inverse[tempindex+1]-IPWmu0inverse[tempindex])*(lgrid-1)
} else if (u<=0) {
 	IPWmu0inversefunction<-0
} else {
	IPWmu0inversefunction<-1
}	
}

### mu0 function
IPWmu0inversefun<-inverse(IPWmu0inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu0
IPWhatmu0<-sapply(grid, IPWmu0inversefun)

### linear interpolation 
outcomemu0inversefunction<-function(u){
	tempindex<-floor(u*(lgrid-1))+1
if (u>0 & u<1){
	outcomemu0inversefunction<-outcomemu0inverse[tempindex]+(u-grid[tempindex])*(outcomemu0inverse[tempindex+1]-outcomemu0inverse[tempindex])*(lgrid-1)
} else if (u<=0) {
 	outcomemu0inversefunction<-0
} else {
	outcomemu0inversefunction<-1
}	
}

### mu0 function
outcomemu0inversefun<-inverse(outcomemu0inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu0
outcomehatmu0<-sapply(grid, outcomemu0inversefun)


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

### mu0 function
DRmu0inversefun<-inverse(DRmu0inversefunction, lower=0, upper=1)

### Estimated Frechet mean hatmu0
DRhatmu0<-sapply(grid, DRmu0inversefun)


# Fit propensity score model
logisticfit<-glm(A~X, family=binomial(link='logit'))
pihat<-logisticfit$fitted.values

IPWLY<-outcomeLY<-DRLY<-IPWcomponent1<-IPWcomponent2<-DRIPWcomponent1<-DRIPWcomponent2<-matrix(NA, n, lgrid)
for (i in 1:n){
Yinversefunction<-function(u){
		tempindex<-floor(u*(lgrid-1))+1
if (u!=1){
	Yinversefunction<-Yinverse[i,tempindex]+(u-grid[tempindex])*(Yinverse[i, tempindex+1]-Yinverse[i, tempindex])*(lgrid-1)
} else {
	Yinversefunction<-1
}	
}

IPWLY[i,]<-sapply(IPWhatmu0, Yinversefunction)
outcomeLY[i,]<-sapply(outcomehatmu0, Yinversefunction)
DRLY[i,]<-sapply(DRhatmu0, Yinversefunction)

IPWcomponent1[i,]<-A[i]*IPWLY[i,]/pihat[i]
IPWcomponent2[i,]<-(1-A[i])*IPWLY[i,]/(1-pihat[i])

DRIPWcomponent1[i,]<-A[i]*DRLY[i,]/pihat[i]
DRIPWcomponent2[i,]<-(1-A[i])*DRLY[i,]/(1-pihat[i])
}

### IPW treatment effect estimate
IPWtxeffect<-apply(IPWcomponent1, 2, mean)-apply(IPWcomponent2, 2, mean)


### regression estimate
outcomelmfit<-lm(outcomeLY~A+X)

### outcome treatment effect estimate
outcometxeffect<-outcomelmfit$coefficients[2,]

DRlmfit<-lm(DRLY~A+X)

newdesign1<-cbind(rep(1,n), rep(1, n), X)
newdesign0<-cbind(rep(1,n), rep(0, n), X)

m1fit<-newdesign1%*%DRlmfit$coefficients
m0fit<-newdesign0%*%DRlmfit$coefficients

DRcomponent1<-DRcomponent2<-matrix(NA, n, lgrid)
for (i in 1:n){
	DRcomponent1[i, ]<-DRIPWcomponent1[i,]-(A[i]/pihat[i]-1)*m1fit[i,]
		DRcomponent2[i, ]<-DRIPWcomponent2[i,]-((1-A[i])/(1-pihat[i])-1)*m0fit[i,]
}

### Doubly robust treatment effect estimate
DRtxeffect<-apply(DRcomponent1, 2, mean)-apply(DRcomponent2, 2, mean)


IPWtransportvalue<-sapply(truemu0, IPWmu0inversefunction)


IPWtxeffectfunction<-function(u){
	tempindex<-floor(u*(lgrid-1))+1
if (u>0 & u<1){
	IPWtxeffectfunction<-IPWtxeffect[tempindex]+(u-grid[tempindex])*(IPWtxeffect[tempindex+1]-IPWtxeffect[tempindex])*(lgrid-1)
} else {
	IPWtxeffectfunction<-0
}
}

IPWparralleleffect<-sapply(IPWtransportvalue, IPWtxeffectfunction)


outcometransportvalue<-sapply(truemu0, outcomemu0inversefunction)


outcometxeffectfunction<-function(u){
	tempindex<-floor(u*(lgrid-1))+1
if (u>0 & u<1){
	outcometxeffectfunction<-outcometxeffect[tempindex]+(u-grid[tempindex])*(outcometxeffect[tempindex+1]-outcometxeffect[tempindex])*(lgrid-1)
} else {
	outcometxeffectfunction<-0
}
}

outcomeparralleleffect<-sapply(outcometransportvalue, outcometxeffectfunction)


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


### The following quantities summarize the bias at each grid point
IPWbias[s,]<-IPWparralleleffect-truetreatmentvalue

outcomebias[s,]<-outcomeparralleleffect-truetreatmentvalue

DRbias[s, ]<-DRparralleleffect-truetreatmentvalue

### The following quantities are MSEs
IPWMSE[s]<-sqrt(mean((IPWbias[s,])^2))
outcomeMSE[s]<-sqrt(mean((outcomebias[s,])^2))
DRMSE[s]<-sqrt(mean((DRbias[s,])^2))

tempsequence<-seq(101, 901, by=100)
### The following quantities summarize the biases from 0.1 to 0.9 quantiles 
IPWmedianbias[s,]<-IPWparralleleffect[tempsequence]-truetreatmentvalue[tempsequence]

outcomemedianbias[s,]<-outcomeparralleleffect[tempsequence]-truetreatmentvalue[tempsequence]

DRmedianbias[s,]<-DRparralleleffect[tempsequence]-truetreatmentvalue[tempsequence]
}

IPWbiasmean<-apply(IPWbias,2, mean)
IPWbiasse<-apply(IPWbias,2, sd)/sqrt(nrep)

outcomebiasmean<-apply(outcomebias[1:100,],2, mean)
outcomebiasse<-apply(outcomebias,2, sd)/sqrt(nrep)

DRbiasmean<-apply(DRbias,2, mean)
DRbiasse<-apply(DRbias,2, sd)/sqrt(nrep)

IPWMSEmean<-mean(IPWMSE)
IPWMSEse<-sd(IPWMSE)/sqrt(nrep)

outcomeMSEmean<-mean(outcomeMSE)
outcomeMSEse<-sd(outcomeMSE)/sqrt(nrep)

DRMSEmean<-mean(DRMSE)
DRMSEse<-sd(DRMSE)/sqrt(nrep)

summarybiasresult<-cbind(IPWbiasmean, IPWbiasse, outcomebiasmean, outcomebiasse, DRbiasmean, DRbiasse)

summaryresult<-c(IPWMSEmean, IPWMSEse, outcomeMSEmean, outcomeMSEse, DRMSEmean, DRMSEse)

summarymedianbiasresult<-rbind(IPWmedianbias, outcomemedianbias, DRmedianbias)

save(summarybiasresult, file=paste("W2hatmu0referencesummarybiasresult_n", n, "_p", p, ".dat",sep=""))

save(summaryresult, file=paste("W2hatmu0referencesummaryresult_n", n, "_p", p, ".dat",sep=""))

save(summarymedianbiasresult, file=paste("W2hatmu0referencesummarymedianbiasresult_n", n, "_p", p, ".dat",sep=""))





























