Treatment=c("non","a","b")
delta=15 # smallest interval between two visits
Time=c(delta,2*delta,4*delta)
nd=length(Treatment)*length(Time)

Decisions=matrix(ncol=2,nrow=nd)
Decisions[,1]<-rep(Treatment,length(Time))
Decisions[,2]<-rep(Time,each=3)

H=2400 # Horizon in days
N=H/delta # Iteration horizon in delta increments

nobs=length(obs)

## Parameters of jump time distribution

#Recall mu1(s)=Pi1[(lambda1*s)^alpha1]+(1-Pi1)N(2000,400)

alpha1=1.420543 # shape parameter for mode 1 
alpha2=1.053497 # shape parameter for mode 2
lambda1=1.729762e-05 # scale parameter for mode 1 
lambda2=4.801774e-06 # scale parameter for mode 2 
lambda=c(lambda1,lambda2)
alpha=c(alpha1,alpha2)
Pi1=0.80 # mixture proportion 
Pi2=0.80
Pi=c(Pi1,Pi2)



D=40 # Death frontier 
x0=1 # process value in mode 0
## Careful, in mode 1 use  mu'2 !!! mu'2(s)=(b2*s)^beta2
beta1=-3 # shape parameter for mu'1
beta2=-3# 

b1=20 # scale parameter for mu'1 
b2=20 # 
b=c(b1,b2)
beta=c(beta1,beta2)

v1not=0.02 # slope for uncontrolled disease 1 
v2not=0.006 # 
v1=0.01 # slope for disease 1 with treatment b
v2=0.003 # slope for disease 2 with treatment a
v=c(v1not,v1,v2not,v2)
vprime1=0.077 # slope for disease 1 with treatment a
vprime2=0.025 # slope for disease 2 with treatment b
vprime=c(vprime1,vprime2)


sigma2=1 # variance for Gaussian noise 
flink="id" # link function between process X and observations 


