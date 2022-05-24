Treatment=c("non","a","b")
delta=15 # smallest interval between two visits
Time=c(delta,2*delta,4*delta)
nd=length(Treatment)*length(Time)

Decisions=matrix(ncol=2,nrow=nd)
Decisions[,1]<-rep(Treatment,length(Time))
Decisions[,2]<-rep(Time,each=3)

H=2400
N=H/delta # Horizon in days
obs=delta*(0:floor(N/delta)) # Iteration horizon in delta increments
nobs=length(obs)

## Parameters of jump time distribution

#intensity in mode 0
tau1=c(750,500)
tau2=5*365
tau3=6*365
nu1=c(-2*log(0.8)/tau1[1],-2*log(0.8)/tau1[2])
nu2=c(-2/(2*H-tau2-tau3)*(nu1[1]/2*(tau3+tau2-tau1[1])+log(0.1)),-2/(2*H-tau2-tau3)*(nu1[2]/2*(tau3+tau2-tau1[2])+log(0.1)))

# intensity in modes 1 and 2 
beta1=-0.8 # shape parameters
beta2=-0.8# 

b1=1000 # scale parameters
b2=1000 # 

b=c(b1,b2)
beta=c(beta1,beta2)

## Flow parameters

D=40 # Death frontier 
x0=1 # process value in mode 0

v1not=0.02 # slope for uncontrolled disease 1
v2not=0.006 # slope for uncontrolled disease 2
v1=0.01 # slope for disease 1 with treatment b
v2=0.003 # slope for disease 2 with treatment a
v=c(v1not,v1,v2not,v2)
vprime1=0.077 # slope for disease 1 with treatment a
vprime2=0.025 # slope for disease 2 with treatment b
vprime=c(vprime1,vprime2)


## Observation parameters

sigma2=1 # variance for Gaussian noise 
flink="id" # link function between process X and observations 


