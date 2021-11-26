set.seed(1)
source('SimulationsFunctions.R')
source('Parameters.R')
library(MCMCpack)




################### Frontiers #########################


F1=sort(unique(c(x0,D*exp(-v1not*Time),D*exp(-v1*Time),x0*exp(vprime1*Time),D))) # all boundaries in mode 1
F2=sort(unique(c(x0,D*exp(-v2not*Time),D*exp(-v2*Time),x0*exp(vprime2*Time),D))) # all boundaries in mode 2
F1=F1[F1<=D]
F2=F2[F2<=D]



################# m=1 #################



################# m=1 #################

dF1=diff(F1)
dF2=diff(F2)

epsilon=min(c(0.05,min(dF1)/4,min(dF2)/4))

ind=c()
for (i in 2:(length(F1)-1))
	ind=c(ind,F1[i]-epsilon,F1[i]+epsilon)

ind=sort(ind)

Minimalx1=ind

ind=c(x0,ind,D)
Morex1=c()

for (j in 2:(length(ind)-1))
{
	temp=c(ind[j-1])
	n=1
	while(temp[n]<ind[j])
	{
		if (n>1)
			temp=c(temp,Step(1,1,temp[n],15,0,15,"b")[2]+(floor(j/2)-1)*epsilon) else temp=c(temp,Step(1,1,temp[n],15,0,15,"b")[2])
		n=n+1
	}
	Morex1=c(Morex1,temp[-c(1,length(temp))])
}


x1=sort(c(Minimalx1,Morex1))

################# m=2 #################

ind=c()
for (i in 2:(length(F2)-1))
	ind=c(ind,F2[i]-epsilon,F2[i]+epsilon)

ind=sort(ind)
u=rep(250,length(ind))

Minimalx2=ind
plot(Minimalx2,u,xlim=c(0,D+1),ylim=c(0,H),main="Initial Grid, m=2",xlab="discretization on x",ylab="discretization on u")
abline(v=c(F2),col="red")
abline(v=c(x0,D),col="blue")


ind=c(x0,ind,D)

Morex2=c()
for (j in 2:(length(ind)-1))
{
	temp=c(ind[j-1])
	n=1
	while(temp[n]<ind[j])
	{
		if (n>1)
			temp=c(temp,Step(1,2,temp[n],15,0,15,"a")[2]+(floor(j/2)-1)*epsilon) else temp=c(temp,Step(1,2,temp[n],15,0,15,"a")[2])
		n=n+1
	}
	Morex2=c(Morex2,temp[-c(1,length(temp))])
}



points(Morex2,rep(500,length(Morex2)),col="grey")



x2=sort(c(Minimalx2,Morex2))

################# m=0 #################
timeline=delta*(0:(N/2))
nm0=length(timeline)


############## Complete grid #####################

nomega=nm0+nm1+nm2
Omega<-matrix(ncol=4,nrow=nomega) # 1st column =mode, 2nd=zeta value, 3rd=time since last jump, 4rth=index of the point in the grid

Omega[,1]<-c(rep(0,nm0),rep(1,nm1),rep(2,nm2))
Omega[,2]<-c(rep(x0,nm0),x1,x2)
Omega[,3]<-c(timeline,rep(delta,nm1+nm2))
Omega[,4]<-1:nomega


xval<-Omega[,2]
nOmega=nomega+1

##################### Transition matrix estimation #########################

NMC = 1000 # number of monte carlo simulations
w=0 # time since beginning: set to 0 for grid computations
P=array(data=0,dim=c(nOmega,nOmega,3*length(Time)))
P[nOmega,nOmega,]<-1


start.time <- Sys.time()
for (m in 0:2)
{
	print(m)
	temp=Omega[Omega[,1]==m,] # only points for the current mode

	for (gind in 1:nrow(temp))
	{
		x=temp[gind,2]
		t=temp[gind,3]
		for (rind in 1:length(Time))
		{
			r=Time[rind]
			for (jind in 1:3)
			{
				j=Treatment[jind]
				Sim=sapply(1:NMC,Step,m=m,x=x,t=t,w=w,r=r,j=j) # simulates next process value
				Neighbor=apply(Sim[-nrow(Sim),],2,NNeighbor) # projects on Omega
				Indfill=sort(unique(Neighbor[4,]))
				for (l in Indfill)
				{
					prob=sum(Neighbor[4,]==l)/NMC
					P[temp[gind,4],l,(rind-1)*3+jind]<-prob # adds fraction of observed occurrences to each point
				}
			}
		}
	}
}


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken 



save(Omega,P,nOmega,nm0,nm1,nm2,xval,time.taken,file="RData/XGrid.RData")

