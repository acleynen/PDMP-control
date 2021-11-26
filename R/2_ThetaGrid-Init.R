set.seed(100100)
source('SimulationsFunctions.R')
source('Parameters.R')
load("RData/XGrid.RData")
library(dplyr)
library(parallel)
library(MCMCpack)



# Needed arguments:
#distance=which distance to use for projection on Gamma. Possible choices are "L2" or "Lm"
#NMC=number of Monte Carlo simulations for transition matrix estimation

args <- commandArgs(trailingOnly = TRUE)
distance=args[1]
NMC=as.numeric(args[2]) # number of Monte Carlo simulations 


################### Construction of the grid #########################


alpha0=rep(0.1,nm0)
alpha0m=rep(1,nOmega-1-nm0)
alpha1=rep(0.1,nm1)
alpha1m=rep(1,nOmega-1-nm1)
alpha2=rep(0.1,nm2)
alpha2m=rep(1,nOmega-1-nm2)


#Constructs a filter with 95% probability on an element of Omega, and randomly dispatches the rest of the mass on all other elements of Omega through flat Dirichlet distribution
alpha=rep(1,nOmega-2)
dirtemp<-rdirichlet(nOmega-1,alpha)
GammaT95<-matrix(0.95,ncol=nOmega-1,nrow=nOmega-1)
for (i in 1:(nOmega-1))
	GammaT95[,i][-i]<-dirtemp[i,]*0.05
	

#Constructs a filter with 80% (resp 50%) probability on an element of Omega, and dispatches the rest of the mass uniformly on all elements of Omega whose zeta value are closer than 3 sigmas to the current element, plus the first element of mode 0 (this ensures the probability of E_0 is non-negative)
sig3=3*sqrt(sigma2)
GammaT80<-matrix(0,ncol=nm1+nm2,nrow=nOmega-1)
GammaT50<-matrix(0,ncol=nm1+nm2,nrow=nOmega-1)
for (i in 1:(nm1+nm2))
{
	realindex=nm0+i
	x12=xval[(nm0+1):(nOmega-1)]
	xtemp=xval[realindex]
	indtemp=c(1) 
	indextemp=which((x12<(xtemp+sig3)) & (x12>(xtemp-sig3)))+nm0
	indtemp=c(indtemp,indextemp[!is.element(indextemp,realindex)])
	GammaT80[,i][(nm0+i)]<-0.8
	GammaT80[,i][indtemp]<-0.2/length(indtemp)
	GammaT50[,i][(nm0+i)]<-0.5
	GammaT50[,i][indtemp]<-0.5/length(indtemp)	
}	

	
#Constructs a filter with 80% probability divided on pairs of elements of Omega, and randomly dispatches the rest of the mass on all other elements of Omega through flat Dirichlet distribution. Pairs are either closest elements of the same mode, or elements with close zeta values of different modes
alpha=rep(1,nOmega-2)
#40+40%
Pairs=c()
for (i in 1:(nm0-1)) # Pairs of closest elements within a same mode 
Pairs=cbind(Pairs,c(i,i+1))
for (i in (nm0+1):(nm0+nm1-1))
Pairs=cbind(Pairs,c(i,i+1))
for (i in (nm0+nm1+1):(nOmega-2))
Pairs=cbind(Pairs,c(i,i+1))
x1=xval[(nm0+1):(nm0+nm1)]
x2=xval[(nm0+nm1+1):(nm0+nm1+nm2)]
for (i in 1:10) # Pairs of elements of close zeta values with different modes
{
	a=which.min(abs(x1[i]-x2))
	if (a>1)
		a=c(a-1,a)
	if (max(a)<nm1)
		a=c(a,max(a)+1)	
	for (j in 1:length(a))	
		Pairs=cbind(Pairs, c(nm0+i,nm0+nm1+a[j]))
}

alphaa=rep(1,nOmega-3)
dirtemp2<-rdirichlet(ncol(Pairs),alphaa)
GammaT4040<-matrix(0.4,ncol=ncol(Pairs),nrow=nOmega-1)
for (i in 1:ncol(Pairs))
	GammaT4040[,i][-Pairs[,i]]<-dirtemp2[i,]*0.2
		
# Final grid
Gamma=cbind(c(1,rep(0,nOmega-2)),GammaT95,GammaT80,GammaT50,GammaT4040)
nT=ncol(Gamma)
ngamma=nT+1









##################### Transition matrix estimation #########################

start.time <- Sys.time()


Restimate<-function(g)
{
	Rtemp<-matrix(0,nrow=ngamma,ncol=nd)
	Epsilon<-rnorm(NMC,0,sd=sqrt(sigma2))	
	PT=Ptheta(Gamma[,g])
	for (e in 1:NMC)
	{
		NoiseM=NoiseMat(Epsilon[e],flink)
		Psi=EvalPsi(NoiseM,PT) # filter matrix: one filter in each row, associated to a noise associated to an element of Omega 
		for (dind in 1:nd)	
		{
			if (distance=="L2")
			{
				NT=apply(Psi[,,dind],1,NNeighborT)[nOmega,] # For each row, the nearest neighbour in Gamma
			} else
			{
				NT=apply(Psi[,,dind],1,NNeighborBlockT)[nOmega,] 
			}
				A=data.frame(Index=NT,Prob=PT[,dind]) # first column: index of closest neighbour in Gamma, 2nd column, probability to arrive at this point 
				B=data.frame(A %>% group_by(Index) %>%  summarise(Proba = sum(Prob))) # Concatenates (and orders) the probabilities for each possible neighbour 
				Rtemp[as.numeric(B[,1]),dind]<-Rtemp[as.numeric(B[,1]),dind]+as.numeric(B[,2]) # updates the transition matrix 
		}
	}
	Rtemp<-Rtemp/NMC
	for (dind in 1:nd)
			Rtemp[ngamma,dind]<-Rtemp[ngamma,dind]+P[1:(nOmega-1),nOmega,dind]%*%Gamma[,g]
	return(Rtemp)		
}

L<-mclapply(1:nT,Restimate,mc.cores=28)


R=array(data=0,dim=c(nT,ngamma,nd))

for (g in 1:nT)
	R[g,,]<-L[[g]]

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken 


save(Gamma,R,Omega,P,nOmega,ngamma,time.taken,nT,NMC,file=paste("RData/ThetaGrid_",nT,"-",distance,".RData",sep=""))


