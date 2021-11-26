set.seed(100100)
source('SimulationsFunctions.R')
source('Parameters.R')
load("RData/XGrid.RData")
library(dplyr)
library(parallel)
library(MCMCpack)

# Needed arguments:
#numadd=layer of the new grid
#version=variation number of grid at this layer (allows to save different files for different values of thresh)
#distance=which distance to use for projection on Gamma. Possible choices are "L2" or "Lm"
#NMC=number of Monte Carlo simulations for transition matrix estimation

args <- commandArgs(trailingOnly = TRUE)
numadd=args[1]
version=args[2]
distance=args[3]
NMC=as.numeric(args[4]) # number of Monte Carlo simulations 


load(paste("RData/GammaAdded",numadd,version,"-",distance,".RData",sep=""))

################### Grid #########################

Gamma=GammaN
nT=ncol(Gamma)
ngamma=nT+1

if (distance=="Lm") Gdist<-apply(Gamma,2,GammaDist)


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


save(NMC,Gamma,R,Omega,P,nOmega,ngamma,time.taken,nT,nsim,thresh,file=paste("RData/ThetaGrid_",nT,"-",distance,".RData",sep=""))


