RNGkind("L'Ecuyer-CMRG")
source('SimulationsFunctions.R')
source('Parameters.R')
load("RData/XGrid.RData")
source('Costs.R')
source('FunctionEvalStrategies.R')
library(parallel)

NbClusters=28

# Needed arguments:
#nT=size of the grid (without cemetery, nT=nOmega-1), 
#nsim=number of simulations to perform from which to compute weights

#distance=which distance to use for projection on Gamma. Possible choices are "L2" or "Lm"


args <- commandArgs(trailingOnly = TRUE)
nT=args[1]
nsimadd=args[2]
distance=args[3]
cost=args[4]


load(paste("RData/ThetaGrid_",nT,"-",distance,".RData",sep=""))


CVs=1; CVi=1
step="all"
load(paste("RData/ProgDyn_",nT,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))			
if (distance=="Lm") Gdist<-apply(Gamma,2,GammaDist);	


SimToExplore<-function(sim)
{

	Ind<-c()
	xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N)) # initialisation of the Process
	yn=0
	Psin=as.matrix(c(1,rep(0,nOmega-2))) # first filter is always dirac in (0,x0,0,0)
	Strat=list(dec=rep("non",N),vis=rep(0,N)) 
	RealTime=0 # Since dynamic programming is run without track of time, counter of time spent is updated separately 
	n=1
	if (cost=="simple")	decision<-Decsimple[1,1] else decision<-Decint[1,1]# decision at first time point: dynamic programming output at time 0 with filter dirac in (0,x0,0,0)
	Strat$dec[1]<-Decisions[decision,1]
	Strat$vis[1]<-Decisions[decision,2]				
	while(RealTime<H)
	{
		xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n]) # process at iteration n+1 when starting from iteration n with decision Strat[n]
		xn$traj[n+1]=xnew[2]
		xn$time[n+1]=xnew[3]
		xn$mode[n+1]=xnew[1]
		xn$realtime[n+1]=xnew[4]
		RealTime<-RealTime+as.numeric(Strat$vis[n])

		if(xnew[1]==100) # xnew=100 if patient dies. In which case time is artificially set to horizon H to send process to cemetery state
		{
			ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
		} else
		{	
			vmult=which(Time==Strat$vis[n])
			dmult=which(Treatment==Strat$dec[n])
			dind=(vmult-1)*3+dmult # index of decision when combining next visit date and treatment
			ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0) # observation is always positive
			# next five lines compute new filter value	
			PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
			eG<-evalGauss(ynew-xval)
			Num=PTemp*eG
			Den=eG%*%PTemp
			Psitemp=Num/as.vector(Den)
			Psin=cbind(Psin,Psitemp)
			
			if (distance=="L2")
			{
				NT=NNeighborT(Psitemp)[nOmega] # closest neighbour in Gamma with L2 distance
				} else
			{
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				}
			
			indGamma=NT
			indvisit=RealTime/delta
			if (RealTime<H)
			{
				if (cost=="simple") decision<-Decsimple[indvisit+1,indGamma] else decision<-Decint[indvisit+1,indGamma]
				Strat$dec[n+1]<-Decisions[decision,1]
				Strat$vis[n+1]<-Decisions[decision,2]
			} else 
			{
				Strat$dec[n+1]<-nd+1	
			}	
			yn=c(yn,ynew)
		}
		Ind=c(Ind,indGamma)
		n<-n+1
	}
	return(Ind)
}


L<-mclapply(1:nsimadd,SimToExplore,mc.cores=NbClusters,mc.set.seed = TRUE)

AllInds=c()
for (i in 1:nsimadd)
	if (!is.null(L[[i]]))
		AllInds<-c(AllInds,L[[i]])

T=table(AllInds)
n1=as.numeric(rownames(T))
n2=(1:nT)[!is.element(1:nT,n1)]

all=c(n1,n2)

w=c(T/sum(T),rep(0,length(n2)))

a=sort(all,index.return=T)

word=w[a$ix]

Gdist<-apply(Gamma,2,GammaDist)
Gdistw=rbind(Gdist,word)

thresh=0.001
rm=which(word<thresh/nT)

if (is.element(1,rm))
	rm=rm[-1]

Gamma=Gamma[,-rm]

save(rm,Gamma,Gdistw,file=paste("RData/EvalG",nT,"-",distance,".RData",sep=""))


