source('SimulationsFunctions.R')
source('Parameters.R')
source('Costs.R')
source('ProgDyn.R')
source('ProgDynRT.R')
source('FunctionEvalStrategies.R')

tmax=max(Time/delta)

# Needed arguments: 
#nT=size of the grid (without cemetery, nT=nOmega-1), 
#distance=which distance to use for projection on Gamma. Possible choices are "L2" or "Lm"
#RunRT=whether to run dynamic programming with transition matrices computed as products of matrices for r=15 days, 
#allcost=whether to run dynamic programming for different values of visit cost. If FALSE, only CV=1 is run

args <- commandArgs(trailingOnly = TRUE)
nT=args[1]
distance=args[2]
RunRT=args[3]
allcost=args[4]



load("RData/XGrid.RData")
load(paste("RData/ThetaGrid_",nT,"-",distance,".RData",sep=""))

		
CVs=1 # visit cost for simple cost
CVi=1 # visit cost for marker-dependent cost


ProgD()

if(allcost)
{
	CVs=0
	CVi=0
	ProgD()
	CVs=2
	CVi=2
	ProgD()
}

if (RunRT)
{
	RT=array(data=0,dim=c(nT,ngamma,9)) # compute transition matrices of return dates 30 and 60 days as products of transition matrices for 15 days
	for (d in 1:3)
	{
		Rtemp<-rbind(R[,,d],c(rep(0,nT),1))
		RT15<-Rtemp	
		RT30<-Rtemp%*%Rtemp
		RT60<-Rtemp%*%Rtemp%*%Rtemp%*%Rtemp
		RT[,,d]<-RT15[-nrow(RT15),]
		RT[,,3+d]<-RT30[-nrow(RT30),]
		RT[,,6+d]<-RT60[-nrow(RT60),]
	}

	if(allcost)
	{
		CVs=0
		CVi=0
		ProgDRT()
		CVs=2
		CVi=2
		ProgDRT()
	}	
	CVs=1
	CVi=1
	ProgDRT()
}



