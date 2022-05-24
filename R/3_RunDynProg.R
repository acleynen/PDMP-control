source('SimulationsFunctions.R')
source('Parameters.R')
source('Costs.R')
source('ProgDyn.R')
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

