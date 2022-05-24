source('SimulationsFunctions.R')
source('Parameters.R')
source('Costs.R')
source('ProgDyn.R')
source('FunctionEvalStrategies.R')

tmax=max(Time/delta)

# Needed arguments: 
#nT=size of the grid (without cemetery, nT=nOmega-1), 
#distance=which distance to use for projection on Gamma. Possible choices are "L2" or "Lm"

args <- commandArgs(trailingOnly = TRUE)
nT=args[1]
distance=args[2]



load("RData/XGrid.RData")
load(paste("RData/ThetaGrid_",nT,"-",distance,".RData",sep=""))

		
CVs=1 # visit cost for simple cost
CVi=1 # visit cost for marker-dependent cost


