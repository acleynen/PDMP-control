# PDMP-control

This repository contains all code necessary to reproduce the simulation study presented in the manuscript ...

## PDMP parameters

All PDMP parameters (flow, intensity, noise), MPD parameters (horizon, decisions) are encoded in file R/Parameters.R
An interested user might modify these parameters and run grids constructions. If flows should not be exponentials, the simulator (encoded in file R/SimulationsFunctions) should be modified accordingly. So far codes are not implemented to easily allow such modifications.

## Cost parameters

Cost functions and parameters are implemented in file R/Costs.R
Three cost functions are proposed : 
 - the time dependent cost : c(x,d,x')=beta_l*r if l is not the appropriate treatment
 - the marker dependent cost : c(x,d,x')=kappa|zeta'-zeta0|* r +beta_l*r  if process is in mode 0
 - the clinic-like cost : c(x,d,x')=kappa|zeta'-zeta0|* r if zeta' is above a given threshold +beta_l*r if process is in mode 0. 
This last cost is combined with a minimal treatment duration

## Grids construction

To construct first discretisation grid, run R/1_XGrid.R Points will be added in the grids automatically to respect frontiers as described in the paper. This code will create the grid, estimate the transition matrix P and store output in an RData file.

To construct the second discretisation grid, run R/2_ThetaGrid-Init.R This code creates a minimal grid with probability measures charging a point of the first discretisation grid with probability 95% and dispatches the rest of the mass to each other point with a Dirichlet probability. Parameters of this initial grid can easily be modified within the code.

## Dynamic programming

The dynamic programming iterations are encoded in file R/ProgDyn.R. This function can then be launched using R/3_RunDynProg.R using parameters
-cost function
-grid identification through its size
-distance function

## Iterative grid construction

The initial second discretisation grid can be enriched through simulation in the following manner:
 - code R/6_EvalGrids.R evaluates through a given number of simulations (entry parameter) the use of each point of the current grid, and removes points whose density is below a given threshold
 - code R/4_AddPoints.R simulates a given number of trajectories(entry parameter) with current optimal strategy  and adds to the current grid each filters whose projection to the grid is above a given threshold (entry parameter). 
 - code R/5_ThetaGrid-Added.R then computes the transition matrix of the updated grid, and dynamic programming can be run.
 
## Simulations in the paper

The grids were computed by running file Scripts-Paper/RunGridConstruction.sh
All RData files were saved in directory RData.
Then simulations and comparisons were run using Sweave codes Scripts-Paper/GridsEvaluation.Rnw to evaluate the iterative procedure, Scripts-Paper/StrategiesComparisons.Rnw to compare our approach with other strategies, and Scripts-Paper/Filter.Rnw to compare the use of different filter projections.
Finally, the codes to produce the figures of the manuscript can be found in directory figures.

