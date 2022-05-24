# PDMP-control

This repository contains all code necessary to reproduce the simulation study presented in the manuscript ...

## PDMP parameters

All PDMP parameters (flow, intensity, noise), MPD parameters (horizon, decisions) are encoded in file R/Parameters.R
An interested user might modify these parameters and run grids constructions. If flows should not be exponentials, the simulator (encoded in file R/SimulationsFunctions) should be modified accordingly. So far codes are not implemented to easily allow such modifications.

## Cost parameters

Cost functions and parameters are implemented in file R/Costs.R
Three cost functions are proposed : 
-the time dependent cost : c(x,d,x')=beta_l*r if l is not the appropriate treatment
-the marker dependent cost : c(x,d,x')=kappa|zeta'-zeta0|*r +beta_l*r  if process is in mode 0
-the clinic-like cost : c(x,d,x')=kappa|zeta'-zeta0|*r if zeta' is above a given threshold +beta_l*r if process is in mode 0. 
This last cost is combined with a minimal treatment duration

## Grids construction
