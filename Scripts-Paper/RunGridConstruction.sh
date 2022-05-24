#!/bin/bash
#SBATCH --job-name=GridConstruction
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 28
#SBATCH -o GridConstruction.out
#SBATCH -e GridConstruction.err


module load R/4.0.2
module load cv-standard
module load gcc/7.5.0

echo "Computing grid Omega and transition matrix"
Rscript Scripts-Paper/1_XGrid.R

########################## L2 ##############################

echo "Computing initial grid Gamma and transition matrix with distance L2"
Rscript Scripts-Paper/2_ThetaGrid-Init.R L2 100
echo "Running dynamic programming on initial grid"
Rscript Scripts-Paper/3_RunDynProg.R 184 L2
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 184	10000 L2	int
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 184 10 0.2 2 a L2	int

echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 2 a L2 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 319 L2
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 319	10000 L2	int
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 319 10 0.2 3 a L2	int

echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 3 a L2 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 549 L2
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 549	10000 L2	int
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 549 100 0.4 4 a L2	int

echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 4 a L2 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 722 L2
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 722	10000 L2	int
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 722 12 0.2 5 b L2	int

echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 5 b L2 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 989 L2
 
 
 ########################## Lm ##############################
 
echo "Computing initial grid Gamma and transition matrix with distance Lm"
Rscript Scripts-Paper/2_ThetaGrid-Init.R Lm 100
echo "Running dynamic programming on initial grid"
Rscript Scripts-Paper/3_RunDynProg.R 184 Lm
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 184	10000 Lm	int
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 184 7 0.2 2 a Lm	int

echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 2 a Lm 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 310 Lm
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 310	10000 Lm	int
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 310 7 0.2 3 a Lm	int

echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 3 a Lm 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 545 Lm
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 545	10000 Lm	int
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 545 40 0.55 4 a Lm	int

echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 4 a Lm 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 709 Lm
echo "Removing unused points"
Rscript Scripts-Paper/6_EvalGrids.R 709	10000 Lm	int
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 709 8 0.2 5 b Lm	int

echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 5 b Lm 100
echo "Running dynamic programming"
Rscript Scripts-Paper/3_RunDynProg.R 1021 Lm

