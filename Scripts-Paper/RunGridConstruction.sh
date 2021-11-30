#!/bin/bash
#SBATCH --job-name=GridConstruction
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
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
Rscript Scripts-Paper/3_RunDynProg.R 642 L2 FALSE FALSE
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 642 25 0.35 2 a L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 2 a L2 100
echo "Running dynamic programming on layer 2"
Rscript Scripts-Paper/3_RunDynProg.R 977 L2 FALSE FALSE
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 977 50 0.35 3 a L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 3 a L2 100
echo "Running dynamic programming on layer 3a"
Rscript Scripts-Paper/3_RunDynProg.R 1260 L2 FALSE FALSE

echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 977 15 0.3 3 b L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 3 b L2 100
echo "Running dynamic programming on layer 3b"
Rscript Scripts-Paper/3_RunDynProg.R 1212 L2 FALSE FALSE
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 1212 50 0.3 4 ba L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 4 ba L2 100
echo "Running dynamic programming on layer 4ba"
Rscript Scripts-Paper/3_RunDynProg.R 1576 L2 FALSE FALSE

echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 1212 20 0.25 4 bb L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 4 bb L2 100
echo "Running dynamic programming on layer 4bb"
Rscript Scripts-Paper/3_RunDynProg.R 1582 L2 FALSE FALSE
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 1582 50 0.25 5 bba L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 5 bba L2 100
echo "Running dynamic programming on layer 5bba"
Rscript Scripts-Paper/3_RunDynProg.R 2100 L2 FALSE FALSE
echo "Adding points with distance L2"
Rscript Scripts-Paper/4_AddPoints.R 1582 20 0.2 5 bbb L2
echo "Computing extended grid  distance L2"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 5 bbb L2 100
echo "Running dynamic programming on layer 5bbb"
Rscript Scripts-Paper/3_RunDynProg.R 2071 L2 FALSE FALSE

########################## Lm ##############################

echo "Computing initial grid Gamma and transition matrix with distance Lm"
Rscript Scripts-Paper/2_ThetaGrid-Init.R Lm 100
echo "Running dynamic programming on initial grid"
Rscript Scripts-Paper/3_RunDynProg.R 642 Lm FALSE FALSE
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 642 15 0.6 2 a Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 2 a Lm 100
echo "Running dynamic programming on layer 2"
Rscript Scripts-Paper/3_RunDynProg.R 919 Lm FALSE FALSE
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 919 50 0.55 3 a Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 3 a Lm 100
echo "Running dynamic programming on layer 3a"
Rscript Scripts-Paper/3_RunDynProg.R 1292 Lm FALSE FALSE


echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 919 20 0.5 3 b Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 3 b Lm 100
echo "Running dynamic programming on layer 3b"
Rscript Scripts-Paper/3_RunDynProg.R 1346 Lm FALSE FALSE
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 1346 50 0.4 4 ba Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 4 ba Lm 100
echo "Running dynamic programming on layer 4ba"
Rscript Scripts-Paper/3_RunDynProg.R 1714 Lm FALSE FALSE
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 1714 125 0.4 5 baa Lm
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 1714 125 0.4 5 baa2 Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 5 baa Lm 100
echo "Running dynamic programming on layer 5baa"
Rscript Scripts-Paper/3_RunDynProg.R 2084 Lm FALSE FALSE


echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 1714 40 0.35 5 bab Lm
echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 1714 40 0.35 5 bab2 Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 5 bab2 Lm 100
echo "Running dynamic programming on layer 5bab"
Rscript Scripts-Paper/3_RunDynProg.R 2126 Lm FALSE FALSE



echo "Adding points with distance Lm"
Rscript Scripts-Paper/4_AddPoints.R 1346 25 0.35 4 bb Lm
echo "Computing extended grid  distance Lm"
Rscript Scripts-Paper/5_ThetaGrid-Added.R 4 bb Lm 100
echo "Running dynamic programming on layer 4bb"
Rscript Scripts-Paper/3_RunDynProg.R 1807 Lm FALSE FALSE



