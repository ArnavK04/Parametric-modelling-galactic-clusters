#!/bin/bash
#PBS -q workq
#PBS -N diagpmsmacsJ0723_aNFW_gold_BCG_pert_vdp10_250_refsig200_400_20261025
#PBS -l nodes=comp1:ppn=1
#PBS -V

LOGNAME="smacsJ0723_aNFW_gold_BCG_pert_vdp10_250_refsig200_400_20261025_analyse"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/SMACS_J0723/Atalanta0723

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv analysesmacsj0723.jl
