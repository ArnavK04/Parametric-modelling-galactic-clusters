#!/bin/bash
#PBS -q workq
#PBS -N pm_SMACS0723test0_aNFW_BCG_2_35
#PBS -l nodes=comp1:ppn=16
#PBS -V

LOGNAME="SMACSJ0723_test0_aNFW_BCG_2_35"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/SMACS_J0723

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 16 fitSMACS.jl
