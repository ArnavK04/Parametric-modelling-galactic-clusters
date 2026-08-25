#!/bin/bash
#PBS -q workq
#PBS -N pm_SMACS0723test0_aNFW
#PBS -l nodes=comp1:ppn=16
#PBS -V

LOGNAME="SMACSJ0723_test0_aNFW"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/SMACS_J0723

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 16 fitSMACS.jl
