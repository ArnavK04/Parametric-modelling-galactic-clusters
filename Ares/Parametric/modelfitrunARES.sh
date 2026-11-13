#!/bin/bash
#PBS -q workq
#PBS -N pm_ARES_test0
#PBS -l nodes=comp1:ppn=32
#PBS -V

LOGNAME="ARES_test0"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Ares/Parametric

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 32 fitARES.jl
