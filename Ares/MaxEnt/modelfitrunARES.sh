#!/bin/bash
#PBS -q workq
#PBS -N pm_ARES_test0MaxEnt
#PBS -l nodes=comp1:ppn=1
#PBS -V

LOGNAME="ARES_test0MaxEnt"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Ares/MaxEnt

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 1 fitARES.jl
