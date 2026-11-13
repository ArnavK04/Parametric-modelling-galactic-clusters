#!/bin/bash
#PBS -q workq
#PBS -N pm_HERA_test0MaxEntOG
#PBS -l nodes=comp1:ppn=1
#PBS -V

LOGNAME="HERA_test0MaxEntOG"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Hera/MaxEnt

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 1 fitHERA.jl
