#!/bin/bash
#PBS -q workq
#PBS -N diagHERA_Parametric_20261110
#PBS -l nodes=comp1:ppn=1
#PBS -V

LOGNAME="Hera_Parametric_20261110_analyse"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Hera/Parametric

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv analyseHERA.jl
