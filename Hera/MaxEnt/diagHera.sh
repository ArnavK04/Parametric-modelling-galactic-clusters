#!/bin/bash
#PBS -q workq
#PBS -N diagHERA_MaxEnt0
#PBS -l nodes=comp1:ppn=8
#PBS -V

LOGNAME="Hera_MaxEnt_20261114OG_analyse"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Hera/MaxEnt

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 8 analyseHERA.jl
