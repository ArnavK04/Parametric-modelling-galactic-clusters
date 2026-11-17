#!/bin/bash
#PBS -q workq
#PBS -N diagARES_MaxEnt_20261117
#PBS -l nodes=comp1:ppn=32
#PBS -V

LOGNAME="Ares_MaxEnt_20261117_analyse"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Ares/MaxEnt

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 32 analyseARES.jl
