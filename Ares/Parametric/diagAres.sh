#!/bin/bash
#PBS -q workq
#PBS -N diagARES_Parametric_20261111
#PBS -l nodes=comp1:ppn=1
#PBS -V

LOGNAME="Ares_Parametric_20261111_analyse"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/Ares/Parametric

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv analyseARES.jl
