#!/bin/bash
#PBS -q workq
#PBS -N diagpmSMACS0723test0
#PBS -l nodes=comp1:ppn=1
#PBS -V

LOGNAME="SMACSJ0723_test0_analyse"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/SMACS_J0723

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv analyseSMACS.jl
