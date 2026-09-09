#!/bin/bash
#PBS -q workq
#PBS -N pm_ATsmacsj0723_SIE_t0imgv2
#PBS -l nodes=comp1:ppn=16
#PBS -V

LOGNAME="ATsmacsj0723_SIE_t0_imgv2"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/SMACS_J0723/Atalanta0723

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 16 fitsmacsj0723.jl
