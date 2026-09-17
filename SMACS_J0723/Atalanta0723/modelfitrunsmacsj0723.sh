#!/bin/bash
#PBS -q workq
#PBS -N pm_ATsmacsj0723_aNFW_t0imgv2_BCG_pert_vdp_10_250_refsig200_400
#PBS -l nodes=comp1:ppn=32
#PBS -V

LOGNAME="ATsmacsj0723_aNFW_t0_imgv2_BCG_pert_vdp_10_250_refsig200_400"

exec > "/home/arnav/Parametric/logfiles/logfiles/${LOGNAME}.log" 2> "/home/arnav/Parametric/logfiles/errfiles/${LOGNAME}.err"

cd /home/arnav/Parametric/Parametric-modelling-galactic-clusters/SMACS_J0723/Atalanta0723

~/.juliaup/bin/julia --project=/home/arnav/UpdatedLensFactoryEnv -t 32 fitsmacsj0723.jl
