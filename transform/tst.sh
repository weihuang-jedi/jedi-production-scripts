#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 2
#SBATCH -n 80
#SBATCH -t 00:15:00
#SBATCH -A gsienkf
#SBATCH --partition=orion
##SBATCH --partition=bigmem
#SBATCH --job-name=upgradeobs
#SBATCH --output=log.upgradeobs.o%j
##SBATCH --mem=0

 set -x

 source ~/intelenv
 export HDF5_USE_FILE_LOCKING=FALSE

 cd /work2/noaa/gsienkf/weihuang/production/run/transform

 for var in sondes_q sondes_tsen sondes_tv sondes_uv
 do
   srun -n 80 python upgradeobs.py \
     --run_dir=/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.36t1n_36p \
     --datestr=2020011006 \
     --varname=${var}
 done

 sacct --format=JobID,CPUTime,Elapsed,AveVMSize,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID

