#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 00:15:00
#SBATCH -A gsienkf
##SBATCH --partition=orion
#SBATCH --partition=bigmem
#SBATCH --job-name=interp
#SBATCH --output=log.interp.o%j

 set -x

#module purge
#module load intel/2020.2   impi/2020.2

 source ~/intelenv

#rm -f fv3interp.exe

#mpiicc t.c -o fv3interp.exe

#srun -n 40 fv3interp.exe

#rm -f fv3interp.exe

#mpiifort t.F90 -o fv3interp.exe

#srun -n 40 fv3interp.exe

#rm -f fv3interp.exe

#mpiifort tst.F90 -o fv3interp.exe

#srun -n 40 fv3interp.exe

 time srun -n 40 fv3interp.exe

 sacct --format=JobID,CPUTime,Elapsed,AveVMSize,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID

