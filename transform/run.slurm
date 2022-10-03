#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 02:15:00
#SBATCH -A gsienkf
##SBATCH --partition=orion
#SBATCH --partition=bigmem
#SBATCH --job-name=interp
#SBATCH --output=log.interp.o%j
##SBATCH --mem=0

 set -x

 source ~/intelenv

 cd /work2/noaa/gsienkf/weihuang/production/run/transform

 rm -f core* std.out.* std.err.*

#python plot-compare-incr.py --output=1
#python make-movie.py 

 srun -n 40 python mpi_trans.py

