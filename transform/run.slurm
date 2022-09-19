#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 02:25:00
#SBATCH -A gsienkf
#SBATCH --partition=orion
##SBATCH --partition=bigmem
#SBATCH --job-name=pool_trans
#SBATCH --output=log.pool_trans
##SBATCH --mem=0

 source ~/intelenv

 ulimit -S unlimited
 ulimit -c unlimited
#--------------------------------------------------------------------------------------------
 cd /work2/noaa/gsienkf/weihuang/production/run/transform

 time python pool_trans.py 

