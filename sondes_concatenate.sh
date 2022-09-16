#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 02:25:00
#SBATCH -A gsienkf
##SBATCH --partition=orion
#SBATCH --partition=bigmem
#SBATCH --job-name=concatenate
#SBATCH --output=log.concatenate
##SBATCH --mem=0

 source ~/intelenv

 ulimit -S unlimited
 ulimit -c unlimited
#--------------------------------------------------------------------------------------------
 cd /work2/noaa/gsienkf/weihuang/production/run

 poolsize=5
 debug=1
 datestring='2020011006'
 workdir='/work2/noaa/gsienkf/weihuang/jedi/per_core_timing/run/soundes/run_80.40t1n_36p/obsout'
 obslist='["sondes_q_obs","sondes_tsen_obs","sondes_uv_obs","sondes_tv_obs"]'

 time python concatenate.py \
   --debug=${debug} \
   --datestring=${datestring} \
   --workdir=${workdir} \
   --obslist=${obslist} \
   --poolsize=${poolsize}

