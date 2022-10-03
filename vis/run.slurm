#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 02:25:00
#SBATCH -A gsienkf
#SBATCH --partition=bigmem
#SBATCH --job-name=plot_fcst
#SBATCH --output=log.plot_fcst
##SBATCH --mem=0

 source ~/intelenv

 ulimit -S unlimited
 ulimit -c unlimited
#--------------------------------------------------------------------------------------------
 cd /work2/noaa/gsienkf/weihuang/production/run/vis

 time python fcst-jedi-gsi.py --output=1 --type=sfg
 time python fcst-jedi-gsi.py --output=1 --type=sanl

 exit 0

 time python compute_mean_incr.py \
   --workdir=/work2/noaa/gsienkf/weihuang/gsi/gsi_C96_lgetkf_sondesonly \
   --datestr=2020010600 \
   --outfile=gsi_2020010600_mean_incr.nc4

 time python compute_mean_incr.py \
   --workdir=/work2/noaa/gsienkf/weihuang/gsi/jedi_C96_lgetkf_sondesonly \
   --datestr=2020010600 \
   --outfile=jedi_2020010600_mean_incr.nc4

 exit 0

 time python fcst-jedi-gsi.py --output=1 --type=sfg
 time python fcst-jedi-gsi.py --output=1 --type=sanl

