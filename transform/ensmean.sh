#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 00:15:00
#SBATCH -A gsienkf
#SBATCH --partition=orion
##SBATCH --partition=bigmem
#SBATCH --job-name=upgradeobs
#SBATCH --output=log.upgradeobs.o%j
##SBATCH --mem=0

#set -x

 module load cdo/1.9.10

 indir=/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/Data/ens
 outdir=${indir}/mem000
 mkdir -p ${outdir}

 if [ ! -f ${outdir}/coupler.res ]
 then
   cp ${indir}/mem001/coupler.res ${outdir}/.
   cp ${indir}/mem001/grid_spec.nc ${outdir}/.
   cp ${indir}/mem001/atm_stoch.res.nc ${outdir}/.
   cp ${indir}/mem001/C96_grid.tile*.nc ${outdir}/.
 fi

 typelist=(fv_core.res fv_srf_wnd.res fv_tracer.res oro_data phy_data sfc_data)

 for i in ${!typelist[@]}
 do
   echo "element $i is ${typelist[$i]}"
   type=${typelist[$i]}
   echo "Working on type: $type"

   tile=0
   while [ $tile -lt 6 ]
   do
     tile=$(( $tile + 1 ))
     echo "\tWorking on tile: $tile"

     ofile=${outdir}/${type}.tile${tile}.nc
     rm -f $ofile

     ifiles=`ls ${indir}/mem*/${type}.tile${tile}.nc`
     cdo ensmean $ifiles $ofile
   done
 done

