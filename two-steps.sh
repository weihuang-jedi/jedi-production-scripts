#!/bin/bash

#set -x
#------------------------------------------------------------------------------
 yyyy=2020
 mm=01
 dd=01
 hh=12

#------------------------------------------------------------------------------
 taskspernode=40
 NUMMEM=80

#totalcpus=36
#nodes=1
#MYLAYOUT="2,3"

 totalcpus=288
 nodes=8
 MYLAYOUT="8,6"

#templatedir=/work2/noaa/gsienkf/weihuang/production/util/templates
 templatedir=/work2/noaa/gsienkf/weihuang/production/run/templates
 topdir=/work2/noaa/gsienkf/weihuang/production/run
 casename=sondes

#------------------------------------------------------------------------------
 datetime=${yyyy}-${mm}-${dd}T${hh}:00:00Z
 yyyymmddhh=${yyyy}${mm}${dd}${hh}

#------------------------------------------------------------------------------
 echo "Case: $casename"
 casedir=${topdir}/${casename}

#cd ${casedir}

 workdir=${casedir}/run_${NUMMEM}.${taskspernode}t${nodes}n_${totalcpus}p
 rm -rf ${workdir}
 mkdir -p ${workdir}
 cd ${workdir}

 sed -e "s?TASKSPERNODE?${taskspernode}?g" \
     -e "s?TOTALNODES?${nodes}?g" \
     -e "s?TOTALCPUS?${totalcpus}?g" \
     -e "s?CASENAME?${casename}?g" \
     -e "s?WORKDIR?${workdir}?g" \
     -e "s?NUMMEM?${NUMMEM}?g" \
     -e "s?MYLAYOUT?${MYLAYOUT}?g" \
     ${templatedir}/slurm.template.2steps > run.slurm

 sed -e "s?LAYOUT?${MYLAYOUT}?g" \
     -e "s?NUMBEROFMEMBERS?${NUMMEM}?g" \
     -e "s?DATETIME?${datetime}?g" \
     ${templatedir}/getkf.yaml.template.rr.observer > getkf.yaml.rr.observer

#sed -e "s?YYYYMMDDHH?${yyyymmddhh}?g" \
#    -e "s?MAXPOOLSIZE?${totalcpus}?g" \
#    ${templatedir}/${casename}.obs.yaml.template.rr.observer >> getkf.yaml.rr.observer

 sed -e "s?YYYYMMDDHH?${yyyymmddhh}?g" \
     -e "s?MAXPOOLSIZE?1?g" \
     ${templatedir}/${casename}.obs.yaml.template.rr.observer >> getkf.yaml.rr.observer

 sed -e "s?LAYOUT?${MYLAYOUT}?g" \
     -e "s?NUMBEROFMEMBERS?${NUMMEM}?g" \
     -e "s?DATETIME?${datetime}?g" \
     ${templatedir}/getkf.yaml.template.solver > getkf.yaml.solver

#sed -e "s?YYYYMMDDHH?${yyyymmddhh}?g" \
#    ${templatedir}/${casename}.obs.yaml.template.solver >> getkf.yaml.solver

 sed -e "s?YYYYMMDDHH?${yyyymmddhh}_0000?g" \
     ${templatedir}/${casename}.obs.yaml.template.solver >> getkf.yaml.solver

 obsdir=${workdir}/observer
 sed -e "s?DATESTRING?${yyyymmddhh}?g" \
     -e "s?OBSDIR?${obsdir}?g" \
     -e "s?TOPDIR?${topdir}?g" \
     ${templatedir}/sondes_concatenate.sh.template > concatenate.sh

 chmod +x concatenate.sh

 sbatch run.slurm

