#!/bin/bash

#set -x
#------------------------------------------------------------------------------
 yyyy=2020
 mm=01
 dd=01
 hh=12

#------------------------------------------------------------------------------
 MYLAYOUT="2,3"
 taskspernode=36
 NUMMEM=81
 nodes=9
 totalcpus=$(( $taskspernode * $nodes ))

 templatedir=/work2/noaa/gsienkf/weihuang/production/run/templates
 topdir=/work2/noaa/gsienkf/weihuang/production/run
 casename=sondes

#------------------------------------------------------------------------------
 datetime=${yyyy}-${mm}-${dd}T${hh}:00:00Z
 yyyymmddhh=${yyyy}${mm}${dd}${hh}

#------------------------------------------------------------------------------
 echo "Case: $casename"
 casedir=${topdir}/${casename}

 workdir=${casedir}/run_${NUMMEM}.${taskspernode}t${nodes}n
 rm -rf ${workdir}
 mkdir -p ${workdir}
 cd ${workdir}

 sed -e "s?TASKSPERNODE?${taskspernode}?g" \
     -e "s?TOTALNODES?${nodes}?g" \
     -e "s?TOTALCPUS?${totalcpus}?g" \
     -e "s?DATATYPE?${casename}?g" \
     -e "s?WORKDIR?${workdir}?g" \
     -e "s?NUMMEM?${NUMMEM}?g" \
     ${templatedir}/slurm.template.each_member > run.slurm

 sed -e "s?LAYOUT?${MYLAYOUT}?g" \
     -e "s?NUMBEROFMEMBERS?${NUMMEM}?g" \
     -e "s?DATETIME?${datetime}?g" \
     ${templatedir}/getkf.yaml.template.1_member.rr.observer > getkf.yaml.observer.template

 sed -e "s?YYYYMMDDHH?${yyyymmddhh}?g" \
     -e "s?MAXPOOLSIZE?1?g" \
     ${templatedir}/${casename}.obs.yaml.template.rr.observer >> getkf.yaml.observer.template

 sed -e "s?LAYOUT?${MYLAYOUT}?g" \
     -e "s?NUMBEROFMEMBERS?${NUMMEM}?g" \
     -e "s?DATETIME?${datetime}?g" \
     ${templatedir}/getkf.yaml.template.1_member.solver > getkf.yaml.solver.template

 sed -e "s?YYYYMMDDHH?${yyyymmddhh}_0000?g" \
     ${templatedir}/${casename}.obs.yaml.template.solver >> getkf.yaml.solver.template

 sbatch run.slurm

