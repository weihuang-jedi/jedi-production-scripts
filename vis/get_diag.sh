#!/bin/bash

#set -x

 plot_stats () {
   argnum=$#
   if [ $argnum -lt 3 ]
   then
     echo "Usage: $0 sdate edate interval, for example: $0 2020010112 2020010812 12"
     exit -1
   fi

   sdate=$1
   edate=$2
   interval=$3
   flag=$4

  #argnum=0
  #echo "@ gives:"
  #for arg in "$@"
  #do
  #  argnum=$(( argnum + 1 ))
  #  echo "Arg $argnum: <$arg>"
  #done

   python diag.py --sdate=$sdate --edate=$edate --interval=$interval > obs_count_${flag}.csv
   python plot-jedi-gsi-diag.py --output=1 >> obs_count_${flag}.csv
   for fl in diag_omf_rmshumid diag_omf_rmstemp diag_omf_rmswind humidity_rms temp_rms wind_rms
   do
     mv ${fl}.png ${fl}_${flag}.png
   done
   for case in gsi jedi
   do
     mv ${case}_stats ${case}_${flag}_stats
     mv ${case}_stats.nc ${case}_${flag}_stats.nc
   done
 }

 plot_stats 2020010118 2020010818 12 at_6h
 plot_stats 2020010112 2020010812 12 at_12h
 plot_stats 2020010112 2020010812 6  all

 exit 0

