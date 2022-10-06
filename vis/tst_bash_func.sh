#!/bin/bash

#set -x

 plot_stats () {
   echo The function location is $0
   echo There are $# arguments
   echo "Argument 1 is $1"
   echo "Argument 2 is $2"
   echo "<$@>" and "<$*>" are the same.
   echo List the elements in a for loop to see the difference!

   argnum=0
   echo "* gives:"
   for arg in "$*"
   do
     argnum=$(( argnum + 1 ))
     echo "Arg $argnum: <$arg>"
   done

   argnum=0
   echo "@ gives:"
   for arg in "$@"
   do
     argnum=$(( argnum + 1 ))
     echo "Arg $argnum: <$arg>"
   done
 }

 plot_stats 2020010118 2020010818 12

 exit 0

diag_omf_rmshumid.png  diag_omf_rmstemp.png  diag_omf_rmswind.png  humidity_rms.png  temp_rms.png  wind_rms.png

 python diag.py --sdate=2020010118 --edate=2020010818 --interval=12
 mv gsi_stats.nc gsi06_stats.nc
 mv gsi_stats gsi06_stats
 mv jedi_stats jedi06_stats
 mv jedi_stats.nc jedi06_stats.nc

 python diag.py --sdate=2020010112 --edate=2020010818 --interval=12
 mv gsi_stats.nc gsi12_stats.nc
 mv gsi_stats gsi12_stats
 mv jedi_stats jedi12_stats
 mv jedi_stats.nc jedi12_stats.nc

 python diag.py --sdate=2020010112 --edate=2020010818 --interval=6
 mv gsi_stats gsi_all_stats
 mv gsi_stats.nc gsi_all_stats.nc
 mv jedi_stats jedi_all_stats
 mv jedi_stats.nc jedi_all_stats.nc

