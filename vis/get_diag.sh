#!/bin/bash

 set -x

 edate=2020010606
#edate=2020012218

 plot_stats () {
   argnum=$#
   if [ $argnum -lt 4 ]
   then
     echo "Usage: $0 sdate edate interval, for example: $0 2020010112 2020010812 12 all"
     exit -1
   fi

   sdate=$1
   edate=$2
   interval=$3
   flag=$4
   case1=$5
   case2=$6

  #argnum=0
  #echo "@ gives:"
  #for arg in "$@"
  #do
  #  argnum=$(( argnum + 1 ))
  #  echo "Arg $argnum: <$arg>"
  #done

   python diag.py --sdate=$sdate --edate=$edate \
     --case1=${case1} --case2=${case2} --interval=$interval \
     --datadir=/work2/noaa/da/weihuang/cycling > obs_count_${flag}.csv

   lbl1=$(echo ${case1} | tr '[:lower:]' '[:upper:]')
   lbl2=$(echo ${case2} | tr '[:lower:]' '[:upper:]')

   python plot-jedi-gsi-diag.py --lbl1=${lbl1} --lbl2=${lbl2} \
	--fexp=${case1} --sexp=${case2} \
	--output=1 >> obs_count_${flag}.csv

   dirname=${case2}-${case1}
   rm -rf ${dirname}
   mkdir -p ${dirname}
   mv -f obs_count_${flag}.csv ${dirname}/.
   for fl in diag_omf_rmshumid diag_omf_rmstemp diag_omf_rmswind humidity_rms temp_rms wind_rms
   do
     mv -f ${fl}.png ${dirname}/${fl}_${flag}.png
   done
   for case in ${case1} ${case2}
   do
     mv -f ${case}_stats ${dirname}/${case}_${flag}_stats
     mv -f ${case}_stats.nc ${dirname}/${case}_${flag}_stats.nc
   done
   mv -f *stats* ${dirname}/.
   mv -f *.csv ${dirname}/.
   mv -f *.png ${dirname}/.
 }

 tar cvf ~/jg.tar plot-jedi-gsi-diag.py get_diag.sh
#------------------------------------------------------------------------------
#firstlist=(orig.jedi  orig.gsi new.gsi  orig.gsi  new.gsi   orig.gsi)
#secondlist=(jedi      jedi     jedi     orig.jedi orig.jedi new.gsi)
 firstlist=(orig.jedi  orig.gsi gsi  orig.gsi  gsi   orig.gsi)
 secondlist=(jedi      jedi     jedi orig.jedi orig.jedi gsi)
 for j in ${!firstlist[@]}
 do
   first=${firstlist[$j]}
   second=${secondlist[$j]}
   echo "first: ${first}, second: ${second}"

   plot_stats 2020010118 ${edate} 12 at_6h  ${first} ${second}
   plot_stats 2020010112 ${edate} 12 at_12h ${first} ${second}
   plot_stats 2020010112 ${edate} 6  all    ${first} ${second}

   tar uvf ~/jg.tar ${second}-${first}
 done

 exit 0

