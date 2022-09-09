#!/bin/bash

#set -x

 SYEAR=2020
 EYEAR=2020

 SMONTH=01
 EMONTH=01

 SDAY=01
 EDAY=01

 SHOUR=12
 EHOUR=12

 SYMDH=${SYEAR}${SMONTH}${SDAY}${SHOUR}

 src_dir=/work2/noaa/gsienkf/weihuang/gsi/C96_lgetkf_sondesonly/${SYMDH}
 tar_dir=/work2/noaa/gsienkf/weihuang/production/run/Data/ens

 mkdir -p ${tar_dir}
 cd ${tar_dir}
 rm -f mem* coupler.res

 SMINUTE=0
 EMINUTE=0

 SSECOND=0
 ESECOND=0

 sed -e "s?SYEAR?${SYEAR}?g" \
     -e "s?SMONTH?${SMONTH}?g" \
     -e "s?SDAY?${SDAY}?g" \
     -e "s?SHOUR?${SHOUR}?g" \
     -e "s?SMINUTE?${SMINUTE}?g" \
     -e "s?SSECOND?${SSECOND}?g" \
     -e "s?EYEAR?${EYEAR}?g" \
     -e "s?EMONTH?${EMONTH}?g" \
     -e "s?EDAY?${EDAY}?g" \
     -e "s?EHOUR?${EHOUR}?g" \
     -e "s?EMINUTE?${EMINUTE}?g" \
     -e "s?ESECOND?${ESECOND}?g" \
     ../coupler.res.template > coupler.res

 number_members=80
 n=1
 while [ $n -le $number_members ]
 do
   if [ $n -lt 10 ]
   then
     member_str=mem00${n}
   elif [ $n -lt 100 ]
   then
     member_str=mem0${n}
   else
     member_str=mem${n}
   fi

   ln -sf ${src_dir}/${member_str}/INPUT ${member_str}
   cp coupler.res ${member_str}/.

   n=$(( $n + 1 ))
 done

 cd /work2/noaa/gsienkf/weihuang/production/run/Data

 for dir in crtm \
   fieldmetadata \
   fieldsets \
   fv3files \
   satbias \
   TauCoeff
 do
   if [ ! \( -e "${dir}" \) ]
   then
     ln -sf /work2/noaa/gsienkf/weihuang/jedi/case_study/Data/$dir .
   fi
 done

