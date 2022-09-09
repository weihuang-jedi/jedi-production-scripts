#!/bin/bash

 set -x

 source /work2/noaa/gsienkf/weihuang/production/util/intelenv

 module rm python/3.9.2

#ioda-bundle build dir:
 export blddir=/work2/noaa/gsienkf/weihuang/production/build/ioda-bundle
#export blddir=/work2/noaa/gsienkf/weihuang/production/build/fv3-bundle
 export PYTHONPATH=${blddir}/lib/python3.9/pyioda:$PYTHONPATH

#output dir.
 output_dir=/work2/noaa/gsienkf/weihuang/production/run/ioda_v2_data

 datetime=2020010112
#input file fir.
 input_dir=/work2/noaa/gsienkf/weihuang/gsi/C96_lgetkf_sondesonly/${datetime}

 temp_dir=/work2/noaa/gsienkf/weihuang/gsi/C96_lgetkf_sondesonly/${datetime}/diag
 mkdir -p ${temp_dir}
 rm -rf ${temp_dir}/*
 cp ${input_dir}/diag* ${temp_dir}
 
 mkdir -p ${output_dir}
#Convert GSI diag 2 ioda2 format
 python ${blddir}/bin/proc_gsi_ncdiag.py -o ${output_dir} ${temp_dir}

 exit 0

