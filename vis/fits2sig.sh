#!/bin/bash

 set -x

 export PROJ_LIB=/work2/noaa/gsienkf/weihuang/anaconda3/share/proj
 export PYTHONPATH=/work/noaa/gsienkf/weihuang/jedi/vis_tools/xESMF/build/lib:/work2/noaa/gsienkf/weihuang/anaconda3/lib
 export LD_LIBRARY_PATH=/work2/noaa/gsienkf/weihuang/anaconda3/lib:${LD_LIBRARY_PATH}
 export PATH=/work2/noaa/gsienkf/weihuang/anaconda3/bin:${PATH}


 cp /work2/noaa/da/weihuang/jedi/case_study/vis/jedi_at_12h_stats.nc jediGSIqc.nc
 cp jedi_at_12h_stats.nc jediFilter.nc

 python plotobfits2_sig.py \
   jediGSIqc jediFilter \
   JediGSIQC-filter-12h hem

 exit 0

