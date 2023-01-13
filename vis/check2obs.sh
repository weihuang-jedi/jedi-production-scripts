#!/bin/bash

 set -x

 datestr=2020010412

 python compare2obs.py \
   --run_dir=/work2/noaa/da/weihuang/cycling \
   --datestr=${datestr} \
   --obstype=observer \
   --basename=med.jedi \
   --casename=jedi | tee obs.diff.${datestr}

 exit 0

