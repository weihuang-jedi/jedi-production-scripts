#!/bin/bash

 set -x

 python incrment-jedi-gsi.py \
   --output=0 \
   --topdir=/work2/noaa/da/weihuang/cycling \
   --datestr=2020010112 \
   --basename=med.jedi \
   --casename=jedi

 exit 0

