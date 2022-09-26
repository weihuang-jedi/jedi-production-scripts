#!/bin/bash

 set -x

 python compute_mean_incr.py \
   --dir=/work2/noaa/gsienkf/weihuang/gsi/gsi_C96_lgetkf_sondesonly \
   --datestr=2020010206 \
   --casename=GSI \
   --outfilename=gsi_mean_incr.nc4

 python compute_mean_incr.py \
   --dir=/work2/noaa/gsienkf/weihuang/gsi/jedi_C96_lgetkf_sondesonly \
   --datestr=2020010206 \
   --casename=JEDI \
   --outfilename=jedi_mean_incr.nc4

