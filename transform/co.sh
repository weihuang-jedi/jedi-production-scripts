#!/bin/bash

 for var in sondes_tsen sondes_tv sondes_q sondes_uv
 do
   time python concanate-observer.py \
      --run_dir=/work2/noaa/gsienkf/weihuang/production/run/sondes/run_81.36t9n \
      --datestr=2020010112 \
      --nmem=80 \
      --varname=${var}
 done

