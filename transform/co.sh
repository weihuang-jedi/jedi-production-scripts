#!/bin/bash

#rundir=/work2/noaa/gsienkf/weihuang/production/run/sondes/1_rr_observer_ens_solver.run_81.36t9n
 rundir=/work2/noaa/gsienkf/weihuang/production/run/sondes/1_rr_observer_whole_solver.run_81.36t9n

 for var in sondes_tsen sondes_tv sondes_q sondes_uv
 do
   time python concanate-observer.py \
      --run_dir=${rundir} \
      --datestr=2020010112 \
      --nmem=81 \
      --varname=${var} &
 done

 wait

