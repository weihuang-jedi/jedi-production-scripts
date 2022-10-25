#!/bin/bash

 workdir=/work2/noaa/gsienkf/weihuang/production/run
 griddir=/work/noaa/global/glopara/fix_nco_gfsv16/fix_fv3_gmted2010/C96

 number_members=80
 n=1
 while [ $n -lt $number_members ]
 do
   cd ${workdir}

      if [ $n -lt 10 ]
      then
         member_str=mem00${n}
      elif [ $n -lt 100 ]
      then
         member_str=mem0${n}
      else
         member_str=mem${n}
      fi

      if [ -d Data/ens/${member_str} ]
      then
        cd Data/ens/${member_str}

        for tile in 1 2 3 4 5 6
        do
          ln -sf ${griddir}/C96_oro_data.tile${tile}.nc oro_data.tile${tile}.nc
          ln -sf ${griddir}/C96_grid.tile${tile}.nc C96_grid.tile${tile}.nc
        done
        ln -sf ${griddir}/grid_spec.nc grid_spec.nc
      fi

      n=$(( $n + 1 ))
 done

