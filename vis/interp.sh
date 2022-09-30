#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 02:25:00
#SBATCH -A gsienkf
#SBATCH --partition=bigmem
#SBATCH --job-name=plot_fcst
#SBATCH --output=log.plot_fcst
##SBATCH --mem=0

 source ~/intelenv

 ulimit -S unlimited
 ulimit -c unlimited
#--------------------------------------------------------------------------------------------
 cd /work2/noaa/gsienkf/weihuang/production/run/vis

#indirname = "/work2/noaa/gsienkf/weihuang/jedi/case_study/allobs_JEDI_full_run/run_80.36t1n_36p/analysis/increment"
#outdirname = "/work2/noaa/gsienkf/weihuang/jedi/case_study/allobs_JEDI_full_run/run_80.36t1n_36p/analysis/increment"

cat > input.nml << EOF
&control_param
 generate_weights = .false.
 output_flnm = "interp2gaussian_grid.nc4"
 wgt_flnm = "/work2/noaa/gsienkf/weihuang/production/run/transform/interp_fv3cube2gaussian/gaussian_weights.nc4"
 indirname = "/work2/noaa/gsienkf/weihuang/jedi/case_study/allobs_RR_develop/run_80.36t1n_36p/analysis/increment"
 outdirname = "/work2/noaa/gsienkf/weihuang/jedi/case_study/allobs_RR_develop/run_80.36t1n_36p/analysis/increment"
 has_prefix = .true.
 prefix = '20200110.030000.'
 use_gaussian_grid = .true.
 gaussian_grid_file = '/work2/noaa/gsienkf/weihuang/production/run/transform/interp_fv3cube2gaussian/gaussian_grid.nc4'
 nlon = 384
 nlat = 192
 nlev = 127
 nilev = 128
 npnt = 4
 total_members = 1
 num_types = 2
 data_types = 'fv_core.res.tile', 'fv_tracer.res.tile',
/
EOF

 /work2/noaa/gsienkf/weihuang/production/run/transform/interp_fv3cube2gaussian/fv3interp.exe

