#=========================================================================
import os
import sys
import getopt
import netCDF4 as nc4
import time
import numpy as np

#-----------------------------------------------------------------------------------------
def comp_rootvar(nc1, nc2):
  v1list = list(nc1.variables)
  v2list = list(nc2.variables)
  n = 0
  for name in v1list:
    print('Var %d: %s' %(n, name))
    n += 1
    v1 = nc1.variables[name][:]
    if(type(v1) == str):
      continue

    if(name not in v2list):
      print('variable name: %s is not in v2list' %(name))
      continue

    v2 = nc2.variables[name][:]
   #dv = v2 - v1

    print('\tv1.max: %f, v1.min: %f' %(np.max(v1), np.min(v1)))
    print('\tv2.max: %f, v2.min: %f' %(np.max(v2), np.min(v2)))
   #print('\tdv.max: %f, dv.min: %f' %(np.max(dv), np.min(dv)))

#-----------------------------------------------------------------------------------------
def comp_var_in_group(ncg1, ncg2):
  v1list = list(ncg1.variables)
  v2list = list(ncg2.variables)
  n = 0
  for name in v1list:
    n += 1
    print('\tVar %d: %s' %(n, name))
    v1 = ncg1.variables[name][:]
   #print('v1.dtype.type = ', v1.dtype.type)
   #if(v1.dtype.type is np.string_):
    if(type(v1[0]) == str):
      continue
    if(name not in v2list):
      print('variable name: %s is not in v2list' %(name))
      continue

    v2 = ncg2.variables[name][:]
    if(type(v2[0]) == str):
      continue
   #if(v2.dtype.type is np.str):
   #  continue
    dv = v2 - v1

    print('\t\tv1.max: %f, v1.min: %f' %(np.max(v1), np.min(v1)))
    print('\t\tv2.max: %f, v2.min: %f' %(np.max(v2), np.min(v2)))
    print('\t\tdv.max: %f, dv.min: %f' %(np.max(dv), np.min(dv)))

#-----------------------------------------------------------------------------------------
def process(f1, f2):
 #grpnamelist = []

  if(os.path.exists(f1)):
    nc1 = nc4.Dataset(f1, 'r')
  else:
    print('f1: %s does not exist.' %(f1))
    sys.exit(-1)

  if(os.path.exists(f2)):
    nc2 = nc4.Dataset(f2, 'r')
  else:
    print('f2: %s does not exist.' %(f2))
    sys.exit(-1)


 #print('list(nc1.variables): ', list(nc1.variables))
 #print('list(nc1.groups): ', list(nc1.groups))

  comp_rootvar(nc1, nc2)

  g1list = list(nc1.groups)
  g2list = list(nc2.groups)

  for name in g1list:
    print('group name: ', name)
    if(name not in g2list):
      print('group name: %s is not in g2list' %(name))
      continue
    ncg1 = nc1.groups[name]
    ncg2 = nc2.groups[name]
    comp_var_in_group(ncg1, ncg2)
   #grpnamelist.append(name)

#-----------------------------------------------------------------------------------------
debug = 1

run_dir = '/work2/noaa/da/weihuang/cycling'
runtype = 'C96_lgetkf_sondesonly'
datestr = '2020010112'
obstype = 'observer'
basename = 'jedi'
casename = 'jedi'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'run_dir=', 'datestr=', 'runtype=',
                                              'obstype=', 'basename=', 'casename='])

for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--run_dir'):
    run_dir = a
  elif o in ('--datestr'):
    datestr = a
  elif o in ('--runtype'):
    runtype = a
  elif o in ('--obstype'):
    obstype = a
  elif o in ('--basename'):
    basename = a
  elif o in ('--casename'):
    casename = a
  else:
    assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
obslist = ['sondes_tsen', 'sondes_tv', 'sondes_uv', 'sondes_q']
basedir = '%s/%s_%s/%s/obsout' %(run_dir, basename, runtype, datestr)
casedir = '%s/%s_%s/%s/observer' %(run_dir, casename, runtype, datestr)
for obs in obslist:
  if(obstype in ['observer', 'solver']):
    f1 = '%s/%s_obs_%s_0000.nc4' %(basedir, obs, datestr)
    f2 = '%s/%s_obs_%s_0000.nc4' %(casedir, obs, datestr)
  else:
    f1 = '%s/%s_obs_%s.nc4' %(basedir, obs, datestr)
    f2 = '%s/%s_obs_%s.nc4' %(casedir, obs, datestr)
  print('f1: ', f1)
  print('f2: ', f2)

  process(f1, f2)

