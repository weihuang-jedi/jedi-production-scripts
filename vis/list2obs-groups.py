#=========================================================================
import os
import sys
import getopt
import netCDF4 as nc4
import time
import numpy as np

#-----------------------------------------------------------------------------------------
class Compare2Obs():
  def __init__(self, debug=0):
    self.debug = debug

  def list_rootvar(self, nc1, nc2):
    v1list = list(nc1.variables)
    v2list = list(nc2.variables)
    n = 0
    for name in v1list:
      print('Var %d: %s' %(n, name))

#-----------------------------------------------------------------------------------------
  def list_var_in_group(self, ncg1, ncg2):
    v1list = list(ncg1.variables)
    v2list = list(ncg2.variables)
    maxdelt = 1.0e-6
    n = 0
    for name in v1list:
      n += 1
      print('\tvariable No %d: %s' %(n, name))
      if(name not in v2list):
        print('variable name: %s is not in v2list' %(name))
        continue

#-----------------------------------------------------------------------------------------
  def process(self, f1, f2):
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
  
    self.list_rootvar(nc1, nc2)
  
    g1list = list(nc1.groups)
    g2list = list(nc2.groups)
  
    n = 0
    for name in g1list:
      n += 1
      print('group %d: %s' %(n, name))
      if(name not in g2list):
        print('group name: %s is not in g2list' %(name))
        continue
      ncg1 = nc1.groups[name]
      ncg2 = nc2.groups[name]
      self.list_var_in_group(ncg1, ncg2)

    print('len(g1list) = ', len(g1list))
    print('len(g2list) = ', len(g2list))
  
#-----------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1

  run_dir = '/work2/noaa/da/weihuang/cycling'
  runtype = 'C96_lgetkf_sondesonly'
  datestr = '2020010618'
  obstype = 'observer'
 #basename = 'separate-obs.jedi'
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

  c2o = Compare2Obs(debug=debug)

  #-----------------------------------------------------------------------------------------
  obslist = ['sondes_tsen', 'sondes_tv', 'sondes_uv', 'sondes_q']
 #basedir = '%s/%s_%s/%s/solver' %(run_dir, basename, runtype, datestr)
  basedir = '%s/%s_%s/%s/postobserver' %(run_dir, basename, runtype, datestr)
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
  
    c2o.process(f1, f2)

