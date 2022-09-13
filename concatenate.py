import os, sys, time

from pathlib import Path
from solo.basic_files import tree
from solo.ioda import Ioda

from multiprocessing import Pool

import getopt

#--------------------------------------------------------------------------------
def concat_files(workdir, obstype, datestring, debugstr):
 #print('workdir: ', workdir)
 #print('obstype: ', obstype)
 #print('datestring: ', datestring)
 #print('debugstr: ', debugstr)

  output_file_path = './%s_%s.nc4' %(obstype, datestring)
  debug = int(debugstr)

  output = Path(output_file_path)

  netcdf_files = []
 
  nc = 0
  search_more = True
  while(search_more):
    filename = '%s/%s_%s_%4.4d.nc4' %(workdir, obstype, datestring, nc)
    if(os.path.exists(filename)):
      netcdf_files.append(filename)
      if(debug):
        print('Processing ... ', filename)
      nc += 1
    else:
      search_more = False

  nc = Ioda('Concatenated From RoundRobin Distribution Obs File')
  nc.concat_files(netcdf_files, str(output.absolute()))
  if(debug):
    print('Finished ', output_file_path)

def multi_run_wrapper(args):
  concat_files(*args)

#--------------------------------------------------------------------------------
if __name__ == '__main__':
  tstart = time.time()

  debug = "1"
  datestring = '2020011006'
  workdir = '/work2/noaa/gsienkf/weihuang/jedi/run/rr_maxpoolsize_tpe/run_80.40t1n_36p/obsout'
  obslist = ['aircraft_q_obs', 'aircraft_tsen_obs', 'aircraft_uv_obs',
             'amsua_n19_obs_m', 'iasi_metop-a_obs_m', 'iasi_metop-b_obs_m',
             'satwind_obs', 'scatwind_obs', 'sfc_ps_obs',
             'sfcship_ps_obs', 'sfcship_q_obs', 'sfcship_tsen_obs',
             'sfcship_tv_obs', 'sfcship_uv_obs',
             'sondes_ps_obs', 'sondes_q_obs', 'sondes_tsen_obs',
             'sondes_tv_obs', 'sondes_uv_obs',
             'vadwind_obs', 'windprof_obs']

 #--------------------------------------------------------------------------------
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'obslist=', 'workdir=', 'datestring='])

  for o, a in opts:
    if o in ('--debug'):
      debug = a
    elif o in ('-slist'):
      obslist = a
    elif o in ('--workdir'):
      workdir = a
    elif o in ('--datestring'):
      datestring = a
    else:
      assert False, 'unhandled option'

 #--------------------------------------------------------------------------------
  optlist = []
  for type in obslist:
    opt = (workdir, type, datestring, debug)
    optlist.append(opt)

 #--------------------------------------------------------------------------------
  print('optlist: ', optlist)
  pool = Pool(8)
  pool.map(multi_run_wrapper, optlist)
  pool.close()

 #--------------------------------------------------------------------------------
  tend = time.time()

  print('concatenate.py done in {:.4f} seconds'.format(tend-tstart))

