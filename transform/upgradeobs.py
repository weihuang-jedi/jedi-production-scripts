#=========================================================================
import os
import sys
import getopt
import netCDF4 as nc4
import time
import numpy as np

from mpi4py import MPI

#-----------------------------------------------------------------------------------------
def mpi_average(comm, x):
  nprocs = comm.Get_size()
  buf = np.zeros_like(x)
  comm.Allreduce(x, buf, op=MPI.SUM)
  buf /= nprocs
  return buf

#-----------------------------------------------------------------------------------------
def replace_var(comm, infile, outfile, grplist):
  if(os.path.exists(infile)):
    print('infile: ', infile)
  else:
    print('infile: %s does not exist.' %(infile))
    sys.exit(-1)
  if(os.path.exists(outfile)):
    os.remove(outfile)
    print('outfile: ', outfile)
 #else:
 #  print('outfile: %s does not exist.' %(outfile))

  ncin = nc4.Dataset(infile, 'r')
  ncout = nc4.Dataset(outfile, 'w')

 #copy global attributes all at once via dictionary
 #ncout.setncatts(ncin.__dict__)
  ncout.source='JEDI observer only ouptut, each with only one member'
  ncout.comment = 'updated variable hofx_y_mean_xb0 and ombg from all observer files'

 #copy dimensions
  for name, dimension in ncin.dimensions.items():
   #ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
    if dimension.isunlimited():
      ncout.createDimension(name, None)
    else:
      ncout.createDimension(name, len(dimension))

 #copy all var in root group.
  for name, variable in ncin.variables.items():
    ncout.createVariable(name, variable.datatype, variable.dimensions)
    ncout[name][:] = ncin[name][:]
   #copy variable attributes all at once via dictionary
   #ncout[name].setncatts(ncin[name].__dict__)

 #check groups
  for name, group in ncin.groups.items():
    print('name: ', name)
    print('group: ', group)
    ncoutgroup = ncout.createGroup(name)
   #copy all var in group.
    for varname, variable in group.variables.items():
      ncoutgroup.createVariable(varname, variable.datatype, variable.dimensions)
     #copy variable attributes all at once via dictionary
      ncoutgroup[name].setncatts(group[name].__dict__)
      if(name in grplist):
        val = group[varname][:]
        avg = mpi_average(comm, val)
        ncoutgroup[varname][:] = avg[:]
      else:
        ncoutgroup[varname][:] = group[varname][:]

  ncin.close()
  ncout.close()

#-----------------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

debug = 1

run_dir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.36t1n_36p'
datestr = '2020011006'
varname = 'sondes_tsen'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'run_dir=',
                                              'datestr=', 'varname='])
for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--gsi_file'):
    gsi_file = a
  elif o in ('--jedifile'):
    jedifile = a
  elif o in ('--newgsifl'):
    newgsifl = a
  else:
    assert False, 'unhandled option'

obsdir = '%s/observer' %(run_dir)
if(0 == rank):
  if(not os.path.exists(obsdir)):
    os.mkdir(obsdir)

#-----------------------------------------------------------------------------------------
grplist = ['hofx_y_mean_xb0', 'ombg']
infile = '%s/obsout/mem%3.3d/%s_obs_%s_0000.nc4' %(run_dir, rank+1, varname, datestr)
outdir = '%s/mem%3.3d' %(obsdir, rank+1)
if(not os.path.exists(outdir)):
  os.mkdir(outdir)
outfile = '%s/%s_obs_%s_0000.nc4' %(outdir, varname, datestr)

replace_var(comm, infile, outfile, grplist)

