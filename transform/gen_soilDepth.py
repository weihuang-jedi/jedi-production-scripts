#=========================================================================
import os
import sys
import getopt
import netCDF4 as nc4
import time
import numpy as np

#-----------------------------------------------------------------------------------------
def gen_dimension(ncout):
  ncout.createDimension('Time', None)
  ncout.createDimension('xaxis_1', 3)
  ncout.createDimension('xaxis_2', 4)
  ncout.createDimension('xaxis_3', 7)

#-----------------------------------------------------------------------------------------
def gen_var(ncout):
  nx1 = 3
  nx2 = 4
  nx3 = 7

  attrdict = {}
  attrdict['long_name'] = 'Time'
  attrdict['units'] = 'time level'
  attrdict['cartesian_axis'] = 'T'
  time = ncout.createVariable('Time', float, ('Time',))
  time.setncatts(attrdict)
  time[0:1] = 0.0

  attrdict['long_name'] = 'xaxis_1'
  attrdict['units'] = 'meters'
  attrdict['cartesian_axis'] = 'X'
  xaxis_1 = ncout.createVariable('xaxis_1', float, ('xaxis_1',))
  xaxis_1.setncatts(attrdict)
  xaxis_1[:] = np.arange(nx1)

  attrdict['long_name'] = 'xaxis_2'
  xaxis_2 = ncout.createVariable('xaxis_2', float, ('xaxis_2',))
  xaxis_2.setncatts(attrdict)
  xaxis_2[:] = np.arange(nx2)

  attrdict['long_name'] = 'xaxis_3'
  xaxis_3 = ncout.createVariable('xaxis_3', float, ('xaxis_3',))
  xaxis_3.setncatts(attrdict)
  xaxis_3[:] = np.arange(nx3)

  attrdict = {}
  attrdict['long_name'] = 'snowDepth'
  attrdict['units'] = 'meters'
  snowDepth = ncout.createVariable('snowDepth', float, ('Time', 'xaxis_1'))
  snowDepth.setncatts(attrdict)
  snowDepth[0, :] = [0.1, 0.2, 0.5]

  attrdict = {}
  attrdict['long_name'] = 'soilDepth'
  attrdict['units'] = 'meters'
  soilDepth = ncout.createVariable('soilDepth', float, ('Time', 'xaxis_2'))
  soilDepth.setncatts(attrdict)
  soilDepth[0, :] = [0.1, 0.5, 1.0, 2.0]

  attrdict = {}
  attrdict['long_name'] = 'snowDepthPlusSoilDepth'
  attrdict['units'] = 'meters'
  snowNsoilDepth = ncout.createVariable('snowNsoilDepth', float, ('Time', 'xaxis_3'))
  snowNsoilDepth.setncatts(attrdict)
  snowNsoilDepth[0, :] = [0.1, 0.2, 0.5, 0.1, 0.5, 1.0, 2.0]

#-----------------------------------------------------------------------------------------
def genFile(outfile):
  if(os.path.exists(outfile)):
    os.remove(outfile)
    print('outfile: ', outfile)

  ncout = nc4.Dataset(outfile, 'w')

  ncout.comment = 'Manually created soil depth, and snow depth'

  gen_dimension(ncout)
  gen_var(ncout)
 
  ncout.close()

#-----------------------------------------------------------------------------------------
outfile = 'land.nc4'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'outfile='])
for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--infile'):
    infile = a
  elif o in ('--outfile'):
    outfile = a
  else:
    assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------

genFile(outfile)

