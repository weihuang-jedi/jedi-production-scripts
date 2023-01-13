#=========================================================================
import os
import sys
import types
import getopt

import numpy as np

from genplot import GeneratePlot as genplot
from readIODA2Obs import ReadIODA2Obs

import tkinter
import matplotlib
matplotlib.use('TkAgg')

#-----------------------------------------------------------------------------------------
class Compare2Obs():
  def __init__(self, debug=0):
    self.debug = debug

#-----------------------------------------------------------------------------------------
  def cal_sorted_index(self, lat1, lon1, lat2, lon2):
    indx = []
    locn = [i for i in range(len(lat1))]
    maxdelt = 1.0e-6
    for i in range(len(lat1)):
      if(i%100 == 0):
        print('Processing No. ', i)
      for n in range(len(locn)):
        j = locn[n]
        if(np.absolute(lat1[i] - lat2[j]) < maxdelt):
          if(np.absolute(lon1[i] - lon2[j]) < maxdelt):
            indx.append(j)
            locn.pop(n)
            break

    print('len(lat1) = ', len(lat1))
    print('len(indx) = ', len(indx))
  
    return indx

#------------------------------------------------------------------------------
if __name__ == '__main__':
  obstype = 'obsout'
  datestr = '2020010112'
  datadir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly/2020010112'
  basename = '%s/%s/sondes_tsen_obs_%s_0000.nc4' %(datadir, obstype, datestr)
  obstype = 'observer'
  filename = '%s/%s/sondes_tsen_obs_%s_0000.nc4' %(datadir, obstype, datestr)
  debug = 1
  output = 0
  grpname = 'hofxm0_10_1'
  varname = 'air_temperature'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'filename=',
                                                'datestr=', 'obstype=', 'datadir=',
                                                'grpname=', 'varname='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--filename'):
      filename = a
    elif o in ('--datestr'):
      datestr = a
    elif o in ('--obstype'):
      obstype = a
    elif o in ('--datadir'):
      datadir = a
    elif o in ('--grpname'):
      grpname = a
    elif o in ('--varname'):
      varname = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('filename = ', filename)

#------------------------------------------------------------------------------
  grpvarname = '/%s/%s' %(grpname, varname)
  baserio = ReadIODA2Obs(debug=debug, filename=basename)
  baselat, baselon, basevar = baserio.get_latlon4var(varname=grpvarname)

  rio = ReadIODA2Obs(debug=debug, filename=filename)
  lat, lon, var = rio.get_latlon4var(varname=grpvarname)

 #print('lat = ', lat)
 #print('lon = ', lon)
 #print('var = ', var)
  print('len(var) = ', len(var))
  print('var min: %f, var max: %f' %(np.min(var), np.max(var)))

#------------------------------------------------------------------------------
  c2o = Compare2Obs(debug=debug)
  idx = c2o.cal_sorted_index(lat, lon, baselat, baselon)
  dv = var - basevar[idx]

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output)
  gp.set_label('Air Temperature (K)')

  imgname = '%s_%s_%s.png' %(grpname, varname, datestr)
  title = '%s %s %s' %(grpname, varname, datestr)

  gp.set_imagename(imgname)
  gp.set_title(title)
  gp.obsonly(lat, lon, dv, title=title)

