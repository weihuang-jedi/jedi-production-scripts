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
#------------------------------------------------------------------------------
if __name__ == '__main__':
  obstype = 'observer'
  datestr = '2020010112'
  datadir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly/2020010112'
  filename = '%s/%s/sondes_tsen_obs_%s_0000.nc4' %(datadir, obstype, datestr)
  debug = 1
  output = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'filename=',
                                                'datestr=', 'obstype=', 'datadir='])

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
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('filename = ', filename)

#------------------------------------------------------------------------------
  nlon = 360
  nlat = nlon/2 + 1
  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

  gp = genplot(debug=debug, output=output)
  clevs = np.arange(-0.5, 0.51, 0.01)
  cblevs = np.arange(-0.5, 0.6, 0.1)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#------------------------------------------------------------------------------
  rio = ReadIODA2Obs(debug=debug, filename=filename)
 #lat, lon = rio.get_latlon()
 #var = rio.get_var('/GsiHofX/surface_pressure')
  lat, lon, var = rio.get_latlon4var(varname='/hofxm0_10_1/air_temperature')

 #print('lat = ', lat)
 #print('lon = ', lon)
 #print('var = ', var)
  print('len(var) = ', len(var))
  print('var min: %f, var max: %f' %(np.min(var), np.max(var)))

  gp.set_obs_latlon(obslat=lat, obslon=lon)

#------------------------------------------------------------------------------
  gp.set_label('Air Temperature (K)')

  imgname = 'radio-sondes'
  title = 'OBS radio-sondes'

  gp.set_imagename(imgname)
  gp.set_title(title)
  gp.obsonly(lat, lon, var, title=title)

