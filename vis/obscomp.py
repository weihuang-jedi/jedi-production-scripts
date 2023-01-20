#=========================================================================
import os
import sys
import types
import getopt

import numpy as np

from genplot import GeneratePlot as genplot
from readIODA2Obs import ReadIODA2Obs

import matplotlib.pyplot as plt

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

#-----------------------------------------------------------------------------------------
  def get_bin_data(self, rawdata):
    vmin = np.min(rawdata)
    vmax = np.max(rawdata)

    bins = 100
    slices = np.linspace(vmin, vmax, bins+1, True).astype(np.float)
    counts = np.diff(slices)

    print('slices = ', slices)
    print('counts = ', counts)

    mean = np.add.reduceat(rawdata, slices[:-1]) / counts
    print('mean=', mean)

    return slices, counts, mean

#-----------------------------------------------------------------------------------------
  def get_histogram(self, data):
    vmin = np.min(data)
    vmax = np.max(data)

    bins = 101
    histogram = np.zeros(bins)
    xp = np.linspace(vmin, vmax, bins, True).astype(np.float)

    vlen = vmax - vmin
    for n in range(len(data)):
      i = int(bins*((data[n]-vmin)/vlen))
      if(i >= bins):
        i -= 1
      histogram[i] += 1

    print('len(xp) = ', len(xp))
    print('len(histogram) = ', len(histogram))
    print('histogram = ', histogram)

    return xp, histogram

#-----------------------------------------------------------------------------------------
  def plot_histogram(self, xp, histogram, title):
    bins = len(histogram)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)

   #ax.bar(bins, histogram)

    ax.plot(xp, histogram, color='r')

    ax.set_title(title)
   #plt.title(title)

    imagename = '%s.png' %(title.replace(' ', '_'))
    plt.savefig(imagename)
    plt.show() 

#------------------------------------------------------------------------------
if __name__ == '__main__':
 #datadir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly/2020010112.no-use-control-member'
 #obstype = 'obsout.1n'
  datadir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly/2020010112'
  obstype = 'observer'
  datestr = '2020010112'
  basename = '%s/%s/sondes_tsen_obs_%s_0000.nc4' %(datadir, obstype, datestr)
  datadir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly/2020010112'
  obstype = 'obsout'
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

#------------------------------------------------------------------------------
  total_obs = len(var)
  nonzero_obs = 0
  npositive = 0
  nnegative = 0
  maxdelt = 1.0e-9
  plat = []
  plon = []
  pvar = []
  for n in range(total_obs):
    if(abs(dv[n]) > maxdelt):
      nonzero_obs += 1
      plat.append(lat[n])
      plon.append(lon[n])
      pvar.append(var[n])
      if(dv[n] > 0.0):
        npositive += 1
      else:
        nnegative += 1

  print('total_obs: %d, non-zero_obs: %d' %(total_obs, nonzero_obs))
  print('positive obs: %d, negative obs: %d' %(npositive, nnegative))
  print('dv min: %f, dv max: %f' %(np.min(dv), np.max(dv)))

  imgname = '%s_%s_%s_nonezero.png' %(grpname, varname, datestr)
  title = '%s %s %s_nonezero' %(grpname, varname, datestr)

  gp.set_imagename(imgname)
  gp.set_title(title)
  gp.obsonly(plat, plon, pvar, title=title)

  title = 'histogram %s %s %s' %(grpname, varname, datestr)
 #slices, counts, mean = c2o.get_bin_data(dv)
 #histogram = float(counts)
 #c2o.plot_histogram(histogram, title)

  xp, histogram = c2o.get_histogram(dv)
  c2o.plot_histogram(xp, histogram, title)

