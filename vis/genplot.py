#=========================================================================
import os
import sys
import types
import getopt

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import cartopy.crs as ccrs

#from matplotlib.ticker import MultipleLocator
#from modelVerticalpressure import ModelVerticalPressure

import tkinter
import matplotlib
matplotlib.use('TkAgg')

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

   #------------------------------------------------------------------------------
   #filename = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/Data/bkg/fv_core.res.nc'
   #self.mvp = ModelVerticalPressure(debug=debug, filename=filename)
   #prs = self.mvp.get_pressure()
   #print('len(prs) = ', len(prs))
   #for n in range(len(prs)):
   #  print('Level %d pressure %f' %(n, prs[n]))
   #self.logp = self.mvp.get_logp()
   #self.markpres = self.mvp.get_markpres()
   #self.marklogp = self.mvp.get_marklogp()

    self.precision = 1

  def set_precision(self, precision=1):
    self.precision = precision

  def get_markpres(self):
    return self.markpres

  def set_markpres(self, markpres=[]):
    self.markpres = markpres
    self.mvp.set_markpres(markpres = markpres)
    self.marklogp = self.mvp.get_marklogp()

  def set_default(self):
    self.image_name = 'sample.png'

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
    self.cmapname = 'bwr'

    self.clevs = np.arange(-2.0, 2.05, 0.05)
    self.cblevs = np.arange(-2.0, 2.5, 0.5)

    self.obslat = []
    self.obslon = []

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Distance (km)'
    self.title = 'Distance to Patch/Tile Center'

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_obs_latlon(self, obslat=[], obslon=[]):
    self.obslat = obslat
    self.obslon = obslon

  def set_imagename(self, imagename):
    self.image_name = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

  def set_label(self, label):
    self.label = label

  def set_title(self, title):
    self.title = title

  def create_image(self, plt_obj, savename):
    msg = ('Saving image as %s.' % savename)
    print(msg)
    kwargs = {'transparent': True, 'dpi': 500}
    plt_obj.savefig(savename, **kwargs)

  def display(self, output=False, image_name=None):
    if(output):
      if(image_name is None):
        image_name=self.image_name
      self.plt.tight_layout()
      kwargs = {'plt_obj': self.plt, 'savename': image_name}
      self.create_image(**kwargs)
    else:
      self.plt.show()

  def plot(self, x, y, pvar, addmark=0, marker='x', size=3, color='green', title='No Title'):
    X, Y = np.meshgrid(x, y)
    Z = pvar

    print('Plotting ', title)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    proj = ccrs.PlateCarree()

   #cs = plt.contourf(X, Y, z, cmap ="jet")
    cs = ax.contourf(X, Y, Z, levels, colors=colors,
                     transform=proj)
   #cs.cmap.set_under('magenta')
   #cs.cmap.set_over('blue')

   #cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.85)
    cbar = plt.colorbar(cs, ax=ax, orientation='horizontal',
                        pad=.1, fraction=0.06, extend='neither')
   #cblabel =    '1. Thermal High           2. Thermal Low           3. Warm High             '
   #cblabel = '%s 4. Cold Low               5. Warm Low              6. Cold High        ' %(cblabel)
   #cbar.set_label(label=cblabel, weight='bold')
   #cbar.set_label(label=cblabel)
   #cbar.ax.tick_params(labelsize=10)
   #cbar.ax.tick_params(labelsize='xx-small')

    ax.set_extent([-180, 180, -90, 90], crs=proj)
    ax.coastlines(resolution='auto', color='k')
    ax.gridlines(color='lightgrey', linestyle='-', draw_labels=True)
    ax.set_global()
    plt.title(title)
    imagename = '%s.png' %(title.replace(' ', '_'))
    plt.savefig(imagename)
   #plt.show()

   #if(addmark):
   #  self.add_obs_marker(marker=marker, size=size, color=color)

   #self.display(output=self.output, image_name=self.image_name)

  def obsonly(self, obslat, obslon, obsvar, title='No Title'):
    print('Plotting ', title)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    proj = ccrs.PlateCarree()

    colors = [str(item/255.) for item in range(len(obsvar))]

   #plt.scatter(x, y, s=500, c=color)

    vm = np.min(obsvar)
    bm = np.max(obsvar)
    sm = bm - vm
    if(sm < 1.0e-9):
      sm = 1.03-9
    size = np.zeros((len(obsvar)), dtype=float)
    for n in range(len(size)):
      size[n] = 1.0 + (obsvar[n]-vm)/sm
    obsplot = ax.scatter(obslon, obslat, s=size, c=obsvar,
                         cmap=self.cmapname, 
                         alpha=self.alpha)
   #obsplot = ax.scatter(obslon, obslat, s=100, c=colors,
   #                     alpha=self.alpha)

   #cb = plt.colorbar(orientation=self.orientation,
   #                  pad=self.pad, ticks=self.cblevs)

   #cb.set_label(label=self.label, size=self.size, weight=self.weight)

   #cb.ax.tick_params(labelsize=self.labelsize)
   #if(self.precision == 0):
   #  cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in self.cblevs], minor=False)
   #elif(self.precision == 1):
   #  cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in self.cblevs], minor=False)
   #elif(self.precision == 2):
   #  cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in self.cblevs], minor=False)
   #else:
   #  cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in self.cblevs], minor=False)

    ax.set_title(self.title)

    ax.set_extent([-180, 180, -90, 90], crs=proj)
    ax.coastlines(resolution='auto', color='k')
    ax.gridlines(color='lightgrey', linestyle='-', draw_labels=True)
    ax.set_global()
    plt.title(title)

    imagename = '%s.png' %(title.replace(' ', '_'))
    plt.savefig(imagename)
    plt.show()

# ----
if __name__ == '__main__':
  debug = 1
  output = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

 #gp = GeneratePlot(debug=debug, output=output, lat=lat, lon=lon)
 #gp.plot(pvar)

