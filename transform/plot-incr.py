import getopt
import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs
from cartopy import config
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

from netCDF4 import Dataset as netcdf_dataset

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0, filename=None):
    self.debug = debug
    self.output = output
    self.filename = filename

    self.set_default()

  def plot(self):
   #fname = os.path.join(config["repo_data_dir"],
   #                 'netcdf', 'HadISST1_SST_update.nc')

   #dataset = netcdf_dataset(fname)
    dataset = netcdf_dataset(self.filename)

    t = dataset.variables['T_inc'][120, :, :]
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]

    cyclic_data, cyclic_lons = add_cyclic_point(t, coord=lons)

   #ax.coastlines(resolution='110m')
   #ax.gridlines()

    nrows=2
    ncols=1

   #set up the plot
    proj = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw=dict(projection=proj),
                            figsize=(11,8.5))
 
   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      axs[i].set_global()
      cs=axs[i].contourf(cyclic_lons, lats, cyclic_data, transform=proj,
                     cmap=self.cmapname)
     #               cmap=self.cmapname, extend='both')

     #gl = axs[i].gridlines(crs=proj, draw_labels=True,
     #                      linewidth=1, color='green', alpha=0.5, linestyle='.')
     #gl.top_labels = False
     #gl.left_labels = False
     #gl.right_labels = False

     #gl.xlines = False
     #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
     #gl.xformatter = LongitudeFormatter()
     #gl.xlabel_style = {'size': 5, 'color': 'green'}
     #gl.xlabel_style = {'color': 'black', 'weight': 'bold'}

     #gl.ylines = False
     #gl.ylocator = LatitudeLocator()
     #gl.yformatter = LatitudeFormatter()
     #gl.ylabel_style = {'size': 5, 'color': 'green'}
     #gl.ylabel_style = {'color': 'black', 'weight': 'bold'}

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
      axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

      axs[i].set_title(self.runname[i])

   #Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.8,
                        wspace=0.02, hspace=0.02)

   #Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.85, 0.1, 0.9, 0.8])

   #Draw the colorbar
    cbar=fig.colorbar(cs, cax=cbar_ax, orientation='vertical')

   #Add a big title at the top
    plt.suptitle(self.title)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(self.image_name is None):
        image_name = 't_aspect.png'
      else:
        image_name = self.image_name
      plt.savefig(image_name)
      plt.close()
    else:
      plt.show()

  def set_default(self):
    self.image_name = 'sample.png'

    self.runname = ['GSI', 'JEDI', 'GSI', 'JEDI', 'GSI', 'JEDI']

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
   #self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
    self.cmapname = 'jet'

    self.obslat = []
    self.obslon = []

   #self.clevs = np.arange(-0.2, 0.21, 0.01)
   #self.cblevs = np.arange(-0.2, 0.3, 0.1)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'T (C)'
    self.title = 'Temperature Increment'

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_imagename(self, imagename):
    self.image_name = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  filename = 'jedi_increment.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'workdir='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--workdir'):
      workdir = a
    else:
      assert False, 'unhandled option'

  gp = GeneratePlot(debug=debug, output=output, filename=filename)

  gp.plot()

