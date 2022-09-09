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


   #set up the plot
    proj = ccrs.PlateCarree()

    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection=proj))
   #h = ax.pcolormesh(cyclic_lons, lats, cyclic_data, transform=proj, cmap='BuRd')
    h = ax.pcolormesh(cyclic_lons, lats, cyclic_data, transform=proj, cmap=self.cmapname)

    ax.coastlines()
   #ax.gridlines()

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
   #gl.left_labels = False
    gl.right_labels = False
    gl.xlines = False
   #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.xformatter = LongitudeFormatter()
    gl.xlabel_style = {'size': 10, 'color': 'green'}
    gl.xlabel_style = {'color': 'black', 'weight': 'bold'}

   #gl.ylocator = LatitudeLocator()
    gl.yformatter = LatitudeFormatter()
    gl.ylabel_style = {'size': 10, 'color': 'green'}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold'}

   #following https://matplotlib.org/2.0.2/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
   #we need to set axes_class=plt.Axes, else it attempts to create
   #a GeoAxes as colorbar

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

    f.add_axes(ax_cb)
    plt.colorbar(h, cax=ax_cb)

    if(self.output):
      if(self.image_name is None):
        image_name = 't_aspect.png'
      else:
        image_name = self.image_name
      plt.tight_layout()
      plt.savefig(image_name)
      plt.close()
    else:
      plt.show()

  def set_default(self):
    self.image_name = 'sample.png'

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

    self.label = 'Time (sec)'
    self.title = 'Time (sec)'

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

