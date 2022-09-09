import getopt
import os, sys

import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import numpy as np

#import cartopy
import cartopy.crs as ccrs
from cartopy import config
from cartopy.util import add_cyclic_point

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0, filename=None):
    self.debug = debug
    self.output = output
    self.filename = filename

  def plot(self):
   #fname = os.path.join(config["repo_data_dir"],
   #                 'netcdf', 'HadISST1_SST_update.nc')

   #dataset = netcdf_dataset(fname)
    dataset = netcdf_dataset(self.filename)

    t = dataset.variables['T_inc'][120, :, :]
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]

    #my preferred way of creating plots (even if it is only one plot)
    ef, ax = plt.subplots(1,1,figsize=(10,5),subplot_kw={'projection': ccrs.PlateCarree()})
    ef.subplots_adjust(hspace=0,wspace=0,top=0.925,left=0.1)

    #get size and extent of axes:
    axpos = ax.get_position()
    pos_x = axpos.x0+axpos.width + 0.01# + 0.25*axpos.width
    pos_y = axpos.y0
    cax_width = 0.04
    cax_height = axpos.height
    #create new axes where the colorbar should go.
    #it should be next to the original axes and have the same height!
    pos_cax = ef.add_axes([pos_x,pos_y,cax_width,cax_height])

    cyclic_data, cyclic_lons = add_cyclic_point(t, coord=lons)
    im = ax.contourf(cyclic_lons, lats, cyclic_data, 60, transform=ccrs.PlateCarree())

    ax.coastlines()

    plt.colorbar(im, cax=pos_cax)

    ax.coastlines(resolution='110m')
    ax.gridlines()
    ax.set_extent([-20, 60, 33, 63])

    #when using this line the positioning of the colorbar is correct,
    #but the image gets distorted.
    #when omitting this line, the positioning of the colorbar is wrong,
    #but the image is well represented (not distorted).
    ax.set_aspect('auto', adjustable=None)

    plt.savefig('t_aspect.png')
    plt.close()

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

