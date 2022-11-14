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

import netCDF4 as nc4

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

  def plot(self, lons, lats, data=[]):
   #ax.coastlines(resolution='110m')
   #ax.gridlines()

    nrows = len(data)
    ncols = 1

   #set up the plot
    proj = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw=dict(projection=proj),
                            figsize=(11,8.5))
 
   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      axs[i].set_global()

      pvar = data[i]

     #print('Plot No. ', i)
     #print('\tpvar.shape = ', pvar.shape)

      vmin = np.min(pvar)
      vmax = np.max(pvar)

     #if((vmax - vmin) > 1.0e-5):
     #  self.clevs, self.cblevs = get_plot_levels(pvar)

      if(i < 2):
       #self.cmapname = 'coolwarm'
       #self.cmapname = 'rainbow'
        self.cmapname = 'jet'
        self.clevs = np.arange(-2.0, 2.1, 0.1)
        self.cblevs = np.arange(-2.0, 3.0, 1.0)
      else:
        self.cmapname = 'bwr'
       #self.clevs = np.arange(-2.0, 2.1, 0.1)
       #self.cblevs = np.arange(-2.0, 3.0, 1.0)
        self.clevs = np.arange(-0.01, 0.011, 0.001)
        self.cblevs = np.arange(-0.01, 0.015, 0.005)

      cyclic_data, cyclic_lons = add_cyclic_point(pvar, coord=lons)

      cs=axs[i].contourf(cyclic_lons, lats, cyclic_data, transform=proj,
                         levels=self.clevs, extend=self.extend,
                         alpha=self.alpha, cmap=self.cmapname)
     #               cmap=self.cmapname, extend='both')

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
      axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

     #axs[i].set_title(self.runname[i])
      title = '%s min: %5.2f, max: %5.2f' %(self.runname[i], vmin, vmax)
      axs[i].set_title(title)

     #minor_ticks_top=np.linspace(0,60,13)
     #minor_ticks_top=np.linspace(-60,0,13)
     #axs[i].set_yticks(minor_ticks_top,minor=True)

     #Draw the colorbar
      cbar=plt.colorbar(cs, ax=axs[i], pad=self.pad,
                        ticks=self.cblevs,
                        orientation='horizontal')

     #cbar.set_label(self.label, rotation=90)
      cbar.set_label(self.label)

   #Add a big title at the top
    plt.suptitle(self.title)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(self.imagename is None):
        imagename = 't_aspect.png'
      else:
        imagename = self.imagename
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

  def set_default(self):
    self.imagename = 'sample.png'

    self.runname = ['First', 'Second', 'Second - First']

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
    self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
   #self.cmapname = 'jet'

    self.clevs = np.arange(-0.2, 0.21, 0.01)
    self.cblevs = np.arange(-0.2, 0.3, 0.1)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Unit (C)'
    self.title = 'Temperature Increment'

  def set_label(self, label='Unit (C)'):
    self.label = label

  def set_title(self, title='Temperature Increment'):
    self.title = title

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_imagename(self, imagename):
    self.imagename = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  topdir = '/work2/noaa/gsienkf/weihuang/production/run/sondes'
  frtdir = '%s/run_80.40t1n_36p' %(topdir)
  snddir = '%s/1_rr_observer_whole_solver.run_81.36t9n' %(topdir)

  frtfile = '%s/observer/sondes_tv_obs_2020010112_0000.nc4' %(frtdir)
  sndfile = '%s/observer/sondes_tv_obs_2020010112_0000.nc4' %(snddir)

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                                                'sndfile=', 'frtfile='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--sndfile'):
      sndfile = a
    elif o in ('--frtfile'):
      frtfile = a
    else:
      assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
  gp = GeneratePlot(debug=debug, output=output)

  ncfrt = nc4.Dataset(frtfile, 'r')
  ncsnd = nc4.Dataset(sndfile, 'r')

  lats = ncsnd.variables['/MetaData/latitude'][:]
  lons = ncsnd.variables['/MetaData/longitude'][:]

#-----------------------------------------------------------------------------------------
  clevs = np.arange(-1.0, 1.01, 0.01)
  cblevs = np.arange(-1.0, 1.1, 0.1)

  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#-----------------------------------------------------------------------------------------
 #varlist = ['tmp', 'ugrd', 'vgrd', 'delp', 'spfh', 'o3mr', 'delz']
 #unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)', 'Unit (Pa)',
 #            'Unit (kg/kg)', 'Unit (ppm)', 'Unit (m)']
 #varlist = ['T_inc', 'delp_inc', 'sphum_inc', 'o3mr_inc', 'delz_inc']
  varlist = ['T', 'delp', 'sphum', 'o3mr', 'DZ', 'ua', 'va']
  unitlist = ['Unit (C)', 'Unit (Pa)',
              'Unit (kg/kg)', 'Unit (kg/kg)', 'Unit (m)', 'Unit (m/s)', 'Unit (m/s)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist)):
    sndvar = ncsnd.variables[varlist[n]][0,:, :, :]
    frtvar = ncfrt.variables[varlist[n]][0,:, :, :]

    print('sndvar.shape = ', sndvar.shape)
   #print('frtvar.shape = ', frtvar.shape)
    nlev, nlat, nlon = sndvar.shape

    gp.set_label(unitlist[n])

    for lev in range(5, nlev, 10):
      v0 = frtvar[lev,:,:]
      v0 = np.where(np.isnan(v0), 0.0, v0)
      v1 = sndvar[lev,:,:]
      v1 = np.where(np.isnan(v1), 0.0, v1)
      v2 = v1 - v0

      data = [v0, v1, v2]

      title = '%s at Level %d' %(varlist[n], lev)
      gp.set_title(title)

      print('Plotting ', title)
      print('\tv0.shape = ', v0.shape)
      print('\tv1.shape = ', v1.shape)
      print('\tv2.shape = ', v2.shape)
 
      print('\tv0.max: %f, v0.min: %f' %(np.max(v0), np.min(v0)))
      print('\tv1.max: %f, v1.min: %f' %(np.max(v1), np.min(v1)))
      print('\tv2.max: %f, v2.min: %f' %(np.max(v2), np.min(v2)))

      imagename = '%s_lev_%3.3d.png' %(varlist[n], lev)
      gp.set_imagename(imagename)

      gp.plot(lons, lats, data=data)

#-----------------------------------------------------------------------------------------
  ncsnd.close()
  ncfrt.close()

