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

     #vmin = np.min(pvar)
     #vmax = np.max(pvar)

     #if((vmax - vmin) > 1.0e-5):
     #  self.clevs, self.cblevs = get_plot_levels(pvar)

      cyclic_data, cyclic_lons = add_cyclic_point(pvar, coord=lons)

      cs=axs[i].contourf(cyclic_lons, lats, cyclic_data, transform=proj,
                         levels=self.clevs, extend=self.extend,
                         alpha=self.alpha, cmap=self.cmapname)
     #               cmap=self.cmapname, extend='both')

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
      axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

      axs[i].set_title(self.runname[i])

   #Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.8,
                        wspace=0.02, hspace=0.02)

   #Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.85])

   #Draw the colorbar
    cbar=fig.colorbar(cs, cax=cbar_ax, pad=self.pad, ticks=self.cblevs,
                      orientation='vertical')

    cbar.set_label(self.label, rotation=90)

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

    self.runname = ['GSI', 'JEDI', 'JEDI - GSI']

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
def get_plot_levels(var):
  vmin = np.min(var)
  vmax = np.max(var)

  print('\tget_plot_lelves: vmin = %f, vmax = %f' %(vmin, vmax))

  vcen = 0.5*(vmax + vmin)

  if(vcen >= 1.0):
    vcen = 10.0*int(vcen/10.0)
    pmax = 10.0*int(0.5 + (vmax-vcen)/10.0) + vcen
    pmin = vcen - 10.0*int(0.5 + (vmax-vcen)/10.0)
  else:
    pmax = vmax
    if(np.abs(vmin) > pmax):
      pmax = np.abs(vmin)
  
    fact = 1.0
    while(pmax < 10.0):
      fact *= 10.0
      pmax *= 10.0
    pmin = fact*vmin
    pcen = 0.5*(pmin + pmax)
    pmax = (10.0*int(0.5 + (pmax-pcen)/10.0) + vcen)/fact
    pmin = (vcen - 10.0*int(0.5 + (pmax-pcen)/10.0))/fact

  print('\tget_plot_lelves: pmin = %f, pcen = %f, pmax = %f' %(pmin, vcen, pmax))

  delt = (pmax - pmin)/200.0
  clevs = np.arange(pmin, pmax + delt, delt)

  delt = (pmax - pmin)/20.0
  cblevs = np.arange(pmin, pmax + delt, delt)

 #print('\tclevs: ', clevs[::10])
 #print('\tcblevs: ', cblevs)

  return clevs, cblevs

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  datestr = '2020010200'
  gsidir = '/work2/noaa/gsienkf/weihuang/gsi/gsi_C96_lgetkf_sondesonly'
  jedidir = '/work2/noaa/gsienkf/weihuang/gsi/jedi_C96_lgetkf_sondesonly'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                                                'jedifile=', 'gsifile='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--jedifile'):
      jedifile = a
    elif o in ('--gsifile'):
      gsifile = a
    else:
      assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
  gp = GeneratePlot(debug=debug, output=output)

  gsifile = './gsi_2020010200_mean_incr.nc4'
  jedifile = './jedi_2020010200_mean_incr.nc4'

  ncgsi = nc4.Dataset(gsifile, 'r')
  ncjedi = nc4.Dataset(jedifile, 'r')

  lats = ncjedi.variables['lat'][:]
  lons = ncjedi.variables['lon'][:]

#-----------------------------------------------------------------------------------------
  clevs = np.arange(-1.0, 1.01, 0.01)
  cblevs = np.arange(-1.0, 1.1, 0.1)

  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#-----------------------------------------------------------------------------------------
  varlist = ['T_inc', 'u_inc', 'v_inc', 'delp_inc', 'sphum_inc', 'o3mr_inc', 'delz_inc']

  unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)', 'Unit (Pa)',
              'Unit (kg/kg)', 'Unit (ppm)', 'Unit (m)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist)):
    jedivar = ncjedi.variables[varlist[n]][:, :, :]
    gsivar = ncgsi.variables[varlist[n]][:, :, :]

    nlev, nlat, nlon = jedivar.shape
   #print('jedivar.shape = ', jedivar.shape)
   #print('gsivar.shape = ', gsivar.shape)

    if('T_inc' == varlist[n]):
      clevs = 2.0*np.arange(-1.0, 1.01, 0.01)
      cblevs = 2.0*np.arange(-1.0, 1.1, 0.1)
    elif('u_inc' == varlist[n]):
      clevs = 10.0*np.arange(-1.0, 1.01, 0.01)
      cblevs = 10.0*np.arange(-1.0, 1.1, 0.1)
    elif('v_inc' == varlist[n]):
      clevs = 10.0*np.arange(-1.0, 1.01, 0.01)
      cblevs = 10.0*np.arange(-1.0, 1.1, 0.1)
    elif('delp_inc' == varlist[n]):
      clevs = np.arange(-1.0, 1.01, 0.01)
      cblevs = np.arange(-1.0, 1.1, 0.1)
    elif('delz_inc' == varlist[n]):
      clevs = np.arange(-1.0, 1.01, 0.01)
      cblevs = np.arange(-1.0, 1.1, 0.1)
    elif('o3mr_inc' == varlist[n]):
      clevs = 1.0-e8*np.arange(-1.0, 1.01, 0.01)
      cblevs = 1.0-e8*np.arange(-1.0, 1.1, 0.1)
    elif('sphum_inc' == varlist[n]):
      clevs = 0.001*np.arange(-1.0, 1.01, 0.01)
      cblevs = 0.001*np.arange(-1.0, 1.1, 0.1)
    
    gp.set_clevs(clevs=clevs)
    gp.set_cblevs(cblevs=cblevs)
    gp.set_label(unitlist[n])

    for lev in range(5, nlev, 50):
      v0 = gsivar[lev,:,:]
      v1 = jedivar[lev,:,:]
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
  ncjedi.close()
  ncgsi.close()

