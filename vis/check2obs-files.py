#=========================================================================
import os
import sys
import getopt
import netCDF4 as nc4
import time
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

#-----------------------------------------------------------------------------------------
def average(vlist):
  buf = np.zeros_like(vlist[0])
  nmem = len(vlist)
  for n in range(1, nmem):
    buf += vlist[n]
  buf /= (nmem-1)
  return buf

#-----------------------------------------------------------------------------------------
def copy_dimension(ncin, ncout):
 #copy dimensions
  for name, dimension in ncin.dimensions.items():
   #ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
    if dimension.isunlimited():
      ncout.createDimension(name, None)
    else:
      ncout.createDimension(name, len(dimension))

#-----------------------------------------------------------------------------------------
def copy_attributes(ncin, ncout):
 #copy the global attributes to the new file
  inattrs = ncin.ncattrs()
  for attr in inattrs:
    if('_FillValue' != attr):
      ncout.setncattr(attr, ncin.getncattr(attr))

#-----------------------------------------------------------------------------------------
def copy_rootvar(ncin, ncout):
 #copy all var in root group.
  for name, variable in ncin.variables.items():
    ncout.createVariable(name, variable.datatype, variable.dimensions)
   #copy variable attributes all at once via dictionary
    ncout[name].setncatts(ncin[name].__dict__)
    ncout[name][:] = ncin[name][:]

#-----------------------------------------------------------------------------------------
def copy_var_in_group(ncingroup, ncoutgroup):
  fvname = '_FillValue'
 #copy all var in group.
  for varname, variable in ncingroup.variables.items():
    if(fvname in variable.__dict__):
      fill_value = variable.getncattr(fvname)
      newvar = ncoutgroup.createVariable(varname, variable.datatype, variable.dimensions, fill_value=fill_value)
    else:
      newvar = ncoutgroup.createVariable(varname, variable.datatype, variable.dimensions)
    copy_attributes(variable, newvar)
    newvar[:] = ncingroup[varname][:]

#-----------------------------------------------------------------------------------------
def copy_grp2newname(name, n, group, ncout):
  item = name.split('_')
  item[-1] = '%d' %(n)
  newname = '_'.join(item)
  print('No %d name: %s, newname: %s' %(n, name, newname))
  ncoutgroup = ncout.createGroup(newname)
  copy_var_in_group(group, ncoutgroup)

#-----------------------------------------------------------------------------------------
class Compare2Files():
  def __init__(self, debug=0, flnm1=None, flnm2=None):
    self.debug = debug
    self.flnm1 = flnm1
    self.flnm2 = flnm2

    if(os.path.exists(flnm1)):
     #print('flnm1: ', flnm1)
      self.nc1 = nc4.Dataset(flnm1, 'r')
    else:
      self.nc1 = None
      print('flnm1: %s does not exist.' %(flnm1))
      sys.exit(-1)

    if(os.path.exists(flnm2)):
     #print('flnm2: ', flnm2)
      self.nc2 = nc4.Dataset(flnm2, 'r')
    else:
      self.nc2 = None
      print('flnm2: %s does not exist.' %(flnm2))
      sys.exit(-1)

    self.exceptlist = ['EffectiveError0', 'EffectiveQC0',
                       'GsiAdjustObsError', 'GsiEffectiveQC', 'GsiFinalObsError',
                       'GsiHofX', 'GsiHofXBc', 'GsiQCWeight', 'GsiUseFlag',
                       'MetaData', 'ObsBias', 'ObsBias0', 'ObsError',
                       'ObsType', 'ObsValue', 'PreUseFlag', 'PreQC',
                       'hofx_y_mean_xb0']

    self.gp = GeneratePlot(debug, output=0)

   #copy attributes
   #for name in ncinlist[0].ncattrs():
   #  ncout.setncattr(name, ncinlist[0].getncattr(name))

   #copy_dimension(ncinlist[0], ncout)
   #copy_rootvar(ncinlist[0], ncout)

    self.process()
 
  def __del__(self):
    if(self.nc1 is not None):
      self.nc1.close()
    if(self.nc2 is not None):
      self.nc2.close()

#-----------------------------------------------------------------------------------------
  def process(self):
    metadata = self.nc2.groups['MetaData']
    lons = metadata['longitude'][:]
    lats = metadata['latitude'][:]

    ng = 0
   #check groups
    for grpname, grp1 in self.nc1.groups.items():
      if(grpname in self.exceptlist):
        continue
      if(grpname.find('hofx0_') == 0):
        continue
      ng += 1
      print('Grp No %d, name: %s' %(ng, grpname))
      grp2 = self.nc2.groups[grpname]
      nv = 0
      for varname, variable in grp1.variables.items():
        nv += 1
        print('\tVar No %d, name: %s' %(nv, varname))
        val1 = grp1[varname][:]
        val2 = grp2[varname][:]
        vald = val1 - val2
        data = [val1, val2, vald]
        title = '%s %s' %(grpname, varname)
        self.gp.set_title(title)
        self.gp.plot(lons, lats, data)
       #print('\t\ttype(val2[0]) = ', type(val2))
        if(type(val2[0]) == str):
          print('\t\tThis variable is a string')
        else:
          print('\t\tval1 min: %f, max: %f' %(np.min(val1), np.max(val1)))
          print('\t\tvar2 min: %f, max: %f' %(np.min(val2), np.max(val2)))

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
 
   #print('axs = ', axs)
   #print('axs.shape = ', axs.shape)

    axs=axs.flatten()

    i = 0
    for ax in axs.flatten():
      ax.set_global()

      pvar = data[i]

     #print('Plot No. ', i)
     #print('\tpvar.shape = ', pvar.shape)

      vmin = np.min(pvar)
      vmax = np.max(pvar)

      if(i < 2):
        scale = 1
        carr = 'orange'
      else:
        carr = np.arange(len(pvar))
        vlen = vmax - vmin
        if(vlen < 0.01):
          vlen = 0.01
        scale = 10.0*(pvar - vmin)/vlen + 1.0

      cs=ax.scatter(x=lons, y=lats,
            c=carr,
            s=scale,
            alpha=0.5,
            transform=proj, cmap=self.cmapname)

      ax.set_extent([-180, 180, -90, 90], crs=proj)
      ax.coastlines(resolution='auto', color='k')
      ax.gridlines(color='lightgrey', linestyle='-', draw_labels=True)

     #ax.set_title(self.runname[i])
      title = '%s min: %5.2f, max: %5.2f' %(self.runname[i], vmin, vmax)
      ax.set_title(title)

     #Draw the colorbar
      if(i == 2):
       #print('ticks=', self.cblevs)
        cbar=plt.colorbar(cs, ax=ax, ticks=self.cblevs)
       #cbar=plt.colorbar(cs, ax=ax, pad=self.pad,
        #                 ticks=self.cblevs)
        #                 orientation='horizontal')

        tickpos = np.arange(len(self.cblevs))/float(len(self.cblevs))
        cbar.set_ticks(tickpos)
        ticklbl = []
        for v in self.cblevs:
          vl = '%5.2f' %(v)
          ticklbl.append(vl)
        cbar.set_ticklabels(ticklbl)

      i += 1

   #Add a big title at the top
    plt.suptitle(self.title)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(self.imagename is None):
        imagename = 't_scatter.png'
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
   #self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
    self.cmapname = 'jet'

    self.clevs = np.arange(-2.0, 2.1, 0.1)
    self.cblevs = np.arange(-2.0, 2.5, 0.5)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Unit (C)'
    self.title = 'Temperature Increment'

  def set_title(self, title):
    self.title = title

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1

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

  c2f = Compare2Files(flnm1 = frtfile, flnm2 = sndfile)

