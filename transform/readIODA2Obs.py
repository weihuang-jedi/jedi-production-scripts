#=========================================================================
#https://unidata.github.io/python-training/workshop/Bonus/netcdf-writing/
import os
import sys
import getopt
import math
import numpy as np
import netCDF4
#import h5py
#from scipy.io import netcdf

#=========================================================================
class ReadIODA2Obs():
  def __init__(self, debug=0, filename=None):
    self.debug = debug

    if(self.debug):
      print('debug = ', debug)
      print('filename = ', filename)

    self.nlocs = 0
    self.ncfile = None

    self.set_filename(filename=filename)

  def __del__(self):
    if(self.ncfile != None):
      self.ncfile.close()

  def set_filename(self, filename=None):
    self.filename = filename
    if(self.ncfile != None):
      self.ncfile.close()
    self.ncfile = netCDF4.Dataset(self.filename, 'a')
   #self.ncfile = netcdf.netcdf_file(self.filename, 'r+')

  def get_groupNvar_name(self, gvstr):
    np = gvstr.rfind('/')
    if (np < 0):
      gname = None
      vname = gvstr
    else:
      gname = gvstr[:np]
      vname = gvstr[np+1:]

   #if(self.debug):
   #  print('gname = ', gname)
   #  print('vname = ', vname)

    return gname, vname

  def get_fileinfo(self, verb=True):
    # NetCDF global attributes
    self.glb_attrs = self.ncfile.ncattrs()
    if self.debug:
        print("NetCDF Global Attributes:")
        for attr in self.glb_attrs:
            print('\t%s:' % attr, repr(self.ncfile.getncattr(attr)))
    self.glb_dims = [dim for dim in self.ncfile.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if self.debug:
        print("NetCDF dimension information:")
    for dim in self.glb_dims:
        print("\tName:", dim)
        print("\t\tsize:", len(self.ncfile.dimensions[dim]))
        if ('nlocs' == dim):
            self.nlocs = len(self.ncfile.dimensions[dim])
        if self.debug:
            print('Dim %s len: %d' %(dim, len(self.ncfile.dimensions[dim])))
    # Variable information.
    self.glb_vars = [var for var in self.ncfile.variables]  # list of nc variables
    if verb:
        print("NetCDF global variable information:")
    for var in self.glb_vars:
        if var not in self.glb_dims:
            print('\tName:', var)
            print("\t\tdimensions:", self.ncfile.variables[var].dimensions)
            print("\t\tsize:", self.ncfile.variables[var].size)
            print(self.glb_attrs(var))

    if self.debug:
        print('\tself.nlocs:', self.nlocs)

    # Group information.
    self.glb_grps = [grp for grp in self.ncfile.groups]  # list of nc variables
    if verb:
        print("NetCDF global group information:")
        for grp in self.glb_grps:
            print('\tGroup Name:', grp)

  def get_group_varlist(self, grpname):
   #varlist = [var for var in self.ncfile.groups[grpname]]
    varlist = [var for var in self.ncfile[grpname]]
    print('grpname: ', grpname)
    print('varlist: ', varlist)
    return varlist

  def get_var_from_group(self, vname, gname):
    ncgroup = self.ncfile[gname]
    var = ncgroup.variables[vname][:]
    return var

  def put_val2groupvar(self, val, vname, gname):
    ncgroup = self.ncfile[gname]
    ncgroup.variables[vname][:] = val

  def get_var(self, ncvarname):
    print('Processing FV3 file %s for variable %s.' % (self.filename, ncvarname))

    gname, vname = self.get_groupNvar_name(ncvarname)

   #print('gname = ', gname)
   #print('vname = ', vname)

    if (gname is None):
      var = self.ncfile.variables[ncvarname][:]
    else:
     #print('gname = ', gname)
     #print('vname = ', vname)

      ncgroup = self.ncfile[gname]
      var = ncgroup.variables[vname][:]
    return var

  def get_latlon(self):
    lat = self.get_var('/MetaData/latitude')
    lon = self.get_var('/MetaData/longitude')

    return lat, lon

#------------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  filename = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.36t1n_36p/obsout/mem001/sondes_tsen_obs_2020011006_0000.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'yaml_file='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--filename'):
      filename = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('filename = ', filename)

  rio = ReadIODA2Obs(debug=debug, filename=filename)
  lat, lon = rio.get_latlon()
 
  print('lat = ', lat)
  print('lon = ', lon)

  rio.get_fileinfo()

