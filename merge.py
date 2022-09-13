#=========================================================================
import os, sys
import getopt
import netCDF4 as nc4

#=========================================================================
class MergeNetcdfFiles():
  def __init__(self, debug=0, obsdir='unknown', filetype='unknown',
                     outfile='combined_files.nc4'):
    self.debug = debug
    self.obsdir = obsdir
    self.outfile = outfile
    self.filetype = filetype

    self.filelist = []
    self.ncarray = []

    self.ncout = nc4.Dataset(outfile, 'w')

   #-----------------------------------------------------------------------------------------
   #copy global attributes all at once via dictionary
   #ncout.setncatts(nc_gsi.__dict__)
    self.ncout.source='Round Robin Observer'
    self.ncout.comment = 'Concatenated Observation file from Round Robin Observer'

    self.openfiles()

  def __del__(self):
   #body of destructor
    self.ncout.close()

    for ncf in self.ncarray:
      ncf.close()

  def openfiles(self):
    nc = 0
    search_more = True
    while(search_more):
      filename = '%s/%s_%4.4d.nc4' %(self.obsdir, self.filetype, nc)
      if(os.path.exists(filename)):
       #print('File No %d: %s' %(nc, filename))
        self.filelist.append(filename)
        ncin = nc4.Dataset(filename, 'r')
        self.ncarray.append(ncin)
        nc += 1
      else:
        search_more = False

    self.create_dimension()

  def create_dimension(self):
    nf = len(self.ncarray)
    nc = 0
    dimlist = {}
    varlist = {}

    ncin = self.ncarray[0]
    for name, dimension in ncin.dimensions.items():
     #print('\tname: ', name)
     #print('\tdimension: ', dimension)

     #ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
      if dimension.isunlimited():
        self.ncout.createDimension(name, None)
      else:
        self.ncout.createDimension(name, len(dimension))

      if dimension.isunlimited():
        dimlist[name] = 0
      else:
        dimlist[name] = len(dimension)

    for name, variable in ncin.variables.items():
      print('name: ', name)
      print('variable.datatype: ', variable.datatype)
      print('variable.dimensions: ', variable.dimensions)
      x = self.ncout.createVariable(name, variable.datatype, variable.dimensions)

     #copy variable attributes all at once via dictionary
      self.ncout.variables[name].setncatts(ncin.variables[name].__dict__)

      val = ncin.variables[name][:]
      self.ncout.variables[name][:] = val[:]
      varlist[name] = len(val)

    for n in range(1, len(self.ncarray)):
      ncin = self.ncarray[n]
      nc = n + 1
      print('Processing No %d: %s' %(nc, self.filelist[n]))

     #-----------------------------------------------------------------------------------------
     #print('ncin.dimensions.keys() = ', ncin.dimensions.keys())
      for name, dimension in ncin.dimensions.items():
       #print('\tname: ', name)
       #print('\tdimension: ', dimension)

        if dimension.isunlimited():
          dimlist[name] = 0
        else:
          dimlist[name] += len(dimension)

     #print('ncin.variables.keys() = ', ncin.variables.keys())
      for name, variable in ncin.variables.items():
       #print('name: ', name)
       #print('variable: ', variable)
       #print('variable.datatype: ', variable.datatype)
       #print('variable.dimensions: ', variable.dimensions)

        val = ncin.variables[name][:]

        nb = varlist[name]
        ne = nb + len(val)
        self.ncout.variables[name][nb:ne] = val[:]

        varlist[name] += len(val)

  def process(self):
    nf = len(self.ncarray)
    nc = 0
    for ncin in self.ncarray:
      nc += 1

      for name, group in ncin.groups.items():
        print('name: ', name)
        print('group: ', group)

      for name, variable in ncin.variables.items():
        print('name: ', name)
        print('variable: ', variable)
       #if name in gsibaselist:
       #  x = nc4out.createVariable(name, variable.datatype, variable.dimensions)
       #  nc4out[name][:] = nc_gsi[name][:]
       # #copy variable attributes all at once via dictionary
       #  nc4out[name].setncatts(nc_gsi[name].__dict__)

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1

  obsdir = '/work2/noaa/gsienkf/weihuang/jedi/per_core_timing/run/anna_roundRobin_aircraft_2/run_80.40t1n_36p/obsout'
  filetype = 'aircraft_tsen_obs_2020011006'
  outfile = 'combined_obs_file.nc4'

#-----------------------------------------------------------------------------------------
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'obsdir=', 'filetype=',
                                                'outfile='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--obsdir'):
      obsdir = a
    elif o in ('--filetype'):
      filetype = a
    elif o in ('--outfile'):
      outfile = a
    else:
      assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
  mnf = MergeNetcdfFiles(debug=debug, obsdir=obsdir, filetype=filetype, outfile=outfile)


