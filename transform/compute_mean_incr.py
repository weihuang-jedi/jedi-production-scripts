#=========================================================================
import os, sys
import getopt
import netCDF4 as nc4

#=========================================================================
class ComputeMeanIncrements():
  def __init__(self, debug=0, dir=None, datestr=None,
               casename='GSI', outfilename='mean_incr.nc4'):
    self.debug = debug
    self.dir = dir
    self.datestr = datestr
    self.casename = casename
    self.totalmembers = 80

   #Base variables
    self.baselist = ['lon', 'lat', 'lev', 'ilev', 'hyai', 'hybi']

   #open all members increments.
    self.filelist = []
    self.ncarray = []
    for n in range(self.totalmembers):
      flnm = '%s/%s/mem%3.3d/INPUT/fv3_increment6.nc' %(dir, datestr, n+1)
      if(os.path.exists(flnm)):
       #print('Processing No %d: %s' %(n, flnm))
        ncf = nc4.Dataset(flnm, 'r')
        self.filelist.append(flnm)
        self.ncarray.append(ncf)
      else:
        print('Increment file: %s does not exist. Stop' %(flnm))
        sys.exit(-1)

    self.ncout = nc4.Dataset(outfilename, 'w')

    self.process()

  def __del__(self):
    self.ncout.close()
    for ncf in self.ncarray:
      ncf.close()
  
  def process(self):
   #copy global attributes all at once via dictionary
    self.ncout.source = '%s increments mean from: %s, date: %s' %(self.casename, self.dir, self.datestr)
    self.ncout.comment = 'Calculated for %d members' %(self.totalmembers)

    ncf = self.ncarray[0]
   #print('ncf: ', ncf)
   #copy dimensions
    for name, dimension in ncf.dimensions.items():
     #print('name: ', name)
     #print('dimension: ', dimension)
     #self.ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
      if dimension.isunlimited():
        self.ncout.createDimension(name, None)
      else:
        self.ncout.createDimension(name, len(dimension))

    self.vardict = {}
   #copy all var in baselist
    for name, variable in ncf.variables.items():
      print('name: ', name)
     #print('variable: ', variable)
      if name in self.baselist:
        print('name %s is in self.baselist' %(name))
        x = self.ncout.createVariable(name, variable.datatype, variable.dimensions)
        self.ncout.variables[name][:] = ncf.variables[name][:]
       #copy variable attributes all at once via dictionary
        self.ncout.variables[name].setncatts(ncf.variables[name].__dict__)
      else:
        print('name %s is not in self.baselist' %(name))
        var = ncf.variables[name][:,:,:]
        self.vardict[name] = var

    print('vardict.keys(): ', self.vardict.keys())
    for n in range(1, len(self.ncarray)):
      ncf = self.ncarray[n]
      for name in self.vardict.keys():
        var = ncf.variables[name][:,:,:]
        self.vardict[name] += var

    for name, variable in ncf.variables.items():
      print('working on: name: ', name)
      if name in self.baselist:
        continue
      else:
        x = self.ncout.createVariable(name, variable.datatype, variable.dimensions)
        var = self.vardict[name] / self.totalmembers
        self.ncout.variables[name][:,:,:] = var
       #copy variable attributes all at once via dictionary
        self.ncout.variables[name].setncatts(ncf.variables[name].__dict__)

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1

  dir = '/work2/noaa/gsienkf/weihuang/gsi/gsi_C96_lgetkf_sondesonly'
  datestr = '2020010206'
  casename = 'GSI'
  outfilename = 'mean_incr.nc4'

 #-----------------------------------------------------------------------------------------
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'dir=', 'datestr=',
                                                'casename=', 'outfilename='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--dir'):
      dir = a
    elif o in ('--datestr'):
      datestr = a
    elif o in ('--outfilename'):
      outfilename = a
    elif o in ('--casename'):
      casename = a
    elif o in ('--outfilename'):
      outfilename = a
    else:
      assert False, 'unhandled option'

 #-----------------------------------------------------------------------------------------
  cmi = ComputeMeanIncrements(debug=debug, dir=dir, datestr=datestr,
                              casename=casename, outfilename=outfilename)

