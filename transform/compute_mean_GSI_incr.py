#=========================================================================
import os, sys
import getopt
import netCDF4 as nc4

#=========================================================================
class ComputeGSImeanIncrements():
  def __init__(self, debug=0, gsidir=None, datestr=None):
    self.debug = debug
    self.gsidir = gsidir
    self.datestr = datestr
    self.totalmembers = 80

   #Base variables
    self.baselist = ['lon', 'lat', 'lev', 'ilev', 'hyai', 'hybi']

   #open all GSI members increments.
    self.filelist = []
    self.ncarray = []
    for n in range(self.totalmembers):
      gsiflnm = '%s/%s/mem%3.3d/INPUT/fv3_increment6.nc' %(gsidir, datestr, n+1)
      if(os.path.exists(gsiflnm)):
       #print('Processing No %d: %s' %(n, gsiflnm))
        ncf = nc4.Dataset(gsiflnm, 'r')
        self.filelist.append(gsiflnm)
        self.ncarray.append(ncf)
      else:
        print('GSI increment file: %s does not exist. Stop' %(gsiflnm))
        sys.exit(-1)

    outfilename = 'gsi_mean_incr.nc4'
    self.ncout = nc4.Dataset(outfilename, 'w')

    self.process()

  def __del__(self):
    self.ncout.close()
    for ncf in self.ncarray:
      ncf.close()
  
  def process(self):
   #copy global attributes all at once via dictionary
    self.ncout.source = 'GSI increments mean from: %s, date: %s' %(self.gsidir, self.datestr)
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
     #print('name: ', name)
     #print('variable: ', variable)
      if name in self.baselist:
        x = self.ncout.createVariable(name, variable.datatype, variable.dimensions)
        self.ncout.variables[name][:] = ncf.variables[name][:]
       #copy variable attributes all at once via dictionary
        self.ncout.variables[name].setncatts(ncf.variables[name].__dict__)
      else:
        var = ncf.variables[name][:,:,:]
        self.vardict[name] = var

    print('vardict.keys(): ', vardict.keys())
    for n in range(1, len(self.ncarray)):
      ncf = self.ncarray[n]
      for name in self.vardict.keys():
        var = ncf.variables[name][:,:,:]
        self.vardict[name] += var

    ncf = self.ncarray[0]
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

  gsidir = '/work2/noaa/gsienkf/weihuang/gsi/gsi_C96_lgetkf_sondesonly'
  datestr = '2020010200'

 #-----------------------------------------------------------------------------------------
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'gsidir=',
                                                'datestr='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--gsidir'):
      gsidir = a
    elif o in ('--datestr'):
      datestr = a
    else:
      assert False, 'unhandled option'

 #-----------------------------------------------------------------------------------------
  cgi = ComputeGSImeanIncrements(debug=debug, gsidir=gsidir, datestr=datestr)
