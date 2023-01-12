#=========================================================================
import os, sys
import getopt
import netCDF4 as nc4

#=========================================================================
class ComputeMeanIncrements():
  def __init__(self, debug=0, totalmembers=80):
    self.debug = debug
    self.totalmembers = totalmembers

    print('debug: ', debug)

   #Base variables
    self.baselist = ['lon', 'lat', 'lev', 'ilev', 'hyai', 'hybi']

  def process(self, name='GSI', workdir=None, outfile=None):
    print('workdir: ', workdir)
    print('outfile: ', outfile)

   #open all members increments.
    filelist = []
    ncarray = []
    for n in range(self.totalmembers):
      filename = '%s/mem%3.3d/INPUT/fv3_increment6.nc' %(workdir, n+1)
      if(os.path.exists(filename)):
        print('Processing No %d: %s' %(n, filename))
        ncf = nc4.Dataset(filename, 'r')
        filelist.append(filename)
        ncarray.append(ncf)
      else:
        print('increment file: %s does not exist. Stop' %(filename))
        sys.exit(-1)

    ncout = nc4.Dataset(outfile, 'w')
   #copy global attributes all at once via dictionary
    ncout.source = '%s increments mean from: %s' %(name, workdir)
    ncout.comment = 'Calculated for %d members' %(self.totalmembers)

    ncf = ncarray[0]
   #print('ncf: ', ncf)
   #copy dimensions
    for name, dimension in ncf.dimensions.items():
     #print('name: ', name)
     #print('dimension: ', dimension)
     #ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
      if dimension.isunlimited():
        ncout.createDimension(name, None)
      else:
        ncout.createDimension(name, len(dimension))

    self.vardict = {}
   #copy all var in baselist
    for name, variable in ncf.variables.items():
     #print('name: ', name)
     #print('variable: ', variable)
      if name in self.baselist:
        x = ncout.createVariable(name, variable.datatype, variable.dimensions)
        ncout.variables[name][:] = ncf.variables[name][:]
       #copy variable attributes all at once via dictionary
       #ncout.variables[name].setncatts(ncf.variables[name].__dict__)
      else:
        var = ncf.variables[name][:,:,:]
        self.vardict[name] = var

    print('self.vardict.keys(): ', self.vardict.keys())
    for n in range(1, len(ncarray)):
      ncf = ncarray[n]
      for name in self.vardict.keys():
        var = ncf.variables[name][:,:,:]
        self.vardict[name] += var

    ncf = ncarray[0]
    for name, variable in ncf.variables.items():
      print('working on: name: ', name)
      if name in self.baselist:
        continue
      else:
        x = ncout.createVariable(name, variable.datatype, variable.dimensions)
        var = self.vardict[name] / self.totalmembers
        ncout.variables[name][:,:,:] = var
       #copy variable attributes all at once via dictionary
       #ncout.variables[name].setncatts(ncf.variables[name].__dict__)

    ncout.close()
    for ncf in ncarray:
      ncf.close()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1

 #topdir = '/work2/noaa/gsienkf/weihuang/gsi'
  topdir = '/work2/noaa/da/weihuang/cycling'
  datestr = '2020010112'
  casename = 'jedi'

 #-----------------------------------------------------------------------------------------
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'topdir=',
                            'datestr=', 'casename='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--topdir'):
      topdir = a
    elif o in ('--datestr'):
      datestr = a
    elif o in ('--casename'):
      casename = a
    else:
      assert False, 'unhandled option'

 #-----------------------------------------------------------------------------------------
  workdir = '%s/%s_C96_lgetkf_sondesonly/%s' %(topdir, casename, datestr)
  outfile = '%s_mean_incr_%s.nc4' %(casename, datestr)

 #-----------------------------------------------------------------------------------------
  cmi = ComputeMeanIncrements(debug=debug)
  cmi.process(name=casename, workdir=workdir, outfile=outfile)

