#=========================================================================
import sys
import getopt
import netCDF4 as nc4

#-----------------------------------------------------------------------------------------
debug = 1

gsi_file = 'fv3_increment6.nc'
gridfile = 'gaussian_grid.nc4'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'gsi_file=', 'gridfile='])
for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--gsi_file'):
    gsi_file = a
  elif o in ('--gridfile'):
    gridfile = a
  else:
    assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------

nc_gsi = nc4.Dataset(gsi_file, 'r')
nc4out = nc4.Dataset(gridfile, 'w')

#-----------------------------------------------------------------------------------------
gsibaselist = ['lon', 'lat', 'lev', 'ilev', 'hyai', 'hybi']

#copy global attributes all at once via dictionary
#nc4out.setncatts(nc_gsi.__dict__)
nc4out.source='GSI Gaussian grid output'
nc4out.comment = 'Gaussian grid only'

#copy dimensions
for name, dimension in nc_gsi.dimensions.items():
 #nc4out.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
  if dimension.isunlimited():
    nc4out.createDimension(name, None)
  else:
    nc4out.createDimension(name, len(dimension))

#copy all var in gsibaselist
for name, variable in nc_gsi.variables.items():
  if name in gsibaselist:
    x = nc4out.createVariable(name, variable.datatype, variable.dimensions)
    nc4out[name][:] = nc_gsi[name][:]
   #copy variable attributes all at once via dictionary
    nc4out[name].setncatts(nc_gsi[name].__dict__)

nc_gsi.close()
nc4out.close()

