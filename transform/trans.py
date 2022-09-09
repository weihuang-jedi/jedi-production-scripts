#=========================================================================
import sys
import netCDF4 as nc4

#-----------------------------------------------------------------------------------------
debug = 1

jedifile = 'xainc.20200101_120000z.nc4'
gsi_file = 'fv3_increment6.nc'
nco4incr = 'jedi_increment.nc4'

ncjedi = nc4.Dataset(jedifile, 'r')
nc_gsi = nc4.Dataset(gsi_file, 'r')
nc4out = nc4.Dataset(nco4incr, 'w')

#-----------------------------------------------------------------------------------------
gsibaselist = ['lon', 'lat', 'lev', 'ilev', 'hyai', 'hybi']

#copy global attributes all at once via dictionary
#nc4out.setncatts(nc_gsi.__dict__)
nc4out.source='JEDI getkf'
nc4out.comment = 'GSI increment converted from JEDI increment'

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

#-----------------------------------------------------------------------------------------
jedi_varlist = ['DZ', 'delp', 'o3mr', 'sphum', 'T', 'ua', 'va']
gsi_varlist = ['delz_inc', 'delp_inc', 'o3mr_inc',
               'sphum_inc', 'T_inc', 'u_inc', 'v_inc']

numvar = len(jedi_varlist)

#copy all var in gsi_varlist
for name, variable in nc_gsi.variables.items():
  if name in gsi_varlist:
    for n in range(numvar):
      gsivarname = gsi_varlist[n]
      if(name == gsivarname):
        jedivarname = jedi_varlist[n]
        print('\tjedivarname: %s, gsivarname: %s ' %(jedivarname, gsivarname))

        x = nc4out.createVariable(name, variable.datatype, variable.dimensions)
        var = ncjedi.variables[jedivarname][0,:,:,:]
        nc4out[name][:,:,:] = var
       #copy variable attributes all at once via dictionary
        nc4out[name].setncatts(nc_gsi[name].__dict__)

ncjedi.close()
nc_gsi.close()
nc4out.close()

