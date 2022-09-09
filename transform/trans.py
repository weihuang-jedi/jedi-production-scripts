#=========================================================================
import sys
import netCDF4 as nc4

#-----------------------------------------------------------------------------------------
debug = 1

jedifile = 'xainc.20200101_120000z.nc4'
gsi_file = 'fv3_increment6.nc'

ncjedi = nc4.Dataset(jedifile, 'r')
nc_gsi = nc4.Dataset(gsi_file, 'r+')

jedi_varlist = ['DZ', 'delp', 'o3mr', 'sphum', 'T', 'ua', 'va']
gsi_varlist = ['delz_inc', 'delp_inc', 'o3mr_inc',
               'sphum_inc', 'T_inc', 'u_inc', 'v_inc']

numvar = len(jedi_varlist)

for n in range(numvar):
  jedivarname = jedi_varlist[n]
  gsivarname = gsi_varlist[n]
  print('\tjedivarname: %s, gsivarname: %s ' %(jedivarname, gsivarname))

  var = ncjedi.variables[jedivarname][0,:,:,:]

  nc_gsi.variables[gsivarname][:,:,:] = var

ncjedi.close()
nc_gsi.close()

