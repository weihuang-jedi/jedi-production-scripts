#=========================================================================
import sys
import numpy as np
import netCDF4 as nc4

#-----------------------------------------------------------------------------------------
debug = 1

jediori = 'xainc.20200101_120000z.nc4'
jedinew = 'jedi_increment.nc4'

ncori = nc4.Dataset(jediori, 'r')
ncnew = nc4.Dataset(jedinew, 'r')

ori_varlist = ['DZ', 'delp', 'o3mr', 'sphum', 'T', 'ua', 'va']
new_varlist = ['delz_inc', 'delp_inc', 'o3mr_inc',
               'sphum_inc', 'T_inc', 'u_inc', 'v_inc']

numvar = len(ori_varlist)

for n in range(numvar):
  orivarname = ori_varlist[n]
  newvarname = new_varlist[n]
  print('\torivarname: %s, newvarname: %s ' %(orivarname, newvarname))

  orivar = ncori.variables[orivarname][0,:,:,:]
  newvar = ncnew.variables[newvarname][:,:,:]

  diff = orivar - newvar

  print('diff max: %f, min: %f' %(np.max(diff), np.min(diff)))

ncori.close()
ncnew.close()

