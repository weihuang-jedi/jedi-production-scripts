import getopt
import os, sys
import numpy as np
import netCDF4 as nc4

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  topdir = '/work2/noaa/gsienkf/weihuang/production/run/sondes'
  incdir = 'analysis/increment'
  frtdir = '%s/new.1_rr_observer_whole_solver.run_81.36t9n/%s' %(topdir, incdir)
  snddir = '%s/whole.run_80.40t1n_36p/%s' %(topdir, incdir)
 #snddir = '%s/observer.run_80.40t1n_36p/%s' %(topdir, incdir)

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

#-----------------------------------------------------------------------------------------
  frtfile = '%s/xainc.20200101_120000z.nc4' %(frtdir)
  sndfile = '%s/xainc.20200101_120000z.nc4' %(snddir)

  ncfrt = nc4.Dataset(frtfile, 'r')
  ncsnd = nc4.Dataset(sndfile, 'r')

#-----------------------------------------------------------------------------------------
 #varlist = ['tmp', 'ugrd', 'vgrd', 'delp', 'spfh', 'o3mr', 'delz']
 #unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)', 'Unit (Pa)',
 #            'Unit (kg/kg)', 'Unit (ppm)', 'Unit (m)']
 #varlist = ['T_inc', 'delp_inc', 'sphum_inc', 'o3mr_inc', 'delz_inc']
  varlist = ['T', 'ua', 'va', 'delp', 'sphum', 'o3mr', 'DZ']
  unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)', 'Unit (Pa)',
              'Unit (kg/kg)', 'Unit (kg/kg)', 'Unit (m)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist)):
   #sndvar = ncsnd.variables[varlist[n]][:, :, :]
   #frtvar = ncfrt.variables[varlist[n]][:, :, :]
    sndvar = ncsnd.variables[varlist[n]][0, :, :, :]
    frtvar = ncfrt.variables[varlist[n]][0, :, :, :]

   #print('sndvar.shape = ', sndvar.shape)
   #print('frtvar.shape = ', frtvar.shape)
   #nlev, nlat, nlon = sndvar.shape

    print('Var %d: %s' %(n, varlist[n]))

    difvar = sndvar - frtvar

    print('\tfirst.max: %f, first.min: %f' %(np.max(frtvar), np.min(frtvar)))
    print('\tsecnd.max: %f, secnd.min: %f' %(np.max(sndvar), np.min(sndvar)))
    print('\tdiffv.max: %f, diffv.min: %f' %(np.max(difvar), np.min(difvar)))

#-----------------------------------------------------------------------------------------
  ncsnd.close()
  ncfrt.close()

