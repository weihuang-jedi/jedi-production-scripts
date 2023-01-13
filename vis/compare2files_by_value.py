import getopt
import os, sys
import numpy as np
import netCDF4 as nc4

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  frtdir = '/work2/noaa/da/weihuang/cycling/med.jedi_C96_lgetkf_sondesonly'
  snddir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly'

  datestr = '2020010112'

  frtfile = '%s/%s/sanl_2020010112_fhr06_ensmean' %(frtdir, datestr)
  sndfile = '%s/%s/sanl_2020010112_fhr06_ensmean' %(snddir, datestr)

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
  ncfrt = nc4.Dataset(frtfile, 'r')
  ncsnd = nc4.Dataset(sndfile, 'r')

#-----------------------------------------------------------------------------------------
  first_varlist = ['tmp', 'ugrd', 'vgrd', 'dpres', 'delz', 'spfh', 'o3mr']
  second_varlist = first_varlist
  unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)', 'Unit (Pa)',
              'Unit (kg/kg)', 'Unit (kg/kg)', 'Unit (m)']

#-----------------------------------------------------------------------------------------
  for n in range(len(first_varlist)):
    frtvar = ncfrt.variables[first_varlist[n]][:, :, :]
    sndvar = ncsnd.variables[second_varlist[n]][:, :, :]

   #print('frtvar.shape = ', frtvar.shape)
   #print('sndvar.shape = ', sndvar.shape)
   #nlev, nlat, nlon = sndvar.shape

    print('Var %d: %s' %(n, second_varlist[n]))

    difvar = sndvar - frtvar

    print('\tfirst.max: %f, first.min: %f' %(np.max(frtvar), np.min(frtvar)))
    print('\tsecnd.max: %f, secnd.min: %f' %(np.max(sndvar), np.min(sndvar)))
    print('\tdiffv.max: %f, diffv.min: %f' %(np.max(difvar), np.min(difvar)))

#-----------------------------------------------------------------------------------------
  ncsnd.close()
  ncfrt.close()

