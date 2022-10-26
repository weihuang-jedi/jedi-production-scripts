#=========================================================================
import os
import sys
import getopt
import netCDF4 as nc4
import time
import numpy as np

#-----------------------------------------------------------------------------------------
def average(vlist):
  buf = np.zeros_like(vlist[0])
  nmem = len(vlist)
  for val in vlist:
    buf += val
  buf /= nmem
  return buf

#-----------------------------------------------------------------------------------------
def copy_dimension(ncin, ncout):
 #copy dimensions
  for name, dimension in ncin.dimensions.items():
   #ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
    if dimension.isunlimited():
      ncout.createDimension(name, None)
    else:
      ncout.createDimension(name, len(dimension))

#-----------------------------------------------------------------------------------------
def copy_rootvar(ncin, ncout):
 #copy all var in root group.
  for name, variable in ncin.variables.items():
    ncout.createVariable(name, variable.datatype, variable.dimensions)
    ncout[name][:] = ncin[name][:]
   #copy variable attributes all at once via dictionary
   #ncout[name].setncatts(ncin[name].__dict__)

#-----------------------------------------------------------------------------------------
def copy_var_in_group(ncingroup, ncoutgroup):
 #copy all var in group.
  for varname, variable in ncingroup.variables.items():
    ncoutgroup.createVariable(varname, variable.datatype, variable.dimensions)
   #copy variable attributes all at once via dictionary
   #ncoutgroup[varname].setncatts(ncingroup[varname].__dict__)
    ncoutgroup[varname][:] = ncingroup[varname][:]

#-----------------------------------------------------------------------------------------
def copy_grp2newname(name, n, group, ncout):
  item = name.split('_')
  item[-1] = '%d' %(n+1)
  newname = '_'.join(item)
  print('No %d name: %s, newname: %s' %(n+1, name, newname))
  ncoutgroup = ncout.createGroup(newname)
  copy_var_in_group(group, ncoutgroup)

#-----------------------------------------------------------------------------------------
def process(ncinlist, ncout, grplist):
  namelist = []
  hofxgrps = []
  commongrps = []
  ensvarinfo = {}

 #check groups
  for name, group in ncinlist[0].groups.items():
    print('name: ', name)
   #print('group: ', group)
    namelist.append(name)
    if(name in grplist):
     #print('0. name: ', name)
      ensvarinfo[name] = {}
    else:
      if(name.find('hofx') < 0):
       #print('1. name: ', name)
        commongrps.append(name)
        ncoutgroup = ncout.createGroup(name)
        copy_var_in_group(group, ncoutgroup)
      else:
       #print('2. name: ', name)
        hofxgrps.append(name)

  print('len(namelist) = %d, len(commongrps) = %d, len(hofxgrps) = %d' %(len(namelist), len(commongrps), len(hofxgrps)))

  for n in range(len(ncinlist)):
    ncin = ncinlist[n]
    for name in hofxgrps:
      group = ncin.groups[name]
      copy_grp2newname(name, n, group, ncout)
    for name in grplist:
      group = ncin.groups[name]
      if(name == 'hofx0_1'):
        copy_grp2newname(name, n, group, ncout)
      for varname, variable in group.variables.items():
        val = group[varname][:]
        if(n == 0):
          ensvarinfo[name][varname] = []
        ensvarinfo[name][varname].append(val)

  varlist = ensvarinfo['hofx0_1'].keys()
  print('varlist = ', varlist)
  meanvars = {}
  for grpname in grplist:
    meanvars[grpname] = {}
    ncoutgroup = ncout.createGroup(grpname)
    for varname in varlist:
      meanval = average(ensvarinfo[grpname][varname])
      print('meanval.shape = ', meanval.shape)
      print('meanval.size = ', meanval.size)
      meanvars[grpname][varname] = meanval

  for grpname in grplist:
    if(grpname != 'hofx0_1'):
      print('Finally processing: ', grpname)
      ncingroup = ncinlist[0].groups[grpname]
      ncoutgroup = ncout.createGroup(grpname)
      for varname, variable in ncingroup.variables.items():
        ncoutgroup.createVariable(varname, variable.datatype, variable.dimensions)
       #copy variable attributes all at once via dictionary
       #ncoutgroup[name].setncatts(ncingroup[name].__dict__)
        if(grpname == 'ombg'):
          val = meanval + meanvars['hofx_y_mean_xb0'][varname] - meanvars['hofx0_1'][varname]
          ncoutgroup[varname][:] = val
        else:
          ncoutgroup[varname][:] = meanval

#-----------------------------------------------------------------------------------------
def replace_var(filelist, outfile, grplist):
  ncinlist = []
  for infile in filelist:
    if(os.path.exists(infile)):
     #print('infile: ', infile)
      ncin = nc4.Dataset(infile, 'r')
      ncinlist.append(ncin)
    else:
      print('infile: %s does not exist.' %(infile))
      sys.exit(-1)

  if(os.path.exists(outfile)):
    os.remove(outfile)
    print('outfile: ', outfile)
 #else:
 #  print('outfile: %s does not exist.' %(outfile))

  ncout = nc4.Dataset(outfile, 'w')

 #copy global attributes all at once via dictionary
 #ncout.setncatts(ncin.__dict__)
  ncout.source='JEDI observer only ouptut, each with only one member'
  ncout.comment = 'updated variable hofx_y_mean_xb0 and ombg from all observer files'

  copy_dimension(ncinlist[0], ncout)
  copy_rootvar(ncinlist[0], ncout)

  process(ncinlist, ncout, grplist)
 
  for ncin in ncinlist:
    ncin.close()
  ncout.close()

#-----------------------------------------------------------------------------------------
debug = 1
nmem = 80

run_dir = '/work2/noaa/gsienkf/weihuang/production/run/sondes/run_81.36t9n'
datestr = '2020010112'
varname = 'sondes_tsen'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'run_dir=',
                                              'datestr=', 'varname=', 'nmem='])
for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--run_dir'):
    run_dir = a
  elif o in ('--datestr'):
    datestr = a
  elif o in ('--varname'):
    varname = a
  elif o in ('--nmem'):
    nmem = int(a)
  else:
    assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
grplist = ['hofx_y_mean_xb0', 'hofx0_1', 'ombg']
filelist = []
for n in range(nmem):
  flnm = '%s/observer/mem%3.3d/%s_obs_%s_0000.nc4' %(run_dir, n+1, varname, datestr)
  filelist.append(flnm)
outfile = '%s/observer/%s_obs_%s_0000.nc4' %(run_dir, varname, datestr)

replace_var(filelist, outfile, grplist)

