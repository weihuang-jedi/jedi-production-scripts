import os
from pathlib import Path
from solo.basic_files import tree
from solo.ioda import Ioda

def concat_files():
    workdir = '/work2/noaa/gsienkf/weihuang/jedi/run/rr_maxpoolsize_tpe/run_80.40t1n_36p/obsout'
    obstype = 'aircraft_tsen_obs_2020011006'
    output_file_path = './%s.nc4' %(obstype)

    output = Path(output_file_path)

    netcdf_files = []
   
    nc = 0
    search_more = True
    while(search_more):
      filename = '%s/%s_%4.4d.nc4' %(workdir, obstype, nc)
      if(os.path.exists(filename)):
        netcdf_files.append(filename)
        nc += 1
      else:
        search_more = False

    nc = Ioda('obs name')
    nc.concat_files(netcdf_files, str(output.absolute()))


if __name__ == '__main__':
    concat_files()

