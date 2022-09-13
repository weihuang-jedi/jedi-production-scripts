import os, sys, time
import multiprocessing

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

def is_perfect(n):
    sum_factors = 0
    for i in range(1, n):
        if(n % i == 0):
            sum_factors = sum_factors + i
    if (sum_factors == n):
        print('{} is a Perfect number'.format(n))

def sleepy_man(sec):
    print('Starting to sleep for {} seconds'.format(sec))
    time.sleep(sec)
    print('Done sleeping for {} seconds'.format(sec))

if __name__ == '__main__':
 #concat_files()

  tic = time.time()
  pool = multiprocessing.Pool()
  pool.map(is_perfect, range(1,100000))
  pool.close()
  toc = time.time()

  print('is_perfect done in {:.4f} seconds'.format(toc-tic))

  tic = time.time()
  pool = multiprocessing.Pool(5)
  pool.map(sleepy_man, range(1,11))
  pool.close()
  toc = time.time()

  print('sleepy_man done in {:.4f} seconds'.format(toc-tic))
