import fsspec
import h5py
import xarray as xr

#fs = fsspec.filesystem('s3', requester_pays=True, profile='esip-qhub')
#files = fs.glob('s3://coawst-public/rsignell/testing/*.nc')
#files

infile = 'sondes_tsen_obs_2020011006_0000.nc4'
outfile = 'new_sondes_tsen_obs_2020011006_0000.nc4'

#url = files[0]
#f = h5py.File(fs.open(url))
#[key for key in f.keys()]

f = h5py.File(infile)
varlist = [key for key in f.keys()]
#print('varlist = ', varlist)

#ds = xr.open_dataset(infile)
ds = xr.open_dataset(infile, group='hofx_y_mean_xb0')
#print('ds: ', ds)
#print('ds.keys: ', ds.keys)

#varlist = [i.name for i in ds.data_vars]
varlist = [i for i in ds.data_vars]
print('varlist = ', varlist)
#print('ds.data_vars = ', ds.data_vars)

#coordlist = list(ds.coords)
#print('coordlist = ', coordlist)

#grp = ds.groups['ombg']
#grp = ds.groups['hofx_y_mean_xb0']

for var in varlist:
  val = ds.variables[var]
  val += 1.0
 #ds.variables[var].value = val
  ds[var] = val

#ds.to_netcdf(outfile, group='hofx_y_mean_xb0', mode='a')
ds.to_netcdf(infile, group='hofx_y_mean_xb0', mode='a')

