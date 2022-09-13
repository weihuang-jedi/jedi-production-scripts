import netCDF4
import numpy
import xarray

#This is the main line of code.
#tell the script to open any netCDF file in my folder that begins with ‘precip.hour’,
#the * representing that whatever comes there is a don’t-care string.
#The way we want to combine the data is using the time coordinate dimension,
#which is the dimension that defines the relations of the data and order
#with respect to one another. And so we say combine = ‘by_coords’
#and the concat_dim = ‘time’. (xarray provides more option to concatenate data
#with more complex needs like merging along two dimensions. For more, see here)

ds = xarray.open_mfdataset('/work2/noaa/gsienkf/weihuang/jedi/per_core_timing/run/anna_roundRobin_aircraft_2/run_80.40t1n_36p/obsout/aircraft_tsen_obs_2020011006_*.nc4',combine = 'by_coords', concat_dim="nlocs")

#I want to export this data into a combined netCDF, so that’s next.

ds.to_netcdf('aircraft_tsen_obs_2020011006.nc4')

