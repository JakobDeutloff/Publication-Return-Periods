"""
Extrect the UK and Germany regions from ERA5 return period data
"""

import xarray as xr

# %% Define functions

def extr_region(da, lat_min, lat_max, lon_min, lon_max):
    """
    Returns specified region subset of ERA5 data
    """
    return da.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))

# %% open ERA5 datasets
rp5 = xr.open_dataset(
    'Data/ERA/Raw/precipitation-at-fixed-return-period_europe_ecad_30-year_5-yrs_1989-2018_v1.nc')
rp10 = xr.open_dataset(
    'Data/ERA/Raw/precipitation-at-fixed-return-period_europe_ecad_30-year_10-yrs_1989-2018_v1.nc')
rp25 = xr.open_dataset(
    'Data/ERA/Raw/precipitation-at-fixed-return-period_europe_ecad_30-year_25-yrs_1989-2018_v1.nc')
rp50 = xr.open_dataset(
    'Data/ERA/Raw/precipitation-at-fixed-return-period_europe_ecad_30-year_50-yrs_1989-2018_v1.nc')
rp100 = xr.open_dataset(
    'Data/ERA/Raw/precipitation-at-fixed-return-period_europe_ecad_30-year_100-yrs_1989-2018_v1.nc')

# %% Cut out UK and GER regions and save datasets
periods = [5, 10, 25, 50, 100]
rps = [rp5, rp10, rp25, rp50, rp100]

# UK
for i in range(len(rps)):
    UK = extr_region(rps[i], lon_min=-11, lon_max=2, lat_max=49.5, lat_min=60)
    UK.to_netcdf(path='Data/ERA/UK/ERA_UK_rp'+str(periods[i])+'.nc4', format='NETCDF4')

# GER
for i in range(len(rps)):
    GER = extr_region(rps[i], lon_min=5, lon_max=15.5, lat_max=47, lat_min=55)
    GER.to_netcdf(path='Data/ERA/GER/ERA_GER_rp'+str(periods[i])+'.nc4', format='NETCDF4')
