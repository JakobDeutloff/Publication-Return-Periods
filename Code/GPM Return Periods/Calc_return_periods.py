import xarray as xr
from xclim.indices.stats import frequency_analysis

# %% Load concatenated and corrected GPM data
GPM = xr.open_dataset('Data/GPM/GPM_climatology_2000_2021.nc4')

#%% Calculate 5, 10 & 25 yr Return Periods UK and GER

# Set regions
lon_UK = slice(-11, 2)
lat_UK = slice(49.5, 60)
lon_GER = slice(5, 15.5)
lat_GER = slice(47, 55)

# RPs
rps = [5, 10, 25]

# UK
for rp in rps:
    retper = frequency_analysis(da=GPM.precipitationCal.sel(lon=lon_UK, lat=lat_UK),
                                mode="max", freq='YS', t=rp, dist="genextreme")
    #Remove invalid values
    retper_corr = retper.where(retper < 1800)
    #save dataset
    retper_corr.to_netcdf(path='Data/GPM/UK/GPM_UK_rp'+str(rp)+'.nc4', format='NETCDF4')

# GER
for rp in rps:
    retper = frequency_analysis(da=GPM.precipitationCal.sel(lon=lon_GER, lat=lat_GER),
                                mode="max", freq='YS', t=rp, dist="genextreme")
    # Remove invalid values
    retper_corr = retper.where(retper < 1800)
    # save dataset
    retper_corr.to_netcdf(path='Data/GPM/GER/GPM_GER_rp'+str(rp)+'.nc4', format='NETCDF4')








