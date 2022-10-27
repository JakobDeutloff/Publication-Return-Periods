import xarray as xr
import matplotlib.pyplot as plt
from xclim.indices.stats import frequency_analysis
import cartopy.crs as ccrs
import cartopy.feature as ft


# %% Load concatenated and corrected GPM data
GPM = xr.open_dataset('GPM/GPM_climatology_2000_2021.nc4')

# %% extract GB region
GPM_GB = GPM.sel(lon=slice(-11, 2), lat=slice(49.5, 60))

#%% plot GB region and compare both GPM measures
fig, axes = plt.subplots(2, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5, 7))

GPM_GB.HQprecipitation.isel(time=1).plot(ax=axes[0], x='lon', y='lat', transform=ccrs.PlateCarree())
GPM_GB.precipitationCal.isel(time=1).plot(ax=axes[1], x='lon', y='lat', transform=ccrs.PlateCarree())

for ax in axes:
    ax.coastlines()
    ax.add_feature(ft.BORDERS)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
plt.savefig('Plots/GPM_data.png', bbox_inches='tight')
plt.show()

#%% Calculate 10 Yr Return Periods
rp_10 = frequency_analysis(da=GPM_GB.precipitationCal, mode="max", freq='YS', t=10, dist="genextreme")

rp_10_HQ = frequency_analysis(da=GPM_GB.HQprecipitation, mode="max", freq='YS', t=10, dist="genextreme")
#%% Remove invalid values
rp_10_corr = rp_10.where(rp_10 < 500)
rp_10_HQ_corr = rp_10_HQ.where(rp_10 < 500)

# %% save datasets
rp_10_corr.to_netcdf(path='Data/GPM/ReturnPeriods/rp10.nc4', format='NETCDF4')
rp_10_HQ_corr.to_netcdf(path='Data/GPM/ReturnPeriods/rp10_HQ.nc4', format='NETCDF4')




