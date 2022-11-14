import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as ft

# %% load GPM data
GPM = xr.open_dataset('Data/GPM/GPM_climatology_2000_2021.nc4')
# Extract UK region
GPM_GB = GPM.sel(lat=slice(49.5, 60), lon=slice(-11, 2))

# %% plot GB region and compare both GPM measures
fig, axes = plt.subplots(2, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5, 7))

GPM_GB.HQprecipitation.isel(time=1).plot(ax=axes[0], x='lon', y='lat', transform=ccrs.PlateCarree())
GPM_GB.precipitationCal.isel(time=1).plot(ax=axes[1], x='lon', y='lat', transform=ccrs.PlateCarree())

for ax in axes:
    ax.coastlines()
    ax.add_feature(ft.BORDERS)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5,
                      linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
plt.savefig('Plots/GPM_data.png', bbox_inches='tight')
plt.show()
