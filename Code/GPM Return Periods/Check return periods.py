"""
This script contains some code to inspect and plot the GPM return periods to see if the calculation procedure is valid.
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as ft
from matplotlib import colors
import cmocean as ccm
from xclim.indices.generic import select_resample_op
from xclim.indices.stats import fit, parametric_quantile
from scipy.stats import genextreme

# %% set font size
matplotlib.rcParams.update({'font.size': 12})

# %% read GMP and ERA5 data
GPM = xr.open_dataset('Data/GPM/GPM_climatology_2000_2021.nc4')
GPM_rp10 = xr.open_dataset('Data/GPM/ReturnPeriods/rp10.nc4')
ERA_rp10 = xr.open_dataset('Data/ERA/UK/ERA_UK_rp10.nc4')

# Drop unnecessary dim from GPM rp
GPM_rp10 = GPM_rp10.sel(return_period=10)
GPM_rp10 = GPM_rp10.rename({'lon':'longitude', 'lat':'latitude'})
GPM_rp10 = GPM_rp10.reset_coords('return_period', drop=True)

#%% Plot RP in comparision to ERA5

# Set Colorbar
max = 120
min = 40
levels = 21
cmap = ccm.cm.rain
# Plot
fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 5.5), sharey=True)
axes[0].add_feature(ft.OCEAN, zorder=1, color='grey')
axes[1].add_feature(ft.OCEAN, zorder=1, color='grey')
GPM_rp10.precipitationCal.plot.contourf(ax=axes[0], x='longitude', y='latitude', transform=ccrs.PlateCarree(),
                                        vmax=max, vmin=min, levels=levels, add_colorbar=False, cmap=cmap, extend='both')
cm = ERA_rp10.r10yrrp.plot.contourf(ax=axes[1], transform=ccrs.PlateCarree(), vmax=max, vmin=min, levels=levels,
                                    add_colorbar=False, cmap=cmap, extend='both')

# Add Grid and Labels
for ax in axes:
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_xlim(-10.95, 1.95)
    ax.set_ylim(49.55, 59.95)
gl.left_labels = False

# Add Colorbar
fig.subplots_adjust(bottom=0.13)
cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Daily Precipitation [mm]')

# Titles
axes[0].set_title('GPM')
axes[1].set_title('ERA5')


# Save and show
plt.savefig('Plots/10yRP_GPM_ERA.png', bbox_inches='tight')
plt.show()


# %% Plot difference

# calculate difference
diff = np.flip(GPM_rp10.precipitationCal.values.T, axis=0) - ERA_rp10.r10yrrp.values
lat, lon = np.meshgrid(ERA_rp10.latitude.values, ERA_rp10.longitude.values)

# norm for colorbar
divnorm = colors.TwoSlopeNorm(vmin=-20., vcenter=0, vmax=120)
levels=np.arange(-100, 120, 20)
# Plot
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(7, 4))
ax.add_feature(ft.OCEAN, zorder=1, color='grey')
cl = ax.contourf(lon, lat, diff.T, transform=ccrs.PlateCarree(), levels=levels, cmap='coolwarm', extend='max')
cb = fig.colorbar(cl, label='Daily Precipitation [mm]')
ax.set_title('10yr RP GPM - ERA5')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Coordinated of point with high diff
ax.plot(-4.95, 56.45, marker='o', markersize=2, color='lime')

plt.savefig('Plots/GPM-ERA_rp10.png', bbox_inches='tight')
plt.show()

#%% get  timeseries of plotted point in GPM and estimates for rp

GPM_point = GPM.sel(lon=-4.95, lat=56.45, method='nearest')
GPM_rp10_point = GPM_rp10.sel(longitude=-4.95, latitude=56.45, method='nearest')
ERA_rp10_point = ERA_rp10.sel(longitude=-4.95, latitude=56.45, method='nearest')

#%% plot timeseries of GPM data with annual maxima

# Calculate annual maxima
sub = select_resample_op(GPM_point.precipitationCal, op='max', freq='Y')

# get date of annual maxima
dates = []
for max in sub.values:
    idx = np.where(GPM_point.precipitationCal.values == max)
    dates.append(GPM_point.time[idx].values)

# Plot
fig, ax = plt.subplots()
GPM_point.precipitationCal.plot(ax=ax, linewidth=0.5, color='grey', alpha=0.5)
ax.plot(dates, sub.values, linestyle='', marker='o', markersize='2', color='red', label='annual max')
ax.set_xlabel('')
ax.set_title('GPM data at 4.95W 56.45N')
ax.set_ylabel('Daily Precipitation [mm]')
ax.legend()
plt.tight_layout()
plt.savefig('Plots/GPM_timeseries.png')
plt.show()

# %% fit distribution

# parameters of genextreme pdf
params = fit(sub, dist='genextreme')

# Calculate quantile
qt = parametric_quantile(params, q=1 - 1.0 / 10)

# %% plot distribution
x = np.arange(50, 230)
y = genextreme.pdf(x, c=-0.18675725, loc=86.17427517, scale=25.58816589)

fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.hist(sub.values, color='grey', alpha=1)
ax.set_ylabel('Frequency', color='grey')
ax2.plot(x, y, color='red')
ax2.set_ylabel('Probability Density', color='red')
ax.set_xlabel('Daily Precipitation [mm]')
ax.set_title('GPM data at 4.95W 56.45N')
plt.tight_layout()
plt.savefig('Plots/GPM_genextreme.png')
plt.show()