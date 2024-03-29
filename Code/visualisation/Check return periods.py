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


# %% def plot function

def plot_RP(ax, rp, title, c_max, c_min, levels, lat_min, lat_max, lon_min, lon_max):
    # Set Colorbar
    cmap = ccm.cm.rain

    # Plot
    cm = rp.plot.contourf(ax=ax, x='longitude', y='latitude', transform=ccrs.PlateCarree(), vmax=c_max, vmin=c_min,
                          levels=levels, add_colorbar=False, cmap=cmap, extend='both')

    # Plot Features
    ax.add_feature(ft.OCEAN, zorder=1, color='grey')
    ax.add_feature(ft.BORDERS)

    # Set Gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    ax.set_title(title)

    return cm, gl


def plot_difference(ax, GPM_rp, ERA_rp, title, lat_min, lat_max, lon_min, lon_max):
    # calculate difference
    diff = np.flip(GPM_rp.values.T, axis=0) - ERA_rp.values
    lat, lon = np.meshgrid(ERA_rp.latitude.values, ERA_rp.longitude.values)

    # set cmap
    levels = np.arange(-100, 120, 20)
    cmap = 'coolwarm'

    # Plot
    cl = ax.contourf(lon, lat, diff.T, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, extend='max')

    # Plot Features
    ax.add_feature(ft.OCEAN, zorder=1, color='grey')
    ax.add_feature(ft.BORDERS)

    # Set grid
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    ax.set_title(title)

    return cl, gl


def plot_ERA_GPM(ERA_rp, GPM_rp, rp, c_max, c_min, levels, land, lat_min, lat_max, lon_min, lon_max):
    # Set Colorbar
    cmap = ccm.cm.rain

    # Plot
    fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 5.5), sharey=True)
    axes[0].add_feature(ft.OCEAN, zorder=1, color='grey')
    axes[1].add_feature(ft.OCEAN, zorder=1, color='grey')
    axes[1].add_feature(ft.BORDERS)
    axes[0].add_feature(ft.BORDERS)
    GPM_rp.precipitationCal.plot.contourf(ax=axes[0], x='longitude', y='latitude', transform=ccrs.PlateCarree(),
                                          vmax=c_max, vmin=c_min, levels=levels, add_colorbar=False, cmap=cmap,
                                          extend='both')
    cm = ERA_rp.plot.contourf(ax=axes[1], transform=ccrs.PlateCarree(), vmax=c_max, vmin=c_min, levels=levels,
                              add_colorbar=False, cmap=cmap, extend='both')

    # Add Grid and Labels
    for ax in axes:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
    gl.left_labels = False

    # Add Colorbar
    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Daily Precipitation [mm]')

    # Titles
    axes[0].set_title('GPM')
    axes[1].set_title('ERA5')

    # Save and show
    plt.savefig('Plots/' + land + '/' + str(rp) + 'yRP_GPM_ERA.png', bbox_inches='tight')
    plt.show()


def read_GPM(rp, land):
    # read GMP and ERA5 data
    GPM = xr.open_dataset('Data/GPM/' + land + '/GPM_' + land + '_rp' + str(rp) + '.nc4')

    # Drop unnecessary dim from GPM rp
    GPM = GPM.sel(return_period=rp)
    GPM = GPM.rename({'lon': 'longitude', 'lat': 'latitude'})
    GPM = GPM.reset_coords('return_period', drop=True)

    return GPM


# %% Read ERA5 and GPM data

GPM = xr.open_dataset('Data/GPM/GPM_climatology_2000_2021.nc4')
rps = [5, 10, 25]
GPM_UK = {}
GPM_GER = {}
ERA_UK = {}
ERA_GER = {}

# read GMP and ERA5 data
for rp in rps:
    # UK
    GPM_UK[str(rp)] = read_GPM(rp, 'UK')
    ERA_UK[str(rp)] = xr.open_dataset('Data/ERA/UK/ERA_UK_rp' + str(rp) + '.nc4')['r' + str(rp) + 'yrrp']
    # GER
    GPM_GER[str(rp)] = read_GPM(rp, 'GER')
    ERA_GER[str(rp)] = xr.open_dataset('Data/ERA/GER/ERA_GER_rp' + str(rp) + '.nc4')['r' + str(rp) + 'yrrp']

#%% Plot GPM and ERA5 data for all rps
for rp in rps:
    # UK
    plot_ERA_GPM(ERA_UK[str(rp)], GPM_UK[str(rp)], rp, c_min=20, c_max=180, levels=21, land='UK',
                 lon_min=-11, lon_max=2, lat_min=49.5, lat_max=60)
    # GER
    plot_ERA_GPM(ERA_GER[str(rp)], GPM_GER[str(rp)], rp, c_min=20, c_max=180, levels=21, land='GER',
                 lon_min=5, lon_max=15.5, lat_min=47, lat_max=55)

# %% Plot return periods and difference

# coords scotland
lon_scot = -4.9
lat_scot = 56.5

# set rp
rp = '25'

# Setup figure and initialize arrays
fig, axes = plt.subplots(2, 3, figsize=(10, 8), sharey='row', subplot_kw={'projection': ccrs.Miller()})
cms = np.empty(axes.shape).tolist()
gls = np.empty(axes.shape).tolist()

# Plot UK
cms[0][0], gls[0][0] = plot_RP(axes[0, 0], ERA_UK[rp], title='EPI', c_min=20, c_max=180, levels=21, lon_min=-8,
                               lon_max=2, lat_min=49.9, lat_max=58.5)
cms[0][1], gls[0][1] = plot_RP(axes[0, 1], GPM_UK[rp].precipitationCal, title='GPM', c_min=20, c_max=180, levels=21,
                               lon_min=-8, lon_max=2, lat_min=49.9, lat_max=58.5)
gls[0][1].left_labels = False
cms[0][2], gls[0][2] = plot_difference(axes[0, 2], GPM_UK[rp].precipitationCal, ERA_UK[rp], title='GPM - EPI',
                                       lon_min=-8, lon_max=2, lat_min=49.9, lat_max=58.5)
gls[0][2].left_labels = False
axes[0, 1].plot(lon_scot, lat_scot, color='red', marker='o', transform=ccrs.PlateCarree())

# Plot GER
cms[1][0], gls[1][0] = plot_RP(axes[1, 0], ERA_GER[rp], title='EPI', c_min=20, c_max=180, levels=21, lon_min=5, lon_max=15.5,
                               lat_min=47, lat_max=55)
cms[1][1], gls[1][1] = plot_RP(axes[1, 1], GPM_GER[rp].precipitationCal, title='GPM', c_min=20, c_max=180, levels=21,
                               lon_min=5, lon_max=15.5, lat_min=47, lat_max=55)
gls[1][1].left_labels = False
cms[1][2], gls[1][2] = plot_difference(axes[1, 2], GPM_GER[rp].precipitationCal, ERA_GER[rp], title='GPM - EPI',
                                       lon_min=5, lon_max=15.5, lat_min=47, lat_max=55)
gls[1][2].left_labels = False

# Add colorbars
fig.subplots_adjust(bottom=0.16)
cbar_ax1 = fig.add_axes([0.13, 0.1, 0.5, 0.04])
cbar_ax2 = fig.add_axes([0.68, 0.1, 0.2, 0.04])
fig.colorbar(cms[0][0], cax=cbar_ax1, orientation='horizontal', pad=0.01, label='Daily Precipitation [mm]')
fig.colorbar(cms[0][2], cax=cbar_ax2, orientation='horizontal', pad=0.01, label='Daily Precipitation [mm]')
plt.savefig('Plots/Comparison_rp' + rp +'.png', bbox_inches='tight')
plt.show()
# %% Extract specific values
val_ERA = ERA_UK['25'].sel(latitude=lat_scot, longitude=lon_scot, method='nearest')
val_GPM = GPM_UK['25'].precipitationCal.sel(latitude=lat_scot, longitude=lon_scot, method='nearest')
diff = val_GPM.values - val_ERA.values


# %% Plot difference

GPM_rp10 = GPM_UK['10']
ERA_rp10 = ERA_UK['10']

# calculate difference
diff = np.flip(GPM_rp10.precipitationCal.values.T, axis=0) - ERA_rp10.values
lat, lon = np.meshgrid(ERA_rp10.latitude.values, ERA_rp10.longitude.values)

# norm for colorbar
divnorm = colors.TwoSlopeNorm(vmin=-20., vcenter=0, vmax=120)
levels = np.arange(-100, 120, 20)
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

# %% get  timeseries of plotted point in GPM and estimates for rp

GPM_point = GPM.sel(lon=-4.95, lat=56.45, method='nearest')
GPM_rp10_point = GPM_rp10.sel(longitude=-4.95, latitude=56.45, method='nearest')
ERA_rp10_point = ERA_rp10.sel(longitude=-4.95, latitude=56.45, method='nearest')

# %% plot timeseries of GPM data with annual maxima

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
