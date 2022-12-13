import xarray as xr
import cartopy.crs as ccrs
from cartopy import feature
import numpy as np
import matplotlib.pyplot as plt
from Code.visualisation.Plot_prob_of_exceedance import Read_ensemble_GER, Read_ensemble_UK
from matplotlib.colors import LinearSegmentedColormap


def Plot_ensemble_ERA(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    if land == 'GER':
        Ens = Read_ensemble_GER('ERA')

    if land == 'UK':
        Ens = Read_ensemble_UK('ERA', event)

    fig, axes = plt.subplots(6, 3, figsize=(10, 20), subplot_kw={'projection': ccrs.Miller()})

    levels = [2.5, 7.5, 12.5, 37.5, 62.5, 137.5]
    for i in range(18):
        ax = axes.flatten()[i]
        cm = Ens.sel(ens_mem=i).return_period.plot.contourf(ax=ax, levels=levels, cmap='YlOrRd',
                                                            transform=ccrs.PlateCarree(), add_colorbar=False,
                                                            vmin=5, vmax=100, extend='neither')
        ax.coastlines()
        ax.add_feature(feature.OCEAN, zorder=1, color='grey')
        ax.add_feature(feature.BORDERS)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25, 50, 100])

    if event:
        plt.savefig('Plots/' + land + '/casestudy/ERA_members' + event + '.png', bbox_inches='tight')
    else:
        plt.savefig('Plots/' + land + '/casestudy/ERA_members.png', bbox_inches='tight')


def Plot_ensemble_GPM(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    if land == 'GER':
        Ens = Read_ensemble_GER('GPM')

    if land == 'UK':
        Ens = Read_ensemble_UK('GPM', event)

    fig, axes = plt.subplots(6, 3, figsize=(10, 20), subplot_kw={'projection': ccrs.Miller()})

    levels = [2.5, 7.5, 12.5, 37.5]
    colors = [(255, 255, 204), (254, 217, 118), (253, 140, 60)]
    cmap = LinearSegmentedColormap.from_list('GPM', colors)
    for i in range(18):
        ax = axes.flatten()[i]
        cm = Ens.sel(ens_mem=i).return_period.plot.contourf(ax=ax, levels=levels, cmap=cmap,
                                                            transform=ccrs.PlateCarree(), add_colorbar=False,
                                                            vmin=5, vmax=25, extend='neither')
        ax.coastlines()
        ax.add_feature(feature.OCEAN, zorder=1, color='grey')
        ax.add_feature(feature.BORDERS)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25])

    if event:
        plt.savefig('Plots/' + land + '/casestudy/GPM_members' + event + '.png', bbox_inches='tight')
    else:
        plt.savefig('Plots/' + land + '/casestudy/GPM_members.png', bbox_inches='tight')


# %%

if __name__ == '__main__':
    Plot_ensemble_ERA('GER')
    Plot_ensemble_GPM('GER')

    plt.show()
