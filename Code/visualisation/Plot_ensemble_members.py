import xarray as xr
import cartopy.crs as ccrs
from cartopy import feature
import numpy as np
import matplotlib.pyplot as plt
from Code.visualisation.Plot_prob_of_exceedance import Read_ensemble_GER, Read_ensemble_UK, Read_ensemble_ciara
from matplotlib.colors import LinearSegmentedColormap


def Plot_ensemble_ERA(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    if land == 'GER':
        Ens = Read_ensemble_GER('ERA', event)
        transform = ccrs.PlateCarree()
        figsize = (10, 20)
    if land == 'UK':
        Ens = Read_ensemble_UK('ERA', event)
        central_latitude = Ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = Ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (10, 18)

    fig, axes = plt.subplots(6, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()})

    levels = [2.5, 7.5, 12.5, 37.5, 62.5, 137.5]
    for i in range(18):
        ax = axes.flatten()[i]
        cm = Ens.sel(ens_mem=i).return_period.plot.contourf(ax=ax, levels=levels, cmap='YlOrRd',
                                                            transform=transform, add_colorbar=False,
                                                            vmin=5, vmax=100, extend='neither')
        ax.coastlines()
        ax.add_feature(feature.OCEAN, zorder=1, color='grey')
        ax.add_feature(feature.BORDERS)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # draw gridlines
    for j in range(5):
        ax = axes[j, 0]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for k in range(2):
        ax = axes[5, k + 1]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    for l in range(5):
        for m in range(2):
            ax = axes[l, m + 1]
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5,
                              linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = False
            gl.bottom_labels = False

    ax = axes[5, 0]
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25, 50, 100])

    plt.savefig('Plots/' + land + '/casestudy/ERA_members' + event + '.png', bbox_inches='tight')


def Plot_ensemble_GPM(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    if land == 'GER':
        Ens = Read_ensemble_GER('GPM', event)
        transform = ccrs.PlateCarree()
        figsize = (10, 20)
    if land == 'UK':
        Ens = Read_ensemble_UK('GPM', event)
        central_latitude = Ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = Ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (10, 18)

    fig, axes = plt.subplots(6, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharex='col',
                             sharey='row')

    levels = [2.5, 7.5, 12.5, 37.5]
    colors = [(0 / 25, (255/255, 255/255, 204/255)), (10 / 25, (254/255, 217/255, 118/255)), (25 / 25, (253/255, 140/255, 60/255))]
    cmap = LinearSegmentedColormap.from_list('GPM', colors)
    for i in range(18):
        ax = axes.flatten()[i]
        cm = Ens.sel(ens_mem=i).return_period.plot.contourf(ax=ax, levels=levels, cmap=cmap,
                                                            transform=transform, add_colorbar=False,
                                                            vmin=5, vmax=25, extend='neither')
        ax.coastlines()
        ax.add_feature(feature.OCEAN, zorder=1, color='grey')
        ax.add_feature(feature.BORDERS)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # draw gridlines
    for j in range(5):
        ax = axes[j, 0]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for k in range(2):
        ax = axes[5, k + 1]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    for l in range(5):
        for m in range(2):
            ax = axes[l, m + 1]
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5,
                              linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = False
            gl.bottom_labels = False

    ax = axes[5, 0]
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25])

    plt.savefig('Plots/' + land + '/casestudy/GPM_members' + event + '.png', bbox_inches='tight')

def Plot_ensemble_GPM_ciara(event='202002080300', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):

    Ens = Read_ensemble_ciara('GPM', event)
    merc = Ens.transverse_mercator
    transform = ccrs.LambertAzimuthalEqualArea(
        central_latitude=merc.latitude_of_projection_origin,
        central_longitude=merc.longitude_of_central_meridian,
        false_easting=merc.false_easting,
        false_northing=merc.false_northing,
        )
    figsize = (9, 12)

    fig, axes = plt.subplots(4, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharex='col',
                             sharey='row')

    levels = [2.5, 7.5, 12.5, 37.5]
    colors = [(0 / 25, (255/255, 255/255, 204/255)), (10 / 25, (254/255, 217/255, 118/255)), (25 / 25, (253/255, 140/255, 60/255))]
    cmap = LinearSegmentedColormap.from_list('', colors)
    for i in range(12):
        ax = axes.flatten()[i]
        cm = Ens.sel(ens_mem=i).return_period.plot.contourf(ax=ax, levels=levels, cmap=cmap,
                                                            transform=transform, add_colorbar=False,
                                                            vmin=5, vmax=25, extend='neither')
        ax.coastlines()
        ax.add_feature(feature.OCEAN, zorder=1, color='grey')
        ax.add_feature(feature.BORDERS)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # draw gridlines
    for j in range(3):
        ax = axes[j, 0]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for k in range(2):
        ax = axes[3, k + 1]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    for l in range(3):
        for m in range(2):
            ax = axes[l, m + 1]
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5,
                              linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = False
            gl.bottom_labels = False

    ax = axes[3, 0]
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25])

    plt.savefig('Plots/UK/ciara/GPM_members.png', bbox_inches='tight')

def Plot_ensemble_ERA_ciara(event='202002080300', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):

    Ens = Read_ensemble_ciara('ERA', event)
    merc = Ens.transverse_mercator
    transform = ccrs.LambertAzimuthalEqualArea(
        central_latitude=merc.latitude_of_projection_origin,
        central_longitude=merc.longitude_of_central_meridian,
        false_easting=merc.false_easting,
        false_northing=merc.false_northing,
        )
    figsize = (9, 12)

    fig, axes = plt.subplots(4, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharex='col',
                             sharey='row')

    levels = [2.5, 7.5, 12.5, 37.5, 62.5, 137.5]
    for i in range(12):
        ax = axes.flatten()[i]
        cm = Ens.sel(ens_mem=i).return_period.plot.contourf(ax=ax, levels=levels, cmap='YlOrRd',
                                                            transform=transform, add_colorbar=False,
                                                            vmin=5, vmax=100, extend='neither')
        ax.coastlines()
        ax.add_feature(feature.OCEAN, zorder=1, color='grey')
        ax.add_feature(feature.BORDERS)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # draw gridlines
    for j in range(3):
        ax = axes[j, 0]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for k in range(2):
        ax = axes[3, k + 1]
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    for l in range(3):
        for m in range(2):
            ax = axes[l, m + 1]
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5,
                              linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = False
            gl.bottom_labels = False

    ax = axes[3, 0]
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25, 50, 100])

    plt.savefig('Plots/UK/ciara/ERA_members.png', bbox_inches='tight')


# %%

if __name__ == '__main__':
    #Plot_ensemble_ERA('GER', event='20210711T1800Z')
    Plot_ensemble_GPM('GER', event='20210711T1800Z')
    #Plot_ensemble_ERA('UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)
    Plot_ensemble_GPM('UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)
    #Plot_ensemble_GPM_ciara(lon_min=-6, lon_max=2, lat_max=50, lat_min=56)
    # Plot_ensemble_ERA_ciara(lon_min=-6, lon_max=2, lat_max=50, lat_min=56)
    plt.show()
