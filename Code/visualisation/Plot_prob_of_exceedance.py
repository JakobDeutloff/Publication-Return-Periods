import xarray as xr
import cartopy.crs as ccrs
from cartopy import feature
import numpy as np
import matplotlib.pyplot as plt


#  define functions
def Read_ensemble_GER(dataset, event):
    Ens = xr.open_mfdataset('Data/' + dataset + '/GER/' + event + '/*.nc', combine='nested',
                            concat_dim=['ens_mem'])
    Ens = Ens.assign_coords({'ens_mem': np.arange(0, 18)})
    Ens = Ens.rename({'unknown': 'return_period'})
    return Ens


def Read_ensemble_UK(dataset, event):
    Ens = xr.open_mfdataset('Data/' + dataset + '/UK/' + event + '/*.nc', combine='nested',
                            concat_dim=['ens_mem'])
    Ens = Ens.assign_coords({'ens_mem': np.arange(0, 18)})
    Ens = Ens.rename({'unknown': 'return_period'})
    return Ens


def Read_ensemble_ciara(dataset, event):
    Ens = xr.open_mfdataset('Data/' + dataset + '/UK/' + event + '/*.nc', combine='nested',
                            concat_dim=['ens_mem'])
    Ens = Ens.assign_coords({'ens_mem': np.arange(0, 12)})
    Ens = Ens.rename({'unknown': 'return_period'})
    return Ens


def Calculate_exc_prob(Ens, rp):
    # Calculate probability of exceedance
    Exc = Ens.return_period >= rp
    Prob = (Exc.sum('ens_mem') / len(Ens.coords['ens_mem'])) * 100
    return Prob


def Calculate_number_of_members(Ens, rp):
    Exc = Ens.return_period >= rp
    Num = Exc.sum('ens_mem')
    return Num


def Plot_exceedance_prob(Prob, ax, title, color, transform, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    levels = np.arange(0, 110, 10)
    levels[0] = 1

    cm = Prob.plot.contourf(ax=ax, transform=transform, levels=levels, cmap='PuRd', add_colorbar=False,
                            extend='neither')

    ax.coastlines()
    ax.add_feature(feature.OCEAN, zorder=1, color='grey')
    ax.add_feature(feature.BORDERS)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.set_title(title + ' Year RP', color=color)
    return cm


def Plot_number_of_members(Num, ax, title, transform, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    levels = np.arange(0.5, 19.5, 1)

    cm = Num.plot.contourf(ax=ax, transform=transform, levels=levels, cmap='PuRd', add_colorbar=False,
                           extend='neither')

    ax.coastlines()
    ax.add_feature(feature.OCEAN, zorder=1, color='grey')
    ax.add_feature(feature.BORDERS)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.set_title(title + ' Year RP')
    return cm


def Plot_exc_prob_GPM(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        GPM_ens = Read_ensemble_GER('GPM', event)
        transform = ccrs.PlateCarree()
        figsize = (11, 5)

    if land == 'UK':
        GPM_ens = Read_ensemble_UK('GPM', event)
        central_latitude = GPM_ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = GPM_ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (11, 3.5)

    # Calculate Prob of Exceedance

    ex_prob = {}
    ex_prob['5'] = Calculate_exc_prob(GPM_ens, 5)
    ex_prob['10'] = Calculate_exc_prob(GPM_ens, 10)
    ex_prob['25'] = Calculate_exc_prob(GPM_ens, 25)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()})
    rps = list(ex_prob.keys())
    for i in range(3):
        cm = Plot_exceedance_prob(ex_prob[rps[i]], axes[i], title=rps[i], transform=transform,
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl = axes[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False

    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability of >= RP [%]')
    plt.savefig('Plots/' + land + '/casestudy/GPM_ex_prob' + event + '.png', bbox_inches='tight')


def Plot_num_memb_GPM(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        GPM_ens = Read_ensemble_GER('GPM', event)
        transform = ccrs.PlateCarree()
        figsize = (11, 5)

    if land == 'UK':
        GPM_ens = Read_ensemble_UK('GPM', event)
        central_latitude = GPM_ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = GPM_ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (11, 3.5)

    # Calculate Prob of Exceedance

    num = {}
    num['5'] = Calculate_number_of_members(GPM_ens, 5)
    num['10'] = Calculate_number_of_members(GPM_ens, 10)
    num['25'] = Calculate_number_of_members(GPM_ens, 25)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharey='row',
                             sharex='col')
    rps = list(num.keys())
    for i in range(3):
        cm = Plot_number_of_members(num[rps[i]], axes[i], title=rps[i], transform=transform,
                                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl = axes[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False

    # Add Colorbar
    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    cbar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Numer of Ensemble Members')
    cbar.set_ticks(np.arange(1, 19))
    plt.savefig('Plots/' + land + '/casestudy/GPM_num_members' + event + '.png', bbox_inches='tight')


def Plot_exc_prob_GPM_ciara(event='202002080300', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):
    # Read data
    GPM_ens = Read_ensemble_ciara('GPM', event)
    merc = GPM_ens.transverse_mercator
    transform = ccrs.LambertAzimuthalEqualArea(
        central_latitude=merc.latitude_of_projection_origin,
        central_longitude=merc.longitude_of_central_meridian,
        false_easting=merc.false_easting,
        false_northing=merc.false_northing,
    )
    figsize = (11, 4.5)

    # Calculate Prob of Exceedance

    ex_prob = {}
    ex_prob['5'] = Calculate_exc_prob(GPM_ens, 5)
    ex_prob['10'] = Calculate_exc_prob(GPM_ens, 10)
    ex_prob['25'] = Calculate_exc_prob(GPM_ens, 25)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()})
    rps = list(ex_prob.keys())
    for i in range(3):
        cm = Plot_exceedance_prob(ex_prob[rps[i]], axes[i], title=rps[i], transform=transform,
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl = axes[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False

    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability of >= RP [%]')
    plt.savefig('Plots/UK/ciara/GPM_ex_prob' + event + '.png', bbox_inches='tight')


def Plot_num_memb_GPM_ciara(event='202002080300', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):
    # Read data
    GPM_ens = Read_ensemble_ciara('GPM', event)
    merc = GPM_ens.transverse_mercator
    transform = ccrs.LambertAzimuthalEqualArea(
        central_latitude=merc.latitude_of_projection_origin,
        central_longitude=merc.longitude_of_central_meridian,
        false_easting=merc.false_easting,
        false_northing=merc.false_northing,
    )
    figsize = (11, 4.5)

    # Calculate Prob of Exceedance

    num = {}
    num['5'] = Calculate_number_of_members(GPM_ens, 5)
    num['10'] = Calculate_number_of_members(GPM_ens, 10)
    num['25'] = Calculate_number_of_members(GPM_ens, 25)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharey='row',
                             sharex='col')
    rps = list(num.keys())
    for i in range(3):
        cm = Plot_number_of_members(num[rps[i]], axes[i], title=rps[i], transform=transform,
                                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl = axes[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False

    # Add Colorbar
    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    cbar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Numer of Ensemble Members')
    cbar.set_ticks(np.arange(1, 19))
    plt.savefig('Plots/UK/ciara/GPM_num_members' + event + '.png', bbox_inches='tight')


def Plot_exc_prob_ERA(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        ERA_ens = Read_ensemble_GER('ERA', event)
        transform = ccrs.PlateCarree()
        figsize = (11, 9)

    if land == 'UK':
        ERA_ens = Read_ensemble_UK('ERA', event)
        central_latitude = ERA_ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = ERA_ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (11, 6)

    # Calculate Prob of Exceedance
    ex_prob = {}
    ex_prob['5'] = Calculate_exc_prob(ERA_ens, 5)
    ex_prob['10'] = Calculate_exc_prob(ERA_ens, 10)
    ex_prob['25'] = Calculate_exc_prob(ERA_ens, 25)
    ex_prob['50'] = Calculate_exc_prob(ERA_ens, 50)
    ex_prob['100'] = Calculate_exc_prob(ERA_ens, 100)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()})
    rps = list(ex_prob.keys())
    for i in range(5):
        cm = Plot_exceedance_prob(ex_prob[rps[i]], axes.flatten()[i], title=rps[i], transform=transform,
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    for ax in axes[:, 0]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for ax in axes[1, :]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    axes[0, 1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')
    axes[0, 2].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')

    axes[1, 2].remove()

    # Add Colorbar
    fig.subplots_adjust(bottom=0.18)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability >= RP [%]')
    plt.savefig('Plots/' + land + '/casestudy/ERA_ex_prob' + event + '.png', bbox_inches='tight')


def Plot_num_memb_ERA(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        ERA_ens = Read_ensemble_GER('ERA', event)
        transform = ccrs.PlateCarree()
        figsize = (11, 9)

    if land == 'UK':
        ERA_ens = Read_ensemble_UK('ERA', event)
        central_latitude = ERA_ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = ERA_ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (11, 6)

    # Calculate Prob of Exceedance
    num = {}
    num['5'] = Calculate_number_of_members(ERA_ens, 5)
    num['10'] = Calculate_number_of_members(ERA_ens, 10)
    num['25'] = Calculate_number_of_members(ERA_ens, 25)
    num['50'] = Calculate_number_of_members(ERA_ens, 50)
    num['100'] = Calculate_number_of_members(ERA_ens, 100)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharex='col',
                             sharey='row')
    rps = list(num.keys())
    for i in range(5):
        cm = Plot_number_of_members(num[rps[i]], axes.flatten()[i], title=rps[i], transform=transform,
                                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of axes
    for ax in axes[:, 0]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for ax in axes[1, :]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    axes[0, 1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')
    axes[0, 2].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')

    axes[1, 2].remove()

    # Add Colorbar
    fig.subplots_adjust(bottom=0.18)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    cbar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Number of Ensemble members')
    cbar.set_ticks(np.arange(1, 19))
    plt.savefig('Plots/' + land + '/casestudy/ERA_num_members' + event + '.png', bbox_inches='tight')

def Plot_exc_prob_ERA_ciara(event='202002080300', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):

    ERA_ens = Read_ensemble_ciara('ERA', event)
    merc = ERA_ens.transverse_mercator
    transform = ccrs.LambertAzimuthalEqualArea(
        central_latitude=merc.latitude_of_projection_origin,
        central_longitude=merc.longitude_of_central_meridian,
        false_easting=merc.false_easting,
        false_northing=merc.false_northing)
    figsize = (10, 8)

    # Calculate Prob of Exceedance
    ex_prob = {}
    ex_prob['5'] = Calculate_exc_prob(ERA_ens, 5)
    ex_prob['10'] = Calculate_exc_prob(ERA_ens, 10)
    ex_prob['25'] = Calculate_exc_prob(ERA_ens, 25)
    ex_prob['50'] = Calculate_exc_prob(ERA_ens, 50)
    ex_prob['100'] = Calculate_exc_prob(ERA_ens, 100)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()})
    rps = list(ex_prob.keys())
    for i in range(5):
        cm = Plot_exceedance_prob(ex_prob[rps[i]], axes.flatten()[i], title=rps[i], transform=transform,
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    for ax in axes[:, 0]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for ax in axes[1, :]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    axes[0, 1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')
    axes[0, 2].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')

    axes[1, 2].remove()

    # Add Colorbar
    fig.subplots_adjust(bottom=0.18)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability >= RP [%]')
    plt.savefig('Plots/UK/ciara/ERA_ex_prob' + event + '.png', bbox_inches='tight')


def Plot_num_memb_ERA_ciara(event='202002080300', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):
    # Read data

    ERA_ens = Read_ensemble_ciara('ERA', event)
    merc = ERA_ens.transverse_mercator
    transform = ccrs.LambertAzimuthalEqualArea(
        central_latitude=merc.latitude_of_projection_origin,
        central_longitude=merc.longitude_of_central_meridian,
        false_easting=merc.false_easting,
        false_northing=merc.false_northing)
    figsize = (10, 8)

    # Calculate Prob of Exceedance
    num = {}
    num['5'] = Calculate_number_of_members(ERA_ens, 5)
    num['10'] = Calculate_number_of_members(ERA_ens, 10)
    num['25'] = Calculate_number_of_members(ERA_ens, 25)
    num['50'] = Calculate_number_of_members(ERA_ens, 50)
    num['100'] = Calculate_number_of_members(ERA_ens, 100)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharex='col',
                             sharey='row')
    rps = list(num.keys())
    for i in range(5):
        cm = Plot_number_of_members(num[rps[i]], axes.flatten()[i], title=rps[i], transform=transform,
                                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of axes
    for ax in axes[:, 0]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for ax in axes[1, :]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    axes[0, 1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')
    axes[0, 2].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')

    axes[1, 2].remove()

    # Add Colorbar
    fig.subplots_adjust(bottom=0.18)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    cbar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Number of Ensemble members')
    cbar.set_ticks(np.arange(1, 19))
    plt.savefig('Plots/UK/ciara/ERA_num_members' + event + '.png', bbox_inches='tight')

def Plot_exc_prob_GPM_and_ERA(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55, figsize=(10, 10)):
    # Read data GPM
    if land == 'GER':
        GPM_ens = Read_ensemble_GER('GPM', event)
        transform = ccrs.PlateCarree()

    if land == 'UK':
        GPM_ens = Read_ensemble_UK('GPM', event)
        central_latitude = GPM_ens.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = GPM_ens.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)

    # Read data ERA
    if land == 'GER':
        ERA_ens = Read_ensemble_GER('ERA', event)

    if land == 'UK':
        ERA_ens = Read_ensemble_UK('ERA', event)


    # Calculate Prob of Exceedance
    ex_prob_GPM = {}
    ex_prob_GPM['5'] = Calculate_exc_prob(GPM_ens, 5)
    ex_prob_GPM['10'] = Calculate_exc_prob(GPM_ens, 10)
    ex_prob_GPM['25'] = Calculate_exc_prob(GPM_ens, 25)

    # Calculate Prob of Exceedance
    ex_prob_ERA = {}
    ex_prob_ERA['5'] = Calculate_exc_prob(ERA_ens, 5)
    ex_prob_ERA['10'] = Calculate_exc_prob(ERA_ens, 10)
    ex_prob_ERA['25'] = Calculate_exc_prob(ERA_ens, 25)
    ex_prob_ERA['50'] = Calculate_exc_prob(ERA_ens, 50)
    ex_prob_ERA['100'] = Calculate_exc_prob(ERA_ens, 100)


    fig, axes = plt.subplots(2, 4, figsize=figsize, subplot_kw={'projection': ccrs.Miller()})
    # Plot ERA
    rps_ERA = list(ex_prob_ERA.keys())
    for i in range(5):
        cm = Plot_exceedance_prob(ex_prob_ERA[rps_ERA[i]], axes.flatten()[i], title=rps_ERA[i], color='green',
                                  transform=transform,
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Plot GPM
    rps_GPM = list(ex_prob_GPM.keys())
    for i in range(3):
        Plot_exceedance_prob(ex_prob_GPM[rps_GPM[i]], axes.flatten()[i+5], title=rps_GPM[i], color='blue', transform=transform,
                             lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    # Set grid of left ax
    for ax in axes[:, 0]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5,
                          linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False

    for ax in axes[1, :]:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5,
                          linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False

    axes[0, 1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')
    axes[0, 2].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')
    axes[0, 3].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='k', alpha=0.5,
                         linestyle='--')

    fig.subplots_adjust(bottom=0.21)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability of >= RP [%]')
    plt.savefig('Plots/' + land + '/' + event + '/GPM_ERA_ex_prob.png', bbox_inches='tight', dpi=500)

# %%
if __name__ == '__main__':
    #Plot_exc_prob_GPM_and_ERA(land='GER', event='20210711T1800Z', figsize=(11, 6))
    Plot_exc_prob_GPM_and_ERA(land='UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53,
                              figsize=(11, 4))
    #Plot_exc_prob_GPM(land='GER', event='20210711T1800Z')
    #Plot_num_memb_GPM(land='GER', event='20210711T1800Z')
    #Plot_exc_prob_ERA(land='GER', event='20210711T1800Z')
    #Plot_num_memb_ERA(land='GER', event='20210711T1800Z')

    # Plot_exc_prob_GPM(land='UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)
    # Plot_num_memb_GPM(land='UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)
    # Plot_exc_prob_ERA(land='UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)
    # Plot_num_memb_ERA(land='UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)

    #Plot_exc_prob_GPM_ciara(lon_min=-6, lon_max=2, lat_max=50, lat_min=56)
    #Plot_num_memb_GPM_ciara(lon_min=-6, lon_max=2, lat_max=50, lat_min=56)
    #Plot_exc_prob_ERA_ciara(lon_min=-6, lon_max=2, lat_max=50, lat_min=56)
    #Plot_num_memb_ERA_ciara(lon_min=-6, lon_max=2, lat_max=50, lat_min=56)
