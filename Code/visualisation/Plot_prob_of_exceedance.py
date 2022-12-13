import xarray as xr
import cartopy.crs as ccrs
from cartopy import feature
import numpy as np
import matplotlib.pyplot as plt


#  define functions

def Read_ensemble_GER(dataset):
    Ens = xr.open_mfdataset('Data/' + dataset + '/GER/RP_output/*.nc', combine='nested',
                            concat_dim=['ens_mem'])
    Ens = Ens.assign_coords({'ens_mem': np.arange(0, 18)})
    return Ens

def Read_ensemble_UK(dataset, event):
    Ens = xr.open_mfdataset('Data/' + dataset + '/UK/' + event + '/*.nc', combine='nested',
                            concat_dim=['ens_mem'])
    Ens = Ens.assign_coords({'ens_mem': np.arange(0, 18)})
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


def Plot_exceedance_prob(Prob, ax, title, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    levels = np.arange(0, 110, 10)
    levels[0] = 1

    cm = Prob.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap='YlOrBr', add_colorbar=False,
                            extend='neither')

    ax.coastlines()
    ax.add_feature(feature.OCEAN, zorder=1, color='grey')
    ax.add_feature(feature.BORDERS)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.set_title(title + ' Year RP')
    return cm


def Plot_number_of_members(Num, ax, title, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    levels = np.arange(0.5, 19.5, 1)

    cm = Num.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap='YlOrBr', add_colorbar=False,
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
        GPM_ens = Read_ensemble_GER('ERA')

    if land == 'UK':
        GPM_ens = Read_ensemble_UK('ERA', event)

    # Calculate Prob of Exceedance

    ex_prob = {}
    ex_prob['5'] = Calculate_exc_prob(GPM_ens, 5)
    ex_prob['10'] = Calculate_exc_prob(GPM_ens, 10)
    ex_prob['25'] = Calculate_exc_prob(GPM_ens, 25)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(11, 5), subplot_kw={'projection': ccrs.Miller()})
    rps = list(ex_prob.keys())
    for i in range(3):
        cm = Plot_exceedance_prob(ex_prob[rps[i]], axes[i], title=rps[i],
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability of >= RP [%]')
    plt.savefig('Plots/' + land + '/casestudy/GPM_ex_prob.png', bbox_inches='tight')


def Plot_num_memb_GPM(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        GPM_ens = Read_ensemble_GER('ERA')

    if land == 'UK':
        GPM_ens = Read_ensemble_UK('ERA', event)

    # Calculate Prob of Exceedance

    num = {}
    num['5'] = Calculate_number_of_members(GPM_ens, 5)
    num['10'] = Calculate_number_of_members(GPM_ens, 10)
    num['25'] = Calculate_number_of_members(GPM_ens, 25)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(11, 5), subplot_kw={'projection': ccrs.Miller()})
    rps = list(num.keys())
    for i in range(3):
        cm = Plot_number_of_members(num[rps[i]], axes[i], title=rps[i],
                                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    fig.subplots_adjust(bottom=0.13)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    cbar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Numer of Ensemble Members')
    cbar.set_ticks(np.arange(1, 19))
    plt.savefig('Plots/' + land + '/casestudy/GPM_num_members.png', bbox_inches='tight')


def Plot_exc_prob_ERA(land, event=None, lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        ERA_ens = Read_ensemble_GER('ERA')

    if land == 'UK':
        ERA_ens = Read_ensemble_UK('ERA', event)

    # Calculate Prob of Exceedance

    ex_prob = {}
    ex_prob['5'] = Calculate_exc_prob(ERA_ens, 5)
    ex_prob['10'] = Calculate_exc_prob(ERA_ens, 10)
    ex_prob['25'] = Calculate_exc_prob(ERA_ens, 25)
    ex_prob['50'] = Calculate_exc_prob(ERA_ens, 50)
    ex_prob['100'] = Calculate_exc_prob(ERA_ens, 100)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(11, 8), subplot_kw={'projection': ccrs.Miller()})
    rps = list(ex_prob.keys())
    for i in range(5):
        cm = Plot_exceedance_prob(ex_prob[rps[i]], axes.flatten()[i], title=rps[i],
                                  lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    axes[1, 2].remove()
    fig.subplots_adjust(bottom=0.18)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Probability >= RP [%]')
    plt.savefig('Plots/' + land + '/casestudy/ERA_ex_prob.png', bbox_inches='tight')


def Plot_num_memb_ERA(land, event=None,  lon_min=5, lon_max=15.5, lat_min=47, lat_max=55):
    # Read data
    if land == 'GER':
        ERA_ens = Read_ensemble_GER('ERA')

    if land == 'UK':
        ERA_ens = Read_ensemble_UK('ERA', event)

    # Calculate Prob of Exceedance

    num = {}
    num['5'] = Calculate_number_of_members(ERA_ens, 5)
    num['10'] = Calculate_number_of_members(ERA_ens, 10)
    num['25'] = Calculate_number_of_members(ERA_ens, 25)
    num['50'] = Calculate_number_of_members(ERA_ens, 50)
    num['100'] = Calculate_number_of_members(ERA_ens, 100)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(11, 8), subplot_kw={'projection': ccrs.Miller()})
    rps = list(num.keys())
    for i in range(5):
        cm = Plot_number_of_members(num[rps[i]], axes.flatten()[i], title=rps[i],
                                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

    axes[1, 2].remove()
    fig.subplots_adjust(bottom=0.18)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.05])
    cbar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Number of Ensemble members')
    cbar.set_ticks(np.arange(1, 19))
    plt.savefig('Plots/' + land + '/casestudy/ERA_num_members.png', bbox_inches='tight')

# %%
if __name__ == '__main__':
    Plot_exc_prob_GPM(land='GER')
    Plot_num_memb_GPM(land='GER')
    Plot_exc_prob_ERA(land='GER')
    Plot_num_memb_ERA(land='GER')


