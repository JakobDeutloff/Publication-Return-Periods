import xarray as xr
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from cartopy import feature
import numpy as np
import matplotlib.pyplot as plt
from Code.visualisation.Plot_prob_of_exceedance import Read_ensemble_GER, Read_ensemble_UK, Read_ensemble_ciara


def get_box(land, event, box):
    if land == 'GER':
        Ens_GPM = Read_ensemble_GER('GPM', event)
        Ens_ERA = Read_ensemble_GER('ERA', event)
        transform = ccrs.PlateCarree()

    if land == 'UK':
        if event == '202002080300':
            Ens_GPM = Read_ensemble_ciara('GPM', event)
            Ens_ERA = Read_ensemble_ciara('ERA', event)
            merc = Ens_ERA.transverse_mercator
            transform = ccrs.LambertAzimuthalEqualArea(
                central_latitude=merc.latitude_of_projection_origin,
                central_longitude=merc.longitude_of_central_meridian,
                false_easting=merc.false_easting,
                false_northing=merc.false_northing,
            )
        else:
            Ens_GPM = Read_ensemble_UK('GPM', event)
            Ens_ERA = Read_ensemble_UK('ERA', event)
            central_latitude = Ens_GPM.lambert_azimuthal_equal_area.latitude_of_projection_origin
            central_longitude = Ens_GPM.lambert_azimuthal_equal_area.longitude_of_projection_origin
            transform = ccrs.LambertAzimuthalEqualArea(
                central_latitude=central_latitude,
                central_longitude=central_longitude)

    box_trans = [None for i in range(4)]
    box_trans[0], box_trans[2] = transform.transform_point(x=box[0], y=box[2], src_crs=ccrs.PlateCarree())
    box_trans[1], box_trans[3] = transform.transform_point(x=box[1], y=box[3], src_crs=ccrs.PlateCarree())

    ERA_box = Ens_ERA.sel(projection_x_coordinate=slice(box_trans[0], box_trans[1]),
                          projection_y_coordinate=slice(box_trans[3], box_trans[2]))

    GPM_box = Ens_GPM.sel(projection_x_coordinate=slice(box_trans[0], box_trans[1]),
                          projection_y_coordinate=slice(box_trans[3], box_trans[2]))

    return ERA_box, GPM_box


def plot_hist(land, event, ERA_box, GPM_box):

    fig, axes = plt.subplots(2, 1)

    bins = [0, 2.5, 7.5, 12.5, 37.5, 75, 125]
    ticks = range(len(bins)-1)
    h, e = np.histogram(ERA_box.return_period.values.flatten(), bins=bins)
    axes[0].bar(ticks, h, width=1, edgecolor='k')
    axes[0].xaxis.set_ticks(ticks)
    axes[0].xaxis.set_ticklabels(['0', '5', '10', '25', '50', '100'])
    axes[0].set_title('ERA')

    bins = [0, 2.5, 7.5, 12.5, 37.5]
    ticks = range(len(bins)-1)
    h, e = np.histogram(GPM_box.return_period.values.flatten(), bins=bins)
    axes[1].bar(ticks, h, width=1, edgecolor='k')
    axes[1].xaxis.set_ticks(ticks)
    axes[1].xaxis.set_ticklabels(['0', '5', '10', '25'])
    axes[1].set_title('GPM')

    plt.savefig('Plots/' + land + '/' + event + '/hist.png', bbox_inches='tight')

    plt.show()

def count_occurences(box,RPs):

    occurences = np.zeros((box.ens_mem.size, len(RPs)))

    for i in range(occurences.shape[0]):
        for j in range(occurences.shape[1]):
            bool_oc = box.sel(ens_mem=i).return_period.values.flatten() == RPs[j]
            occurences[i, j] = np.sum(bool_oc)

    return occurences



def plot_wiskers(land, event, ERA_oc, GPM_oc, zeros):

    # Sort data so that equal RPs are next to each other
    data = np.zeros((ERA_oc.shape[0], ERA_oc.shape[1] + GPM_oc.shape[1]))
    for i in np.arange(0, GPM_oc.shape[1]*2, 2):
        data[:, i] = GPM_oc[:, int(i/2)]
        data[:, i+1] = ERA_oc[:, int(i/2)]
    data[:, -1] = ERA_oc[:, -1]
    data[:, -2] = ERA_oc[:, -2]

    # Normalize data with respect to grid cells
    data_norm = data/np.sum(data[0, :]) * 100

    # Plot
    fig, ax = plt.subplots()
    if zeros:
        bp = ax.boxplot(data_norm, patch_artist=True)
        # Set xticks
        labels = ax.set_xticklabels(['0', '0', '5', '5', '10', '10', '25', '25', '50', '100'])
    else:
        bp = ax.boxplot(data_norm[:, 2:], patch_artist=True)
        # Set xticks
        labels = ax.set_xticklabels(['5', '5', '10', '10', '25', '25', '50', '100'])

    # Set color of boxes
    colors = {'ERA': 'lightgreen', 'GPM': 'lightblue'}
    colors_lab = {'ERA': 'green', 'GPM': 'blue'}
    keys = ['boxes']
    for key in keys:
        for i in np.arange(0, GPM_oc.shape[1] * 2, 2):
            labels[i].set_color(colors_lab['GPM'])
            bp[key][i].set_facecolor(colors['GPM'])
            bp[key][i+1].set_facecolor(colors['ERA'])
            labels[i+1].set_color(colors_lab['ERA'])
        labels[-1].set_color(colors_lab['ERA'])
        labels[-2].set_color(colors_lab['ERA'])
        bp[key][-1].set_facecolor(colors['ERA'])
        bp[key][-2].set_facecolor(colors['ERA'])

    # Labels
    ax.set_xlabel('Return Period [Yr]')
    ax.set_ylabel('Fraction of Gridcells [%]')

    plt.savefig('Plots/' + land + '/' + event + '/boxplot.png', bbox_inches='tight', dpi=500)
    plt.show()
    return data_norm

# %%
box_london = [-0.6, 0.4, 51.7, 51.3]
box_london_small = [-0.3, 0.2, 51.6, 51.4]
box_cardif = [-3.3, -2.9, 51.75, 51.45]
ERA_cardif, GPM_cardif = get_box(land='UK', event='202002080300', box=box_cardif)
ERA_london, GPM_london = get_box(land='UK', event='20210723T2100Z', box=box_london)
# %%
# plot_hist(land='UK', event='202002080300', ERA_box=ERA_box, GPM_box=GPM_box)

# %%
oc_GPM_cardif = count_occurences(GPM_cardif, RPs=[0, 5, 10, 25])
oc_ERA_cardif = count_occurences(ERA_cardif, RPs=[0, 5, 10, 25, 50, 100])
oc_GPM_london = count_occurences(GPM_london, RPs=[0, 5, 10, 25])
oc_ERA_london = count_occurences(ERA_london, RPs=[0, 5, 10, 25, 50, 100])


# %%
data_cardif = plot_wiskers(land='UK', event='202002080300', ERA_oc=oc_ERA_cardif, GPM_oc=oc_GPM_cardif, zeros=False)
data_london = plot_wiskers(land='UK', event='20210723T2100Z', ERA_oc=oc_ERA_london, GPM_oc=oc_GPM_london,
                           zeros=False)





