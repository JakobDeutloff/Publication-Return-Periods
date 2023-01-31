import xarray as xr
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from cartopy import feature
import matplotlib.pyplot as plt
from Code.visualisation.Plot_prob_of_exceedance import Read_ensemble_GER, Read_ensemble_UK, Read_ensemble_ciara



def plot_max(land, event, lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53, box=None):
    if land == 'GER':
        Ens_GPM = Read_ensemble_GER('GPM', event)
        Ens_ERA = Read_ensemble_GER('ERA', event)
        transform = ccrs.PlateCarree()
        figsize = (5, 6)

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

        figsize = (5, 4)


    # calculate max rp
    max_ERA = Ens_ERA.max('ens_mem')
    max_GPM = Ens_GPM.max('ens_mem')
    dummy = xr.concat([max_ERA, max_GPM], dim='dataset')
    max_total = dummy.max('dataset')

    # plot
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.Miller()}, sharex='col', sharey='row')
    levels = [2.5, 7.5, 12.5, 37.5, 62.5, 137.5]
    cm = max_total.return_period.plot.contourf(ax=ax, cmap='YlOrRd', transform=transform, levels=levels,
                                               add_colorbar=False,
                                               vmin=5, vmax=100, extend='neither')
    ax.coastlines()
    ax.add_feature(feature.OCEAN, zorder=1, color='grey')
    ax.add_feature(feature.BORDERS)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # draw gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # add colorbar
    fig.subplots_adjust(bottom=0.2)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.04])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25, 50, 100])

    # Draw box
    if box:
        ax.add_patch(mpatches.Rectangle(xy=(box[0], box[-1]), width=box[1] - box[0], height=box[2] - box[3],
                                        transform=ccrs.PlateCarree(), fill=False, edgecolor='green'))

    plt.savefig('Plots/' + land + '/' + event + '/max_rp.png', bbox_inches='tight', dpi=500)

    plt.show()

if __name__ == '__main__':
    box_london_small = [0, 0.4, 51.5, 51.3]
    box_london = [-0.6, 0.4, 51.7, 51.3]
    box_cardif = [-3.3, -2.9, 51.75, 51.45]
    plot_max('UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53, box=box_london)
    # plot_max('UK', event='202002080300', lon_min=-6, lon_max=2, lat_max=50, lat_min=56, box=box_cardif)
    #plot_max('GER', event='20210711T1800Z', lon_min=5.7, lon_max=15, lat_min=47.3, lat_max=55)
