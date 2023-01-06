import xarray as xr
import cartopy.crs as ccrs
from cartopy import feature
import matplotlib.pyplot as plt
from Code.visualisation.Plot_prob_of_exceedance import Read_ensemble_GER, Read_ensemble_UK


def plot_max(land, event, lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53):
    if land == 'GER':
        Ens_GPM = Read_ensemble_GER('GPM', event)
        Ens_ERA = Read_ensemble_GER('ERA', event)
        transform = ccrs.PlateCarree()
        figsize = (5, 6)

    if land == 'UK':
        Ens_GPM = Read_ensemble_UK('GPM', event)
        Ens_ERA = Read_ensemble_UK('ERA', event)
        central_latitude = Ens_GPM.lambert_azimuthal_equal_area.latitude_of_projection_origin
        central_longitude = Ens_GPM.lambert_azimuthal_equal_area.longitude_of_projection_origin
        transform = ccrs.LambertAzimuthalEqualArea(
            central_latitude=central_latitude,
            central_longitude=central_longitude)
        figsize = (5, 5)

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
    ax.set_title('Maximum RP  ' + event)

    # draw gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # add colorbar
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.04])
    c_bar = fig.colorbar(cm, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Return Period')
    c_bar.set_ticks([5, 10, 25, 50, 100])

    plt.savefig('Plots/' + land + '/casestudy/max_rp' + event + '.png', bbox_inches='tight')

    plt.show()


if __name__ == '__main__':
    plot_max('UK', event='20210723T2100Z', lon_min=-3, lon_max=2, lat_max=50.5, lat_min=53)
    #plot_max('GER', event='20210711T1800Z', lon_min=5, lon_max=15.5, lat_min=47, lat_max=55)
