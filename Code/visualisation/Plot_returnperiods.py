#!/bin/bash
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np 
from scipy import interpolate 
from scipy.interpolate import NearestNDInterpolator
import cartopy.crs as ccrs
import cartopy.feature as ft

def read_ensemble_member(file, varname):
    '''
    Reads in the forecast ensemble netcdf files
    file: path to the ensemble forecast file 
    varname: variable name to read in (amount_of_precipitation)
    '''
    data = Dataset(file)[varname]
    
    x = Dataset(file)["longitude"]
    y = Dataset(file)["latitude"]
    data = np.array(data)
    data = data.astype(np.float32)
    return data, x, y

def regrid(data, lon, lat, x, y, order):
    f = interpolate.interp2d(x, y, data, kind=order)
    inter = f.__call__(lon, lat)
    return inter


def return_period_GPM(period, max_value, return_periods, times, x, y):
    '''
    Compares the maximum rainfall over the sliding 24 hour window against the rainfall thresholds
    for different return periods at each coordinate position
    
    period: zero array which will store the integer value of the return period threshold which has been met
    max_value: The calculated maximum rainfall over the sliding 24 hour windows in the forecast
    times: This can be the 5,10,25,50,100 Return period
    '''
    # plotting(max_value, lon, lat)
    for year in times:
        GPM = os.path.join(return_periods, str("GPM_GER_rp"+str(year)+".nc4"))
        var = "precipitationCal"
        lon, lat, return_estimate = read_GPM(GPM, var)
        if year == 5:
            max_value = regrid(max_value, lon, lat, x, y, 'cubic')
            period = np.zeros((len(lat),len(lon)))
        for i in range(max_value.shape[0]):
            for j in range(max_value.shape[1]):
                if return_estimate[i,j] < 10000.0:
                    if max_value[i,j] > return_estimate[i,j]:
                        period[i,j] = year
    return lon, lat, period

def return_period_ERA(period, max_value, return_periods, times, x, y):
    '''
    Compares the maximum rainfall over the sliding 24 hour window against the rainfall thresholds
    for different return periods at each coordinate position
    
    period: zero array which will store the integer value of the return period threshold which has been met
    max_value: The calculated maximum rainfall over the sliding 24 hour windows in the forecast
    times: This can be the 5,10,25,50,100 Return period
    '''
    for year in times:
        ERA = os.path.join(return_periods, str("ERA_GER_rp"+str(year)+".nc4"))
        var = "r"+str(year)+"yrrp"
        lon, lat, return_estimate = read_ERA(ERA, var)
        if year == 5:
            max_value = regrid(max_value, lon, lat, x, y, 'cubic')
            period = np.zeros((len(lat),len(lon)))
        for i in range(max_value.shape[0]):
            for j in range(max_value.shape[1]):
                if return_estimate[i,j] < 10000.0:
                    if max_value[i,j] > return_estimate[i,j]:
                        period[i,j] = year
    return lon, lat, period

def read_ERA(file, varname):
    '''
    Read in the ERA return netCDF file
    file: path to the file
    varname: variable to be read in ("r_%year_yrrp")
    '''
    v = Dataset(file)[varname]
    data = v[:,:]
    # Convert kg/m^2 to m
    # 1 kg/m^2 is equal to 0.001m
    data = data/1000
    x = Dataset(file)["longitude"]
    y = Dataset(file)["latitude"]
    x = np.array(x)
    y = np.array(y)
    x = x.astype(np.float32)
    y = y.astype(np.float32)
    
    return x, y, data

def read_GPM(file, varname):
    '''
    Read in the GPM return netCDF file
    file: path to the file
    varname: variable to be read in ("r_%year_yrrp")
    '''
    v = Dataset(file)[varname]
    data = v[:,:]
    # Convert kg/m^2 to m
    # 1 kg/m^2 is equal to 0.001m
    data = data/1000
    data = data[0,:,:].T
    x = Dataset(file)["lon"]
    y = Dataset(file)["lat"]
    x = np.array(x)
    y = np.array(y)
    x = x.astype(np.float32)
    y = y.astype(np.float32)
    
    return x, y, data

def plotting(data, x, y, i, type):

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    con = ax.contourf(x,y, data)
    ax.coastlines()
    ax.add_feature(ft.BORDERS)
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    fig.colorbar(con)
    plt.savefig('../../Plots/Ensemble_'+type+'_%i'%i, bbox_inches='tight')
    plt.show()

    return

def write2nc(coords, period, outputfile):
    '''
    Write calculated return period to a netCDF file
    coords:  longitude and latitude
    return_period: Calculate cumulative maximum rainfall in a 24 hour window
    outputfile: File name for output
    '''

    # realization = coords[0]
    latitude =  coords[1]
    longitude = coords[0]
    dataset_new = Dataset(outputfile, 'w', format='NETCDF4_CLASSIC')

    # ens = dataset_new.createDimension("realization", len(realization))
    lat = dataset_new.createDimension("latitude", len(latitude))
    lon = dataset_new.createDimension("longitude", len(longitude))

    # ens_members = dataset_new.createVariable("realization", np.float64, ("realization",))
    latitudes = dataset_new.createVariable("latitude", np.float64, ("latitude",))
    longitudes = dataset_new.createVariable("longitude", np.float64, ("longitude",))
    
    latitudes.units = 'degrees'
    longitudes.units = 'degrees'
    latitudes[:] = latitude
    longitudes[:] = longitude
    # ens_members[:] = realization
    
    return_period = dataset_new.createVariable("return_period", np.int32, ("latitude", "longitude"))
    return_period.units = "years"
    return_period[:] = period
    
    # Close the netcdf file
    dataset_new.close()

    return

def calculating_return_periods(file, returnpath, i):
    '''
    Calculate the spatially return priod for the ensemble forecasts
    Scenario: Ensemble member number
    
    '''
    print("Calculating ensemble member %i"%i)
    # Return periods
    ERA_times = [5,10,25,50,100]
    GPM_times = [5,10,25]
    # Set up some local paths and fixed variables
    varname = "max_precipitation"
    # Load in the ensemble member
    data, x, y = read_ensemble_member(file, varname)
    # print(np.amax(data), np.amin(data))

    # Plotting Ensemble Member
    # Set up dummy array for the return period array
    period_ERA = []
    period_GPM = []
    # # Calculate return period from ERA file
    type = 'ERA5'
    return_periods = os.path.join(returnpath, type, "GER")
    lon_ERA, lat_ERA, period_ERA = return_period_ERA(period_ERA, data, return_periods, ERA_times, x, y)
    # plotting(period_ERA, lon_ERA, lat_ERA, i, 'ERA')
    write2nc((lon_ERA, lat_ERA), period_ERA, '../../Plots/Ensemble_'+type+'_%i.nc'%i)

    type = 'GPM'
    return_periods = os.path.join(returnpath, "GPM", "GER")
    lon_GPM, lat_GPM, period_GPM = return_period_GPM(period_GPM, data, return_periods, GPM_times,  x, y)
    # plotting(period_GPM, lon_GPM, lat_GPM, i, 'GPM')
    write2nc((lon_GPM, lat_GPM), period_GPM, '../../Plots/Ensemble_'+type+'_%i.nc'%i)
    
    return period_ERA, period_GPM

def main():
    filepath = "/Users/dangiles/Documents/MetOffice/DataChallenge/Forecast_Data/"
    flood = {"event": "Germany", "forecast" : "20210711T1800Z", "ens_folder": "mogrepsg"}
    folder = os.path.join(filepath, flood['ens_folder'], flood['forecast'])

    returnpath = "/Users/dangiles/Documents/MetOffice/DataChallenge/ReturnPeriods/"

    for i in range(18):
        ens_mem = os.path.join(folder, "ens_member_%i.nc"%i)
        period_ERA, period_GPM = calculating_return_periods(ens_mem, returnpath, i)
        
    plt.show()

    return

if __name__ == '__main__':
    main()