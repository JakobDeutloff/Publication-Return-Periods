#!/bin/bash
import os
import iris
from iris.coords import DimCoord
from iris.cube import Cube
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np 


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
    
    return x, y, data


def plotting_xy(data, x, y, i, type):

    fig, ax = plt.subplots()
    con = ax.contourf(x,y, data)
    plt.gca().coastlines()
    plt.savefig('../Plots/Ensemble_'+type+'_%i'%i, bbox_inches='tight')
    plt.show()

    return

def return_period_GPM(period, max_value, return_periods, times, ens, outputfolder):
    '''
    Compares the maximum rainfall over the sliding 24 hour window against the rainfall thresholds
    for different return periods at each coordinate position
    
    period: zero array which will store the integer value of the return period threshold which has been met
    max_value: The calculated maximum rainfall over the sliding 24 hour windows in the forecast
    times: This can be the 5,10,25,50,100 Return period
    '''
    for year in times:
        # Read in the GPM file and regrid to the forecast domain

        GPM = os.path.join(return_periods, str("GPM_UK_rp"+str(year)+".nc4"))
        var = "precipitationCal"
        lon, lat, return_estimate = read_GPM(GPM, var)
        # Set up geographical cube 
        pp_coord_system = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        latitude = DimCoord(np.array(lat),standard_name='latitude',units='degrees', coord_system=pp_coord_system)
        longitude = DimCoord(np.array(lon), standard_name='longitude', units='degrees', coord_system=pp_coord_system)
        cube = Cube(return_estimate, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
        # Regrid the return period estimate
        regrid_return = cube.regrid(max_value, iris.analysis.Linear())

        if year == 5:
            period = np.zeros((regrid_return.shape[0], regrid_return.shape[1]))
        for i in range(regrid_return.shape[0]):
            for j in range(regrid_return.shape[1]):
                if max_value.data[i,j] > regrid_return.data[i,j]:
                    period[i,j] = year

    # Output a return period cube
    return_cube = Cube(period, dim_coords_and_dims=[(max_value.coord('projection_y_coordinate'), 0), (max_value.coord('projection_x_coordinate'), 1)], units = "year")
    output_file = os.path.join(outputfolder, 'Ensemble_GPM_%i.nc'%ens)
    # Save the file
    iris.save(return_cube, output_file)
    return 

def return_period_ERA(period, max_value, return_periods, times, ens, outputfolder):
    '''
    Compares the maximum rainfall over the sliding 24 hour window against the rainfall thresholds
    for different return periods at each coordinate position
    
    period: zero array which will store the integer value of the return period threshold which has been met
    max_value: The calculated maximum rainfall over the sliding 24 hour windows in the forecast
    times: This can be the 5,10,25,50,100 Return period
    '''
    for year in times:
        # Read in the ERA file and regrid to the forecast domain
        ERA = os.path.join(return_periods, str("ERA_UK_rp"+str(year)+".nc4"))
        var = "r"+str(year)+"yrrp"
        lon, lat, return_estimate = read_ERA(ERA, var)
        # Set up geographical cube 
        pp_coord_system = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        latitude = DimCoord(np.array(lat),standard_name='latitude',units='degrees', coord_system=pp_coord_system)
        longitude = DimCoord(np.array(lon), standard_name='longitude', units='degrees', coord_system=pp_coord_system)
        cube = Cube(return_estimate, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
        # Regrid the return period estimate
        regrid_return = cube.regrid(max_value, iris.analysis.Linear())

        if year == 5:
            period = np.zeros((regrid_return.shape[0], regrid_return.shape[1]))
        for i in range(regrid_return.shape[0]):
            for j in range(regrid_return.shape[1]):
                if max_value.data[i,j] > regrid_return.data[i,j]:
                    period[i,j] = year

    # Output a return period cube
    return_cube = Cube(period, dim_coords_and_dims=[(max_value.coord('projection_y_coordinate'), 0), (max_value.coord('projection_x_coordinate'), 1)], units = "year")
    output_file = os.path.join(outputfolder, 'Ensemble_ERA_%i.nc'%ens)
    # Save the file
    iris.save(return_cube, output_file)
    return 

def calculating_return_periods(file, returnpath, flood, i):
    '''
    Calculate the spatially return priod for the ensemble forecasts
    Scenario: Ensemble member number
    
    '''
    print("Calculating ensemble member %i"%i)
    # Return periods
    ERA_times = [5,10,25,50,100]
    GPM_times = [5,10,25]
    # Load in the ensemble member
    data = iris.load_cube(file)
    # Plotting Ensemble Member
    # Set up dummy array for the return period array
    period_ERA = []
    period_GPM = []
    # # Calculate return period from ERA file
    type = 'ERA'
    return_periods = os.path.join(returnpath, type, "UK")
    outputfolder = os.path.join(return_periods, flood['forecast'])
    return_period_ERA(period_ERA, data, return_periods, ERA_times, i, outputfolder)
    print("ERA5 done")
    type = 'GPM'
    return_periods = os.path.join(returnpath, type, "UK")
    outputfolder = os.path.join(return_periods, flood['forecast'])
    return_period_GPM(period_GPM, data, return_periods, GPM_times, i, outputfolder)
    print("GPM done \n")

    return

def main():
    filepath = ""
    # flood = {"event": "UK", "ens_folder": "mogrepsuk", "forecast" : "20210723T2100Z"}
    flood = {"event": "UK", "ens_folder": "mogrepsuk", "forecast" : "20210709T2100Z"}

    folder = os.path.join(filepath, flood['ens_folder'], flood['forecast'])
    n_ens = 18
    returnpath = ""

    for i in range(n_ens):
        ens_mem = os.path.join(folder, "ens_member_%i.nc"%i)
        calculating_return_periods(ens_mem, returnpath, flood, i)
    return

if __name__ == '__main__':
    main()