#!/bin/bash
import iris

from iris.coords import DimCoord
from iris.cube import Cube
import os
from netCDF4 import Dataset
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from datetime import datetime, date, timedelta
import numpy as np 

def step_datetime(idate, date_format):
    '''
    Move the dateime forward in time
    '''
    delta = timedelta(hours=1)
    new_idate = datetime.strptime(idate, date_format) + delta
    return new_idate.strftime(date_format)

def strip(date):
    '''
    Extract the year, month, day and hour from the DateTime string
    '''
    date_format = "%Y%m%d%H"
    year = datetime.strptime(date, date_format).year
    mon = datetime.strptime(date, date_format).month
    day = datetime.strptime(date, date_format).day
    hour = datetime.strptime(date, date_format).hour
    return year, mon, day, hour

def sliding_time(data):
    '''
    Calculates the maximum rainfal over sliding 24 hour windows of the
    48 hour lead time forecast
    
    data: precipitation array at all the time steps
    '''
    max_value = data[0,:,:]
    summed_value = data[0,:,:]
    for t in range((len(data[:,0,0])-24)):
        summed_value = np.sum(data[t:t+24,:,:], axis=0)
        max_value = np.maximum(summed_value[:,:], max_value[:,:])

    return max_value

def write2nc_cube(coords, max_value, outputfile):
    '''
    Write calculated max precipitation to a netCDF file
    coords:  ensmeble realisation, longitude and latitude
    max_value: Calculate cumulative maximum rainfall in a 24 hour window
    outputfile: File name for output
    '''
    y = coords[1]
    x = coords[2]
    cube = Cube(max_value, dim_coords_and_dims=[(y, 0), (x, 1)], standard_name="thickness_of_rainfall_amount", units = "m")
    iris.save(cube, outputfile)
    return

def read_data(folder, i_date, precip_types, lead_time, ens_type):
    '''
    Read the data from the four different precipitation files.
    Combine the precipitation and output to an array.
    folder: Path to directory with all the forecast files
    i_date: Start date at which the event occurs
    precip_types: The four types of precipitation files
    lead_time: The hourly length of the forecast window
    '''
    # Some extra strings to match formatting of the UM files
    trail = '00Z-PT'
    extra = 'H00M-'
    # Date Time format
    date_format = "%Y%m%d%H"
    start_date = i_date
    # Corresponding start time integer for the file
    # This value should be calculated explicity from the difference between forecast time and i_date 
    integ = 1
    # Setting up empty array
    total = ()
    for i in range(int(lead_time)):
        print(start_date)
        rain = ()
        for precip in precip_types:
            year, mon, day, hour = strip(start_date)
            name = str(year) + str(mon).zfill(2) + str(day).zfill(2) + 'T' + str(hour).zfill(2) + trail  + str(integ).zfill(4) + extra + precip + '.nc'
            file = os.path.join(folder, name)            
            precipitation = iris.load_cube(file)
            if ens_type == "mogrepsuk":
                coords = precipitation.coord('realization'), precipitation.coord('projection_y_coordinate'), precipitation.coord('projection_x_coordinate')
            else: 
                GER_lat = iris.Constraint(latitude=lambda v: v > 47 and v <= 55 )
                GER_lon = iris.Constraint(longitude=lambda v: v > 5.04 and v <= 15.5 )
                precipitation =  precipitation.extract(GER_lat & GER_lon)
                coords = precipitation.coord('realization'), precipitation.coord('latitude'), precipitation.coord('longitude')

            if len(rain) == 0:
                # Defining array dimension
                rain = np.zeros((len(coords[0].points), len(coords[1].points), len(coords[2].points)))
            rain = rain + precipitation.data
    
        if len(total) == 0:
             # Defining array dimension
            total = np.zeros((int(lead_time), len(coords[0].points), len(coords[1].points), len(coords[2].points)))

        total[i,:,:,:] = rain[:,:,:]
        start_date = step_datetime(start_date, date_format)
        integ += 1
    return coords, total

def loop_over_ensemble(coords, total, output_dir, returnpath, flood):
    '''
    Loop over the ensemble members. Calculate the maximum culimative precipitation
    over 24 hour window and then write to seperate .nc file
    '''
    num_ens = len(total[0,:,0,0])

    for i in range(num_ens):
        # Call to calculate the maximum precipitation
        max_precp = sliding_time(total[:,i,:,:])
        output_file = os.path.join(output_dir, "ens_member_%i.nc"%i)
        # Write to max cumulative precip to file
        write2nc_cube(coords, max_precp, output_file)
        # Calculate the return period
        calculating_return_periods(output_file, returnpath, flood, i)
    return


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


def return_period_GPM(period, cum_precip, return_periods, times, ens, flood, outputfolder):
    '''
    Compares the maximum rainfall over the sliding 24 hour window against the rainfall thresholds
    for different return periods at each coordinate position
    
    period: zero array which will store the integer value of the return period threshold which has been met
    max_value: The calculated maximum rainfall over the sliding 24 hour windows in the forecast
    times: This can be the 5,10,25,50,100 Return period
    '''
    for year in times:
        # Read in the GPM file and regrid to the forecast domain

        GPM = os.path.join(return_periods, str("GPM_"+str(flood['event'])+"_rp"+str(year)+".nc4"))
        var = "precipitationCal"
        lon, lat, return_estimate = read_GPM(GPM, var)
        # Set up geographical cube 
        pp_coord_system = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        latitude = DimCoord(np.array(lat),standard_name='latitude',units='degrees', coord_system=pp_coord_system)
        longitude = DimCoord(np.array(lon), standard_name='longitude', units='degrees', coord_system=pp_coord_system)
        return_from_file = Cube(return_estimate, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])

        # Regrid the return period estimate if UK

        regrid_return = return_from_file.regrid(cum_precip, iris.analysis.Nearest())

        if year == 5:
            period = np.zeros((regrid_return.shape[0], regrid_return.shape[1]))
        for i in range(regrid_return.shape[0]):
            for j in range(regrid_return.shape[1]):
                if cum_precip.data[i,j] > regrid_return.data[i,j]:
                    period[i,j] = year

    output_file = os.path.join(outputfolder, 'Ensemble_ReturnPeriod_GPM_%i.nc'%ens)
    # Output a return period cube
    if flood['event'] == "UK":
        iris.save(Cube(period, dim_coords_and_dims=[(cum_precip.coord('projection_y_coordinate'), 0), (cum_precip.coord('projection_x_coordinate'), 1)], units = "year"), output_file)
    else:
        iris.save(Cube(period, dim_coords_and_dims=[(cum_precip.coord('latitude'), 0), (cum_precip.coord('longitude'), 1)], units = "year"), output_file)
    return 

def return_period_ERA(period, cum_precip, return_periods, times, ens, flood,  outputfolder):
    '''
    Compares the maximum rainfall over the sliding 24 hour window against the rainfall thresholds
    for different return periods at each coordinate position
    
    period: zero array which will store the integer value of the return period threshold which has been met
    max_value: The calculated maximum rainfall over the sliding 24 hour windows in the forecast
    times: This can be the 5,10,25,50,100 Return period
    '''
    for year in times:
        # Read in the ERA file and regrid to the forecast domain
        ERA = os.path.join(return_periods, str("ERA_"+str(flood['event'])+"_rp"+str(year)+".nc4"))
        var = "r"+str(year)+"yrrp"
        lon, lat, return_estimate = read_ERA(ERA, var)
        # Set up geographical cube 
        pp_coord_system = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        latitude = DimCoord(np.array(lat),standard_name='latitude',units='degrees', coord_system=pp_coord_system)
        longitude = DimCoord(np.array(lon), standard_name='longitude', units='degrees', coord_system=pp_coord_system)
        return_from_file = Cube(return_estimate, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
        # Regrid the return period estimate
        regrid_return = return_from_file.regrid(cum_precip, iris.analysis.Nearest())
        

        if year == 5:
            period = np.zeros((regrid_return.shape[0], regrid_return.shape[1]))
        for i in range(regrid_return.shape[0]):
            for j in range(regrid_return.shape[1]):
                if cum_precip.data[i,j] > regrid_return.data[i,j]:
                    period[i,j] = year

    output_file = os.path.join(outputfolder, 'Ensemble_ReturnPeriod_ERA_%i.nc'%ens)
    # Output a return period cube
    if flood['event'] == "UK":
        iris.save(Cube(period, dim_coords_and_dims=[(cum_precip.coord('projection_y_coordinate'), 0), (cum_precip.coord('projection_x_coordinate'), 1)], units = "year"), output_file)
    else:
        iris.save(Cube(period, dim_coords_and_dims=[(cum_precip.coord('latitude'), 0), (cum_precip.coord('longitude'), 1)], units = "year"), output_file)
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
    return_periods = os.path.join(returnpath, type, flood['event'])
    outputfolder = os.path.join(return_periods, flood['forecast'])
    return_period_ERA(period_ERA, data, return_periods, ERA_times, i, flood, outputfolder)
    print("ERA5 done")
    type = 'GPM'
    return_periods = os.path.join(returnpath, type, flood['event'])
    outputfolder = os.path.join(return_periods, flood['forecast'])
    return_period_GPM(period_GPM, data, return_periods, GPM_times, i, flood, outputfolder)
    print("GPM done \n")

    return

def main():
    # Define users paths
    filepath = ""
    returnpath = ""

    # Defining the event and its folder
    flood = {"event": "GER", "ens_folder": "mogrepsg", "forecast" : "20210711T1800Z", "start_date" : "2021071119", "lead_time" : "132", "precip_type" : ["rainfall_accumulation_from_convection-PT01H", "rainfall_accumulation-PT01H", "snowfall_accumulation_from_convection-PT01H", "snowfall_accumulation-PT01H"]}
    # flood = {"event": "UK", "ens_folder": "mogrepsuk", "forecast" : "20210723T2100Z", "start_date" : "2021072322", "lead_time" : "60", "precip_type" : ["rainfall_accumulation-PT01H", "snowfall_accumulation-PT01H"]}
    # flood = {"event": "UK", "ens_folder": "mogrepsuk", "forecast" : "20210709T2100Z", "start_date" : "2021070922", "lead_time" : "60", "precip_type" : ["rainfall_accumulation-PT01H", "snowfall_accumulation-PT01H"]}
    folder = os.path.join(filepath, flood['ens_folder'], flood['forecast'])
    coords, total = read_data(folder, flood["start_date"], flood["precip_type"], flood["lead_time"], flood['ens_folder'])
    loop_over_ensemble(coords, total, folder, returnpath, flood)
    return


if __name__ == '__main__':
    main()

