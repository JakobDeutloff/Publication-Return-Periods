#!/bin/bash
import iris
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

def write2nc(coords, max_value, outputfile):
    '''
    Write calculated max precipitation to a netCDF file
    coords:  ensmeble realisation, longitude and latitude
    max_value: Calculate cumulative maximum rainfall in a 24 hour window
    outputfile: File name for output
    '''

    # realization = coords[0]
    latitude =  coords[1]
    longitude = coords[2]
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
    
    max_precipitation = dataset_new.createVariable("max_precipitation", np.float64, ("latitude", "longitude"))
    max_precipitation.units = "m"
    max_precipitation[:] = max_value
    
    # Close the netcdf file
    dataset_new.close()

    return

def read_data(folder, i_date, precip_types, lead_time):
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
    # Corresponding starttime integer for the file
    # This value should be calculated explicity from the difference between forecast time and i_date 
    integ = 54
    # Setting up empty array
    total = ()
    for i in range(lead_time):
        print(start_date)
        rain = ()
        for precip in precip_types:
            year, mon, day, hour = strip(start_date)
            name = str(year) + str(mon).zfill(2) + str(day) + 'T' + str(hour).zfill(2) + trail  + str(integ).zfill(4) + extra + precip + '.nc'
            file = os.path.join(folder, name)            
            precipitation = iris.load_cube(file)
            coords = precipitation.coord('realization').points, precipitation.coord('latitude').points, precipitation.coord('longitude').points

            if len(rain) == 0:
                # Defining array dimension
                rain = np.zeros((len(coords[0]), len(coords[1]), len(coords[2])))
            rain = rain + precipitation.data
    
        if len(total) == 0:
             # Defining array dimension
            total = np.zeros((lead_time, len(coords[0]), len(coords[1]), len(coords[2])))

        total[i,:,:,:] = rain[:,:,:]
        start_date = step_datetime(start_date, date_format)
        integ += 1
    return coords, total

def loop_over_ensemble(coords, total, output_dir):
    '''
    Loop over the ensemble members. Calculate the maximum culimative precipitation
    over 24 hour window and then write to seperate .nc file
    '''
    num_ens = len(total[0,:,0,0])

    for i in range(num_ens):
        # Call to calculate the maximum precipitation
        max_precp = sliding_time(total[:,i,:,:])
        output_file = os.path.join(output_dir, "ens_member_%i.nc"%i)
        # Write to file
        write2nc(coords, max_precp, output_file)
    return

def main():
    filepath = "/Users/dangiles/Documents/MetOffice/DataChallenge/Forecast_Data/"
    flood = {"event": "Germany", "forecast" : "20210711T1800Z", "ens_folder": "mogrepsg"}
    folder = os.path.join(filepath, flood['ens_folder'], flood['forecast'])
    start_date = "2021071400"
    precip_type = ['rainfall_accumulation_from_convection-PT01H', 'rainfall_accumulation-PT01H', 'snowfall_accumulation_from_convection-PT01H', 'snowfall_accumulation-PT01H']
    lead_time = 48
    coords, total = read_data(folder, start_date, precip_type, lead_time)
    loop_over_ensemble(coords, total, folder)
    return


if __name__ == '__main__':
    main()

