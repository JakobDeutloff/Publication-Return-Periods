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
    # Corresponding start time integer for the file
    # This value should be calculated explicity from the difference between forecast time and i_date 
    integ = 1
    # Setting up empty array
    total = ()
    for i in range(lead_time):
        print(start_date)
        rain = ()
        for precip in precip_types:
            year, mon, day, hour = strip(start_date)
            name = str(year) + str(mon).zfill(2) + str(day).zfill(2) + 'T' + str(hour).zfill(2) + trail  + str(integ).zfill(4) + extra + precip + '.nc'
            file = os.path.join(folder, name)            
            precipitation = iris.load_cube(file)
            coords = precipitation.coord('realization'), precipitation.coord('projection_y_coordinate'), precipitation.coord('projection_x_coordinate')

            if len(rain) == 0:
                # Defining array dimension
                rain = np.zeros((len(coords[0].points), len(coords[1].points), len(coords[2].points)))
            rain = rain + precipitation.data
    
        if len(total) == 0:
             # Defining array dimension
            total = np.zeros((lead_time, len(coords[0].points), len(coords[1].points), len(coords[2].points)))

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
        write2nc_cube(coords, max_precp, output_file)
    return

def main():
    filepath = ""
    # flood = {"event": "UK", "ens_folder": "mogrepsuk", "forecast" : "20210723T2100Z"}
    flood = {"event": "UK", "ens_folder": "mogrepsuk", "forecast" : "20210709T2100Z"}
    folder = os.path.join(filepath, flood['ens_folder'], flood['forecast'])
    start_date = "2021070922"
    # start_date = "2021072322"
    precip_type = ['rainfall_accumulation-PT01H', 'snowfall_accumulation-PT01H']
    lead_time = 60
    coords, total = read_data(folder, start_date, precip_type, lead_time)
    loop_over_ensemble(coords, total, folder)
    return


if __name__ == '__main__':
    main()

