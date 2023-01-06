import xarray as xr

# %%
path = 'Data/Postcode_Data_Gridded'
file1 = 'pcode_weighted_avg.tif'
file2 = 'pcode_weighted_sum.tif'

pcode_avg = xr.open_dataset(path + file1, engine='rasterio')
