# --------------------------------------------------------
# Regrid altimetry into annual data on 0.5x0.5 degree grid
# Read NASA MEASURES gridded altimetry and average
# in space and time
# --------------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os


def main():
    settings = {}
    settings['dir_data']     = os.getenv('HOME') + '/Data/'
    settings['dir_measures'] = settings['dir_data'] + 'Altimetry/MEASURES/'
    settings['years']        = np.arange(1993,2019)
    settings['dir_save']     = settings['dir_data'] + '/Budget_20c/vlm/Altimetry_annual.nc'

    lat,lon,ssh_annual_hires = annual_average(settings)
    lat_new,lon_new,ssh_annual_regrid = reduce_resolution(lat, lon, ssh_annual_hires, settings)
    save_data(lat_new, lon_new, ssh_annual_regrid, settings)
    return

def annual_average(settings):
    # ------------------------------------------
    # Average high-frequency altimetry data into
    # annual-mean values
    # ------------------------------------------

    # Make list of all files
    dirlist = os.listdir(settings['dir_measures'])
    dirlist = [s for s in dirlist if "ssh_grids_v1812_" in s]
    file_array = np.zeros([len(dirlist), 3], dtype=object)
    for n in range(len(dirlist)):
        file_array[n, 0] = dirlist[n]
        file_array[n, 1] = np.int(dirlist[n][16:20])
        file_array[n, 2] = np.int(dirlist[n][20:22])

    # Read MEASURES lat/lon
    lon = Dataset(settings['dir_measures'] + file_array[0, 0], 'r').variables['Longitude'][:]._get_data()
    lat = Dataset(settings['dir_measures'] + file_array[0, 0], 'r').variables['Latitude'][:]._get_data()

    # Towards annual values
    ssh_annual_hires = np.zeros([len(settings['years']),len(lat),len(lon)])
    for idx,yr in enumerate(settings['years']):
        print(yr)
        # Find all grids within a year and compute the averaged
        file_idx = (file_array[:,1] == yr)
        ssh_in_year = np.zeros([file_idx.sum(),len(lat),len(lon)])
        file_years = file_array[file_idx, 0]
        for num,fname in enumerate(file_years):
            read_data = Dataset(settings['dir_measures'] + fname, 'r').variables['SLA'][0, ...]._get_data().T
            read_data[read_data>1e10] = np.nan
            ssh_in_year[num,:,:] = 1000 * read_data
        mask = np.sum(np.isfinite(ssh_in_year),axis=0)>0.75*file_idx.sum()
        ssh_annual_hires[idx,:,:] = mask*np.nanmean(ssh_in_year,axis=0)
    return(lat,lon,ssh_annual_hires)

def reduce_resolution(lat,lon,ssh_annual_hires,settings):
    # -------------------------------------------
    # Reduce the spatial resolution to 0.5 by 0.5
    # -------------------------------------------
    lat_new = np.arange(-89.75,90.25,0.5)
    lon_new = np.arange(0.25,360.25,0.5)
    ssh_annual_regrid = np.zeros([len(settings['years']),len(lat_new),len(lon_new)])*np.nan
    for i in range(len(lat_new)):
        print(i)
        for j in range(len(lon_new)):
            idx_i = (np.abs((lat - lat_new[i]))) < 0.2
            idx_j = (np.abs((lon - lon_new[j]))) < 0.2
            if (idx_i.sum()>0) & (idx_j.sum()>0):
                idx_i = np.where(idx_i)[0]
                idx_j = np.where(idx_j)[0]
                ssh_sampled = ssh_annual_hires[:,idx_i[0]:idx_i[-1]+1,idx_j[0]:idx_j[-1]+1]
                if np.isfinite(ssh_sampled[-1, :, :]).sum()>2:
                    ssh_annual_regrid[:,i,j] = np.nanmean(ssh_sampled,axis=(1,2))

    # Only retain grid cells for which all years are available
    acc_mask = np.sum(np.isfinite(ssh_annual_regrid),axis=0)==len(settings['years'])
    ssh_annual_regrid = ssh_annual_regrid * acc_mask[np.newaxis,:,:]
    return(lat_new,lon_new,ssh_annual_regrid)

def save_data(lat_new,lon_new,ssh_annual_regrid,settings):
    # Save data into netCDF file
    file_handle = Dataset(settings['dir_save'], 'w')
    file_handle.createDimension('t', len(settings['years']))
    file_handle.createDimension('x', len(lon_new))
    file_handle.createDimension('y', len(lat_new))
    file_handle.createVariable('t', 'f4', ('t',),zlib=True)[:] = settings['years']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = lat_new
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = lon_new
    file_handle.createVariable('z', 'f4', ('t', 'y', 'x'),zlib=True)[:] = ssh_annual_regrid
    file_handle.close()
    return
