# ------------------------------------------------------
# Merge ERA20c and ERA5 by adjusting the mean over 1979.
# Average the results into annual means
# Compute wind stress from wind fields
# ------------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import mod_gentools as gentools
def main():
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['fn_ERA5']    = settings['dir_data'] + 'Reanalyses/ERA5/ERA5.nc'
    settings['fn_ERA20c']  = settings['dir_data'] + 'Reanalyses/ERA20c/ERA20c_1900_1980_sfc.nc'
    settings['fn_save']    = settings['dir_data'] + 'Budget_20c/tg/ERA.nc'
    settings['years'] = np.arange(1900,2019)

    ERA5   = read_dataset(settings['fn_ERA5'])
    ERA20c = read_dataset(settings['fn_ERA20c'])
    ERA    = merge_datasets(ERA5, ERA20c, settings)
    save_data(ERA,settings)
    return

def read_dataset(fn):
    Data = {}
    file_handle = Dataset(fn)
    file_handle.set_auto_mask(False)
    Data['time'] = 1900 + file_handle.variables['time'][:]/365/24
    Data['lon'] = file_handle.variables['longitude'][:]
    Data['lat'] = np.flipud(file_handle.variables['latitude'][:])
    Data['uwind'] = np.fliplr(file_handle.variables['u10'][:])
    Data['vwind'] = np.fliplr(file_handle.variables['v10'][:])
    Data['mslp'] = np.fliplr(file_handle.variables['msl'][:])
    file_handle.close()
    return(Data)

def merge_datasets(ERA5,ERA20c,settings):
    # overlap indices
    ERA5_ovl = (ERA5['time'] >=1979) & (ERA5['time'] < 1980)
    ERA20c_ovl = (ERA20c['time'] >=1979) & (ERA20c['time'] < 1980)

    # Remove bias in 1979
    ERA20c['mslp'] = ERA20c['mslp'] - ERA20c['mslp'][ERA20c_ovl,:,:].mean(axis=0)[np.newaxis,:,:] + ERA5['mslp'][ERA5_ovl,:,:].mean(axis=0)[np.newaxis,:,:]
    ERA20c['uwind'] = ERA20c['uwind'] - ERA20c['uwind'][ERA20c_ovl,:,:].mean(axis=0)[np.newaxis,:,:] + ERA5['uwind'][ERA5_ovl,:,:].mean(axis=0)[np.newaxis,:,:]
    ERA20c['vwind'] = ERA20c['vwind'] - ERA20c['vwind'][ERA20c_ovl,:,:].mean(axis=0)[np.newaxis,:,:] + ERA5['vwind'][ERA5_ovl,:,:].mean(axis=0)[np.newaxis,:,:]

    total_time = gentools.monthly_time(1900,2018)
    mslp  = np.vstack([ERA20c['mslp'][:-12,:,:],ERA5['mslp']])
    uwind = np.vstack([ERA20c['uwind'][:-12,:,:],ERA5['uwind']])
    vwind = np.vstack([ERA20c['vwind'][:-12,:,:],ERA5['vwind']])

    # Annual means from monthly means
    mslp_annual  = np.zeros([len(settings['years']),len(ERA5['lat']),len(ERA5['lon'])])
    uwind_annual = np.zeros([len(settings['years']),len(ERA5['lat']),len(ERA5['lon'])])
    vwind_annual = np.zeros([len(settings['years']),len(ERA5['lat']),len(ERA5['lon'])])
    for idx,yr in enumerate(settings['years']):
        yr_idx = (total_time>=yr)&(total_time<yr+1)
        mslp_annual[idx,:,:] = mslp[yr_idx,:,:].mean(axis=0)
        uwind_annual[idx,:,:] = uwind[yr_idx,:,:].mean(axis=0)
        vwind_annual[idx,:,:] = vwind[yr_idx,:,:].mean(axis=0)

    # From wind speed to wind stress
    uws_annual = (0.8 + 0.065 * np.sqrt(uwind_annual ** 2 + vwind_annual ** 2)) * uwind_annual * np.sqrt(uwind_annual ** 2 + vwind_annual ** 2)
    vws_annual = (0.8 + 0.065 * np.sqrt(uwind_annual ** 2 + vwind_annual ** 2)) * vwind_annual * np.sqrt(uwind_annual ** 2 + vwind_annual ** 2)

    ERA = {}
    ERA['lat']  = ERA5['lat']
    ERA['lon']  = ERA5['lon']
    ERA['time'] = settings['years']
    ERA['mslp'] = mslp_annual
    ERA['uws']  = uws_annual
    ERA['vws']  = vws_annual
    return(ERA)

def save_data(ERA,settings):
    file_handle = Dataset(settings['fn_save'], 'w')
    file_handle.createDimension('x', len(ERA['lon']))
    file_handle.createDimension('y', len(ERA['lat']))
    file_handle.createDimension('t', len(ERA['time']))
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = ERA['lon']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = ERA['lat']
    file_handle.createVariable('t', 'i4', ('t',),zlib=True)[:] = ERA['time']
    file_handle.createVariable('mslp', 'f4', ('t', 'y', 'x',),zlib=True)[:] = ERA['mslp']
    file_handle.createVariable('uws', 'f4', ('t', 'y', 'x',),zlib=True)[:] = ERA['uws']
    file_handle.createVariable('vws', 'f4', ('t', 'y', 'x',),zlib=True)[:] = ERA['vws']
    file_handle.close()
    return
