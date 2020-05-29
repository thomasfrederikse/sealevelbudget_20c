# -------------------------
# Read GRACE/GFO mascons
# - Restore GIA
# - Remove seasonal cycle
# - Compute annual means
# - Save for quick import
# -------------------------

import numpy as np
import os
from netCDF4 import Dataset
import pyshtools as pysh

def main():
    set_settings()
    GRACE = read_grace_data()
    GRACE = remove_gia(GRACE)
    GRACE = remove_seasonal_cycle(GRACE)
    GRACE = compute_annual_mean(GRACE)
    save_data(GRACE)
    return

def set_settings():
    global settings
    settings = {}
    settings['remove_gia'] = False
    settings['dir_data']       = os.getenv('HOME') + '/Data/'
    settings['fn_grace_read']  = settings['dir_data'] + 'GRACE/JPL_mascon/GRCTellus.JPL.200204_201911.GLO.RL06M.MSCNv02CRI.nc'
    settings['fn_grace_write'] = settings['dir_data'] + 'Budget_20c/grd_prep/ewh_GRACE_annual.nc'
    settings['fn_mask']        = settings['dir_data'] + 'GRACE/JPL_mascon/mask.nc'
    settings['fn_mscn_coords'] = settings['dir_data'] + 'GRACE/JPL_mascon/mascon_coords.npy'
    settings['fn_love']        = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    settings['fn_ice6gd_stokes']= settings['dir_data'] + 'GIA/ICE6G_D/ICE6G_D_stokes_GRACE.txt'
    return

def read_grace_data():
    # Read original GRACE data
    print('   Reading GRACE data...')
    global settings
    GRACE = {}
    fh = Dataset(settings['fn_grace_read'],'r')
    fh.set_auto_mask(False)
    GRACE['time'] = 2002 + fh.variables['time'][:]/365.25
    GRACE['lat'] =fh.variables['lat'][:]
    GRACE['lon'] =fh.variables['lon'][:]
    GRACE['ewh'] = 10 * fh.variables['lwe_thickness'][:]
    GRACE['ewh_ste'] = 10 * fh.variables['uncertainty'][:]
    fh.close()

    # Mascon coordinates
    GRACE['mscn_coords'] = np.load(settings['fn_mscn_coords'],allow_pickle=True)
    return(GRACE)

def remove_gia(GRACE):
    # Restore the ICE6G_D correction to the GFO-data
    global settings
    if settings['remove_gia']:
        print('   Restore ICE6G_D GIA correction...')
        love = np.load(settings['fn_love'],allow_pickle=True).all()
        love['k'][0]=1
        love['k'][1]=1
        stokes,lmax  = pysh.shio.shread(settings['fn_ice6gd_stokes'])
        rho_earth  = 5517
        rad_earth  = 6.3710088e6
        gia_ice6gd_sh       = rad_earth*rho_earth*(2*np.arange(lmax+1)[np.newaxis, :, np.newaxis] + 1) / (3 * (1 + love['k'][np.newaxis, :lmax+1, np.newaxis])) * stokes
        gia_ice6gd_spat     = np.flipud(pysh.shtools.MakeGrid2D(gia_ice6gd_sh, interval=0.5, north=GRACE['lat'][-1], west=GRACE['lon'][0]))
        # Masconize
        gia_ice6gd_mscn = np.zeros(gia_ice6gd_spat.shape)
        for k in range(len(GRACE['mscn_coords'])):
            lat_acc = np.where((GRACE['lat'] >= GRACE['mscn_coords'][k, 0]) & (GRACE['lat'] < GRACE['mscn_coords'][k, 1]))[0]
            lon_acc = np.where((GRACE['lon'] >= GRACE['mscn_coords'][k, 2]) & (GRACE['lon'] < GRACE['mscn_coords'][k, 3]))[0]
            weight = np.cos(np.deg2rad(GRACE['lat'][lat_acc])) / np.mean(np.cos(np.deg2rad(GRACE['lat'][lat_acc])))  # Weight by cos lat
            gia_ice6gd_mscn[lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1] = np.nanmean(weight[:, np.newaxis] * gia_ice6gd_spat[lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1])
        # Restoring GIA:
        GRACE['ewh'] += gia_ice6gd_mscn[np.newaxis,...] * (GRACE['time']-GRACE['time'].mean())[:,np.newaxis,np.newaxis]
    return(GRACE)

def remove_seasonal_cycle(GRACE):
    # Remove seasonal cycle from field
    GRACE['ewh_noseas'] = np.zeros(GRACE['ewh'].shape)
    amat = np.ones([len(GRACE['time']),6])
    amat[:,1] = np.sin(2*np.pi*GRACE['time'])
    amat[:,2] = np.cos(2*np.pi*GRACE['time'])
    amat[:,3] = np.sin(4*np.pi*GRACE['time'])
    amat[:,4] = np.cos(4*np.pi*GRACE['time'])
    amat[:,5] = GRACE['time'] - GRACE['time'].mean()
    amat_T  = amat.T
    amat_sq = np.linalg.inv(np.dot(amat_T, amat))
    for i in range(GRACE['ewh'].shape[1]):
        for j in range(GRACE['ewh'].shape[2]):
            sol = np.dot(amat_sq, np.dot(amat_T, GRACE['ewh'][:, i, j]))
            sol[5] = 0
            GRACE['ewh'][:,i,j]-=np.matmul(amat,sol)
    return(GRACE)

def compute_annual_mean(GRACE):
    global settings
    GRACE['years'] = np.arange(2003,2020)
    GRACE['ewh_annual'] = np.zeros([len(GRACE['years']),len(GRACE['lat']),len(GRACE['lon'])])
    GRACE['ewh_ste_annual'] = np.zeros([len(GRACE['years']),len(GRACE['lat']),len(GRACE['lon'])])
    for idx, year in enumerate(GRACE['years']):
        acc_idx = (np.floor(GRACE['time']).astype(int) == year)
        GRACE['ewh_annual'][idx,...]     = GRACE['ewh'][acc_idx,...].mean(axis=0)
        GRACE['ewh_ste_annual'][idx,...] = np.sqrt(np.sum(GRACE['ewh_ste'][acc_idx, :, :] ** 2, axis=0)) / acc_idx.sum()
    return(GRACE)

def save_data(GRACE):
    global settings
    file_handle = Dataset(settings['fn_grace_write'], 'w')
    file_handle.createDimension('lon', len(GRACE['lon']))
    file_handle.createDimension('lat', len(GRACE['lat']))
    file_handle.createDimension('years', len(GRACE['years']))
    file_handle.createVariable('lon', 'f4', ('lon',),zlib=True)[:] = GRACE['lon']
    file_handle.createVariable('lat', 'f4', ('lat',),zlib=True)[:] = GRACE['lat']
    file_handle.createVariable('years', 'f4', ('years',),zlib=True)[:] = GRACE['years']
    file_handle.createVariable('ewh', 'i4', ('years', 'lat', 'lon',),zlib=True)[:] = GRACE['ewh_annual'] * 10
    file_handle.createVariable('ewh_ste', 'i4', ('years', 'lat', 'lon',),zlib=True)[:] = GRACE['ewh_ste_annual'] * 10
    file_handle.variables['ewh'].setncattr('scale_factor', 0.1)
    file_handle.variables['ewh_ste'].setncattr('scale_factor', 0.1)
    file_handle.close()
    return

if __name__ == '__main__':
    main()
