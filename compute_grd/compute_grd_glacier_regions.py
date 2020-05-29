# -------------------------------------------
# Compute fingerprint for each glacier region
# -------------------------------------------
from netCDF4 import Dataset
import numpy as np
import os
import pySLE

def main():
    print('Glacier regions GRD...')
    settings = {}
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    settings['fn_mask'] = settings['dir_data'] + 'Budget_20c/grd_prep/mask.npy'
    settings['fn_love'] = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    settings['fn_grd_glacier_regions'] = settings['dir_data'] + 'Budget_20c/grd_prep/grd_glacier_regions.nc'

    love = np.load(settings['fn_love'],allow_pickle=True).all()
    mask = np.load(settings['fn_mask'], allow_pickle=True).all()

    rsl_glacier_regions = np.zeros([mask['glacier_num'].shape[0],mask['lat'].shape[0],mask['lon'].shape[0]])
    rad_glacier_regions = np.zeros([mask['glacier_num'].shape[0],mask['lat'].shape[0],mask['lon'].shape[0]])

    for region in range(mask['glacier_num'].shape[0]):
        sle = pySLE.solver(lat=mask['lat'], lon=mask['lon'], load=mask['glacier_scale'][region,:,:], slm=1.0-mask['land'], love=love, lmax=359)
        sle.solve()
        rsl_glacier_regions[region,:,:] = sle.result['rsl']/sle.result['barystatic']
        rad_glacier_regions[region,:,:] = sle.result['rad']/sle.result['barystatic']

    # Save data
    file_handle = Dataset(settings['fn_grd_glacier_regions'], 'w')
    file_handle.createDimension('x', len(sle.result['lon']))
    file_handle.createDimension('y', len(sle.result['lat']))
    file_handle.createDimension('rgi', mask['glacier_num'].shape[0])
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = sle.result['lon']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = sle.result['lat']
    file_handle.createVariable('rgi', 'i4', ('rgi',),zlib=True)[:] = mask['glacier_num']
    file_handle.createVariable('rsl', 'f4', ('rgi', 'y', 'x',),zlib=True,complevel=4,least_significant_digit=3)[:] = rsl_glacier_regions
    file_handle.createVariable('rad', 'f4', ('rgi', 'y', 'x',),zlib=True,complevel=4,least_significant_digit=3)[:] = rad_glacier_regions
    file_handle.close()
    return

if __name__ == '__main__':
    main()
