# ---------------------------------------------
# Compute and save the GRD fingerprints for:
# Dam retention (1xtime)
# Groundwater depletion (1xtime)
# Natural TWS changes (100xtime)
# And svae results for ensemble GRD generation
# ---------------------------------------------
from netCDF4 import Dataset
import numpy as np
import os
import pySLE
import multiprocessing as mp
import datetime as dt

def main():
    print('TWS GRD...')
    global settings
    settings = {}
    settings['time']    = np.arange(1900,2004)
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 2
        settings['dir_data']    = os.getenv('HOME') + '/Data/'
        settings['dir_scratch'] = os.getenv('HOME') + '/Storage/'
    else:
        settings['nproc'] = 32
        settings['dir_data']    = os.getenv('HOME') + '/Data/'
        settings['dir_scratch'] = os.getenv('HOME') + '/Storage/'
    # Directories
    settings['dir_gwd_prep'] = settings['dir_data'] + 'Budget_20c/grd_prep/'
    settings['fn_load_tws_dam_gwd'] = settings['dir_gwd_prep'] + 'load_tws_dam_gwd.nc'

    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    settings['fn_love'] = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    settings['tstart'] = dt.datetime.now().replace(microsecond=0)
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_grd_tws_natural_indiv, range(100))
    compute_grd_dam_gwd(settings)
    return

def compute_grd_tws_natural_indiv(ens):
    # For each Humphrey ensemble member, solve SLE and save
    global settings
    print('      Ensemble '+str(ens+1)+'/100 Elapsed time: ',(dt.datetime.now().replace(microsecond=0) - settings['tstart']))
    fn = settings['dir_gwd_prep'] + 'load_tws_natural_' + str(ens) + '.nc'
    file_handle = Dataset(fn)
    file_handle.set_auto_mask(False)
    lat = file_handle.variables['y'][:]
    lon = file_handle.variables['x'][:]
    time = file_handle.variables['t'][:]
    load = file_handle.variables['z'][:]
    file_handle.close()

    # Love numbers
    love = np.load(settings['fn_love'],allow_pickle=True).all()
    mask = np.load(settings['fn_mask'], allow_pickle=True).all()

    load = load - load.mean(axis=0)[np.newaxis,:,:]
    sle = pySLE.solver(lat=lat, lon=lon, time=time, load=load, slm=1.0-mask['land'], love=love, lmax=359)
    sle.solve()

    # Save data
    fn = settings['dir_gwd_prep'] + 'grd_tws_natural_' + str(ens) + '.nc'
    save_grd_tws_indiv(fn, sle)
    return

def compute_grd_dam_gwd(settings):
    # Read dam impoundment and GWD load and compute SLE
    file_handle = Dataset(settings['fn_load_tws_dam_gwd'])
    file_handle.set_auto_mask(False)
    lat = file_handle.variables['y'][:]
    lon = file_handle.variables['x'][:]
    time = file_handle.variables['t'][:]
    load_gwd_wada = file_handle.variables['GWD_wada'][:]
    load_gwd_doll = file_handle.variables['GWD_doll'][:]
    load_dam = file_handle.variables['Dams'][:]

    love = np.load(settings['fn_love'], allow_pickle=True).all()
    mask = np.load(settings['fn_mask'], allow_pickle=True).all()

    sle_gwd_wada = pySLE.solver(lat=lat, lon=lon, time=time, load=load_gwd_wada, slm=1.0 - mask['land'], love=love, lmax=359)
    sle_gwd_doll = pySLE.solver(lat=lat, lon=lon, time=time, load=load_gwd_doll, slm=1.0 - mask['land'], love=love, lmax=359)
    sle_dam = pySLE.solver(lat=lat, lon=lon, time=time, load=load_dam, slm=1.0 - mask['land'], love=love, lmax=359)

    sle_gwd_wada.solve()
    sle_gwd_doll.solve()
    sle_dam.solve()

    # Save data
    save_grd_tws_indiv(settings['dir_gwd_prep'] + 'grd_tws_dam.nc', sle_dam)
    save_grd_tws_indiv(settings['dir_gwd_prep'] + 'grd_tws_gwd_wada.nc', sle_gwd_wada)
    save_grd_tws_indiv(settings['dir_gwd_prep'] + 'grd_tws_gwd_doll.nc', sle_gwd_doll)
    return

# Helper functions
def save_grd_tws_indiv(fn, sle):
    file_handle = Dataset(fn, 'w')
    file_handle.createDimension('x', len(sle.result['lon']))
    file_handle.createDimension('y', len(sle.result['lat']))
    file_handle.createDimension('t', len(sle.result['time']))
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = sle.result['lon']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = sle.result['lat']
    file_handle.createVariable('t', 'i4', ('t',),zlib=True)[:] = sle.result['time']
    file_handle.createVariable('rsl', 'f4', ('t', 'y', 'x',),zlib=True,complevel=4,least_significant_digit=3)[:] = 1000*sle.result['rsl']
    file_handle.createVariable('rad', 'f4', ('t', 'y', 'x',),zlib=True,complevel=4,least_significant_digit=3)[:] =1000*sle.result['rad']
    file_handle.createVariable('barystatic', 'f4', ('t',),zlib=True,complevel=4,least_significant_digit=4)[:] = 1000*sle.result['barystatic']
    file_handle.close()
    return

if __name__ == '__main__':
    main()
