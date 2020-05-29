# -------------------------------------------
# Read basin-scale and global steric products
# I17 1955-2018
# C16 1955-2018
# L11 1957-2018
# Z18 1900-2018
# -------------------------------------------
import os
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp2d
import mod_gentools as gentools

def main():
    set_settings()
    read_mask_mp()
    compute_steric_ensembles()
    save_global_indiv_prod()
    save_data()
    return

def set_settings():
    print('Defining settings...')
    global settings
    settings = {}
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    settings['fn_CZ16']   = settings['dir_data'] + 'Steric/Cheng/cheng_steric_1940_2019.nc'
    settings['fn_I17']   = settings['dir_data'] + 'Steric/I17/I17_1955_2018.nc'
    settings['fn_WOA']   = settings['dir_data'] + 'Steric/Levitus/Levitus_1957_2018.nc'
    settings['fn_Zanna']   = settings['dir_data'] + 'Steric/Zanna/ThSL_GF_1870_2018.nc'
    settings['fn_steric_ens']   = settings['dir_data'] + 'Budget_20c/results/steric_basin_global_ens.npy'
    settings['fn_steric_global_indiv_products']   = settings['dir_data'] + 'Budget_20c/results/steric_global_indiv_products.npy'

    settings['grids'] = ['CZ16','WOA','I17']
    settings['years'] = np.arange(1900,2019)
    settings['years_grid'] = np.arange(1957,2019)
    settings['num_ens'] = 5000
    return

def save_data():
    print('Saving...')
    global steric_ens, settings
    steric = {}
    steric['basin'] = np.zeros(6,dtype=object)
    for basin in range(6):
        steric['basin'][basin] = np.zeros([settings['num_ens'],len(settings['years'])],dtype=np.float32)
        steric['basin'][basin] = steric_ens['basin'][basin]
    steric['global'] = np.zeros([settings['num_ens'],len(settings['years'])],dtype=np.float32)
    steric['global'][:] = steric_ens['global'][:]
    np.save(settings['fn_steric_ens'],steric)
    return


def compute_steric_ensembles():
    print('Processing grids...')
    global steric_ens, steric_grid, mask, settings
    steric_ens = {}

    # Read files
    steric_grid = {}
    for grid in settings['grids']:
        steric_grid[grid] = proc_grid_basin(grid)
    steric_grid['Zanna'] = read_zanna()

    # Random numbers for sampling
    rnd_ptb_zanna = (np.random.normal(0, 1, settings['num_ens'])[:,np.newaxis] + np.random.normal(0, 1, [settings['num_ens'],len(steric_grid['Zanna']['global_sterr'])]))/np.sqrt(2)

    random_prod   = np.random.randint(0, len(settings['grids']), settings['num_ens'])
    steric_ens['basin'] = np.zeros(6,dtype=object)
    grid_time_acc = np.in1d(settings['years'],settings['years_grid'])
    for basin in range(6):
        steric_ens['basin'][basin] = np.zeros([settings['num_ens'],len(settings['years'])])*np.nan
        for ens in range(settings['num_ens']):
            steric_ens['basin'][basin][ens,grid_time_acc] = steric_grid[settings['grids'][random_prod[ens]]]['basin'][basin]
            # Add deep steric contribution
            steric_ens['basin'][basin][ens,:] +=  steric_grid['Zanna']['deep'] + rnd_ptb_zanna[ens,:] * steric_grid['Zanna']['deep_sterr']
    # Global
    steric_ens['global'] = np.zeros([settings['num_ens'], len(settings['years'])]) * np.nan
    for ens in range(settings['num_ens']):
        # 1. 1955-2016 observed steric + deep steric Zanna
        steric_ens['global'][ens, grid_time_acc] = steric_grid[settings['grids'][random_prod[ens]]]['global']
        steric_ens['global'][ens, :] += steric_grid['Zanna']['deep'] + rnd_ptb_zanna[ens,:] * steric_grid['Zanna']['deep_sterr']

        # 2. 1900-1955: Zanna full
        #zanna_idx = np.in1d(settings['years'],np.arange(settings['years'][0],settings['years_grid'][10]))
        zanna_idx = np.in1d(settings['years'],1957)
        ovl_idx   = grid_time_acc & zanna_idx
        zanna_glb = steric_grid['Zanna']['global']
        zanna_glb = zanna_glb - zanna_glb[ovl_idx].mean() + steric_ens['global'][ens, ovl_idx].mean()+ rnd_ptb_zanna[ens,:] * steric_grid['Zanna']['global_sterr']
        steric_ens['global'][ens, ~grid_time_acc] = zanna_glb[~grid_time_acc]

    # Remove baseline 2000-2016
    for basin in range(6):
        steric_ens['basin'][basin] -= steric_ens['basin'][basin][:,-17:].mean(axis=1)[:,np.newaxis]
    steric_ens['global'] -= steric_ens['global'][:, -17:].mean(axis=1)[:, np.newaxis]
    return


def proc_grid_basin(pname):
    print('   Processing grid '+pname+'...')
    global settings, mask
    steric = {}
    file_handle = Dataset(settings['fn_'+pname])
    file_handle.set_auto_mask(False)
    time =  file_handle.variables['t'][:]
    lon  =  file_handle.variables['x'][:]
    lat  =  file_handle.variables['y'][:]
    slm  =  file_handle.variables['slm'][:]
    time_acc = (time < settings['years_grid'][-1]+1)&(time>=settings['years_grid'][0])
    time = time[time_acc]
    steric_monthly = file_handle.variables['h_totalsteric'][time_acc,:,:]
    halo_glb_monthly = file_handle.variables['ts_halosteric'][time_acc]
    file_handle.close()
    steric_monthly -= halo_glb_monthly[:,np.newaxis,np.newaxis]

    # To annual data
    steric['years'] = settings['years_grid']
    steric_annual = np.zeros([len(steric['years']),len(lat),len(lon)])
    for idx,yr in enumerate(steric['years']):
        acc_idx = (time>=yr) & (time<yr+1)
        steric_annual[idx,:,:] = steric_monthly[acc_idx,:,:].mean(axis=0)

    area = gentools.grid_area(lat,lon)
    steric['basin'] = np.zeros(6,dtype=object)
    # Basin mean
    for basin in range(6):
        msk_lcl = (mask['basin']==basin)
        msk_interp = np.rint(interp2d(mask['lon'],mask['lat'],msk_lcl,kind='linear',bounds_error=False,fill_value=0)(lon,lat))*slm
        steric['basin'][basin] = np.nansum((area[np.newaxis,:,:]*msk_interp[np.newaxis,:,:]*steric_annual),axis=(1,2)) / (area*msk_interp).sum()

    # Global steric
    msk_lcl = 1.0 * (np.isfinite(mask['basin']))
    msk_interp = np.rint(interp2d(mask['lon'], mask['lat'], msk_lcl, kind='linear', bounds_error=False, fill_value=0)(lon, lat)) * slm
    steric['global'] = np.nansum((area[np.newaxis, :, :] * msk_interp[np.newaxis, :, :] * steric_annual), axis=(1, 2)) / (area * msk_interp).sum()
    return(steric)


def save_global_indiv_prod():
    global settings, steric_grid
    np.save(settings['fn_steric_global_indiv_products'],steric_grid)
    return

def read_zanna():
    print('Processing grid Zanna...')
    global settings
    grid = {}
    file_handle = Dataset(settings['fn_Zanna'], 'r')
    file_handle.set_auto_mask(False)
    grid['years'] = np.arange(1900,2019)
    grid['global'] = file_handle.variables['ThSL_full_depth'][30:]
    grid['deep'] = file_handle.variables['ThSL_below_2000m'][30:]
    grid['global_sterr'] = file_handle.variables['error_ThSL_full_depth'][30:]
    grid['deep_sterr'] = file_handle.variables['error_ThSL_below_2000'][30:]
    grid['global'] = grid['global'] - grid['global'][-10:].mean()
    grid['deep'] = grid['deep'] - grid['deep'][-10:].mean()
    file_handle.close()
    return(grid)

def read_mask_mp():
    print('   Reading mask...')
    global mask, settings
    mask_raw = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask = {}
    mask['lat'] = mask_raw['lat']
    mask['lon'] = mask_raw['lon']
    mask['basin'] = mask_raw['basin']
    return

if __name__ == '__main__':
    main()
