# --------------------------------------------
# Compute GIA and GRD ensembles on basin scale
# and global scale
# GRD for each individual process
# for global: also dams, gwd and natural
# --------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import multiprocessing as mp
import ctypes as ct
import mod_gentools as gentools

def main():
    set_settings()
    read_mask_mp()
    sample_gia()
    sample_grd_basin()
    sample_grd_global()
    save_data()
    return

def set_settings():
    print('Defining settings...')
    global settings
    settings = {}
    settings['test_run_ICE6G_D'] = True
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
    else:
        settings['nproc'] = 128
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    if settings['test_run_ICE6G_D']:
        settings['dir_grd'] = settings['dir_budget'] + 'grd_ICE6G/'
        settings['fn_gia_rsl'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['fn_gia_ens'] = settings['dir_data'] + 'Budget_20c/results/gia_basin_global_ens_ice6g.npy'
        settings['fn_grd_ens'] = settings['dir_data'] + 'Budget_20c/results/grd_basin_global_ens_ice6g.npy'
    else:
        settings['dir_grd'] = settings['dir_budget'] + 'grd/'
        settings['fn_gia_rsl'] = settings['dir_data']+'GIA/Caron/Ensemble/rsl_ens_05.nc'
        settings['fn_gia_ens'] = settings['dir_data'] + 'Budget_20c/results/gia_basin_global_ens.npy'
        settings['fn_grd_ens'] = settings['dir_data'] + 'Budget_20c/results/grd_basin_global_ens.npy'
    settings['years'] = np.arange(1900, 2019)
    settings['num_ens'] = 100
    settings['terms'] = ['glac','GrIS','AIS','tws','total']
    settings['terms_tws'] = ['natural','gwd','dam']
    return

def save_data():
    global gia_basin, gia_global, grd_basin, grd_global, settings
    gia = {}
    gia['basin'] = np.zeros(6,dtype=object)
    for basin in range(6):
        gia['basin'][basin] = np.zeros(settings['num_ens'],dtype=np.float32)
        gia['basin'][basin][:] = gia_basin[:,basin]
    gia['global'] = np.zeros(settings['num_ens'], dtype=np.float32)
    gia['global'][:] = gia_global[:]
    grd = {}
    grd['basin'] = np.zeros(6,dtype=object)
    for basin in range(6):
        grd['basin'][basin] = {}
        for term in settings['terms']:
            grd['basin'][basin][term] = np.zeros([settings['num_ens'],len(settings['years'])],dtype=np.float32)
            grd['basin'][basin][term][:] = grd_basin[term][:,basin,:].squeeze()
    grd['global'] = {}
    for term in grd_global.keys():
        grd['global'][term] = np.zeros([settings['num_ens'],len(settings['years'])],dtype=np.float32)
        grd['global'][term][:] = grd_global[term][:]
    np.save(settings['fn_gia_ens'],gia)
    np.save(settings['fn_grd_ens'],grd)
    return

def read_mask_mp():
    print('Reading mask...')
    global mask, settings
    mask_raw = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask = {}
    mask['lat'] = mp_filled_float(mask_raw['lat'])
    mask['lon'] = mp_filled_float(mask_raw['lon'])
    mask['basin'] = mp_filled_float(mask_raw['basin'])
    return

def sample_gia():
    print('Sampling GIA...')
    global gia_basin, gia_global, settings, mask
    gia_basin = mp_empty_float([settings['num_ens'],6])
    gia_global = mp_empty_float([settings['num_ens']])
    # read gia
    file_handle = Dataset(settings['fn_gia_rsl'], 'r')
    file_handle.set_auto_mask(False)
    if settings['test_run_ICE6G_D']:
        rsl = file_handle.variables['RSL'][:][np.newaxis,:,:]
    else:
        rsl = file_handle.variables['rsl'][:settings['num_ens'],:,:]
    file_handle.close()
    # Compute basin-mean GIA
    area = gentools.grid_area(mask['lat'],mask['lon'])
    for basin in range(6):
        mask_basin = (mask['basin']==basin)
        gia_basin[:,basin] = area_average(rsl, mask_basin, area)
    mask_glb = (np.isfinite(mask['basin']))
    gia_global[:] = area_average(rsl, mask_glb, area)
    return

def sample_grd_basin():
    global grd_basin, settings, mask
    grd_basin = {}
    for term in settings['terms']: grd_basin[term] = mp_empty_float([settings['num_ens'],6,len(settings['years'])])
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(sample_grd_basin_indiv, range(settings['num_ens']))
    return

def sample_grd_basin_indiv(ens):
    global settings, mask, grd_basin
    print(ens)
    area = gentools.grid_area(mask['lat'],mask['lon'])
    for term in settings['terms']:
        if term == 'total': fn = settings['dir_grd']+'grd_'+str(ens)+'.nc'
        else:               fn = settings['dir_grd']+'grd_'+term+'_'+str(ens)+'.nc'
        file_handle = Dataset(fn, 'r')
        file_handle.set_auto_mask(False)
        rsl = file_handle.variables['rsl'][:]
        file_handle.close()
        for basin in range(6):
            mask_lcl = (mask['basin'] == basin)
            grd_basin[term][ens, basin, :] = area_average(rsl, mask_lcl, area)
    return

def sample_grd_global():
    global grd_global, settings
    grd_global = {}
    for term in settings['terms']: grd_global[term] = mp_empty_float([settings['num_ens'],len(settings['years'])])
    for term_tws in settings['terms_tws']: grd_global['tws_'+term_tws] = mp_empty_float([settings['num_ens'],len(settings['years'])])
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(sample_grd_global_indiv, range(settings['num_ens']))
    return

def sample_grd_global_indiv(ens):
    global settings, mask, grd_global
    print(ens)
    for term in settings['terms']:
        if term == 'total': fn = settings['dir_grd']+'grd_'+str(ens)+'.nc'
        else:               fn = settings['dir_grd']+'grd_'+term+'_'+str(ens)+'.nc'
        file_handle = Dataset(fn, 'r')
        file_handle.set_auto_mask(False)
        grd_global[term][ens, :] = file_handle.variables['qglb'][:]
        if term=='tws':
            for term_tws in settings['terms_tws']:
                grd_global['tws_'+term_tws][ens, :] = file_handle.variables['qglb_'+term_tws][:]
                grd_global['tws_'+term_tws][ens, 104:] = np.nan
        file_handle.close()
    return

def area_average(field,mask_lcl,area):
    field_average = (field * (mask_lcl * area)[np.newaxis, :, :]).sum(axis=(1, 2)) / (mask_lcl*area).sum()
    return field_average

# Parallel processing routines
def mp_empty_float(shape):
    shared_array_base = mp.RawArray(ct.c_float, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_empty_int(shape):
    shared_array_base = mp.RawArray(ct.c_int, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_float(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_float, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_bool(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_bool, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

if __name__ == '__main__':
    main()
