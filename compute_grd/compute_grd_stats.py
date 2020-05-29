# -------------------------------------------------
# Compute mean and standard error from GRD ensemble
# For barystatic, rad and RSL in all four terms
# -------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import multiprocessing as mp
import ctypes as ct

def main():
    set_settings()
    compute_mean_std()
    save_data()
    return

def set_settings():
    print('Settings...')
    global settings
    settings = {}
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
    else:
        settings['nproc'] = 10
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    settings['dir_budget']  = settings['dir_data'] + 'Budget_20c/'
    settings['fn_gia_rsl']  = settings['dir_data']+'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_save_grd'] = settings['dir_budget']+'results/grd_stats.nc'
    settings['years'] = np.arange(1900,2019)
    settings['num_ens']  = 5000
    settings['probability'] = Dataset(settings['fn_gia_rsl'],'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability']/settings['probability'].sum()
    settings['lat'] = mp_filled_float(Dataset(settings['fn_gia_rsl'],'r').variables['y'][:]._get_data())
    settings['lon'] = mp_filled_float(Dataset(settings['fn_gia_rsl'],'r').variables['x'][:]._get_data())
    settings['processes'] = ['glac','GrIS','AIS','tws','total']
    # Split task
    settings['tasks'] = []
    for process in settings['processes']:
        settings['tasks'].append([process,'rsl'])
        settings['tasks'].append([process,'rad'])
    return

def save_data():
    print('Saving data...')
    global settings, result_list
    file_handle = Dataset(settings['fn_save_grd'], 'w')
    file_handle.createDimension('lon', len(settings['lon']))
    file_handle.createDimension('lat', len(settings['lat']))
    file_handle.createDimension('time', len(settings['years']))
    file_handle.createVariable('lon', 'f4', ('lon',),zlib=True)[:] = settings['lon']
    file_handle.createVariable('lat', 'f4', ('lat',),zlib=True)[:] = settings['lat']
    file_handle.createVariable('time', 'i2', ('time',),zlib=True)[:] = settings['years']
    for process in settings['processes']:
        file_handle.createVariable(process+'_rsl_mean', 'i2', ('time', 'lat', 'lon',),zlib=True,complevel=4)[:] = 50*result_list[process]['rsl_mean']
        file_handle.createVariable(process+'_rsl_sterr','i2', ('time', 'lat', 'lon',),zlib=True,complevel=4)[:] = 50*result_list[process]['rsl_sterr']
        file_handle.createVariable(process+'_rad_mean', 'i2', ('time', 'lat', 'lon',),zlib=True,complevel=4)[:] = 50*result_list[process]['rad_mean']
        file_handle.createVariable(process+'_rad_sterr','i2', ('time', 'lat', 'lon',),zlib=True,complevel=4)[:] = 50*result_list[process]['rad_sterr']
        file_handle.variables[process+'_rsl_mean'].setncattr('scale_factor', 0.02)
        file_handle.variables[process+'_rsl_sterr'].setncattr('scale_factor', 0.02)
        file_handle.variables[process+'_rad_mean'].setncattr('scale_factor', 0.02)
        file_handle.variables[process+'_rad_sterr'].setncattr('scale_factor', 0.02)
    file_handle.close()
    return

def compute_mean_std():
    print('Computing mean and std...')
    # 1. Individual tasks
    global settings, result_list, field_mean, field_ste
    result_list = {}
    for process in settings['processes']:
        result_list[process] = {}
        result_list[process]['rsl_mean']  = mp_empty_float([119,360,720])
        result_list[process]['rsl_sterr'] = mp_empty_float([119,360,720])
        result_list[process]['rad_mean']  = mp_empty_float([119,360,720])
        result_list[process]['rad_sterr'] = mp_empty_float([119,360,720])
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_mean_indiv_task, range(len(settings['tasks'])))
    out  = pool.map(compute_ste_indiv_task, range(len(settings['tasks'])))
    for process in settings['processes']:
        result_list[process]['rsl_sterr'] = np.sqrt(result_list[process]['rsl_sterr'])
        result_list[process]['rad_sterr'] = np.sqrt(result_list[process]['rad_sterr'])
    return

def compute_mean_indiv_task(task):
    global settings, file_handle, result_list
    for ens in range(settings['num_ens']):
        print('   MEAN task '+str(task)+' Ensemble '+str(ens))
        if settings['tasks'][task][0]=='total': fn= settings['dir_budget']+'grd/grd_'+str(ens)+'.nc'
        else: fn= settings['dir_budget']+'grd/grd_'+settings['tasks'][task][0]+'_'+str(ens)+'.nc'
        fh = Dataset(fn, 'r')
        fh.set_auto_mask(False)
        result_list[settings['tasks'][task][0]][settings['tasks'][task][1]+'_mean'] += fh.variables[settings['tasks'][task][1]][:] * settings['probability'][ens]
        fh.close()
    return

def compute_ste_indiv_task(task):
    global settings, file_handle, result_list
    for ens in range(settings['num_ens']):
        print('   STERR task '+str(task)+' Ensemble '+str(ens))
        if settings['tasks'][task][0]=='total': fn= settings['dir_budget']+'grd/grd_'+str(ens)+'.nc'
        else: fn= settings['dir_budget']+'grd/grd_'+settings['tasks'][task][0]+'_'+str(ens)+'.nc'
        fh = Dataset(fn, 'r')
        fh.set_auto_mask(False)
        result_list[settings['tasks'][task][0]][settings['tasks'][task][1]+'_sterr'] += settings['probability'][ens] * (fh.variables[settings['tasks'][task][1]][:]-result_list[settings['tasks'][task][0]][settings['tasks'][task][1]+'_mean'])**2
        fh.close()
    return

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

