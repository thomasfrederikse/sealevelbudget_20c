# ----------------------------------------
# Read NASA MEASURES annual altimetry
# Correct for GIA and GRD effects
# and compute basin-mean and global-mean
# sea-level changes
# ----------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import mod_gentools as gentools
import multiprocessing as mp
import ctypes as ct

def main():
    set_settings()
    read_altimetry()
    read_mask()
    generate_gsl_samples()
    compute_basin_rsl()
    save_data()
    return

def set_settings():
    print('Define settings...')
    global settings
    settings = {}
    settings['test_run_ICE6G_D'] = True
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
    else:
        settings['nproc'] = 128
    settings['dir_budget']  = settings['dir_data'] + 'Budget_20c/'
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    settings['fn_altimetry'] = settings['dir_budget']+'vlm/Altimetry_annual.nc'
    if settings['test_run_ICE6G_D']:
        settings['dir_grd'] = settings['dir_budget'] + 'grd_ICE6G/'
        settings['fn_gia_rad'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['fn_gia_rsl'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['fn_alt_ens'] = settings['dir_budget']+'results/alt_basin_global_ens_ice6g.npy'
    else:
        settings['dir_grd'] = settings['dir_budget'] + 'grd/'
        settings['fn_gia_rad'] = settings['dir_data']+'GIA/Caron/Ensemble/rad_ens_05.nc'
        settings['fn_gia_rsl'] = settings['dir_data']+'GIA/Caron/Ensemble/rsl_ens_05.nc'
        settings['fn_alt_ens'] = settings['dir_budget']+'results/alt_basin_global_ens.npy'
    settings['years'] = np.arange(1993,2019)
    settings['years_total'] = np.arange(1900,2019)
    settings['num_ens']  = 100
    return

def save_data():
    global alt_basin, settings
    altimetry = {}
    altimetry['basin'] = np.zeros(6,dtype=object)
    for basin in range(6):
        altimetry['basin'][basin] = np.zeros([settings['num_ens'],len(settings['years_total'])],dtype=np.float32)*np.nan
        altimetry['basin'][basin][:,np.in1d(settings['years_total'],settings['years'])] = alt_basin['ensembles'][:,basin,:].squeeze()
    altimetry['global'] = np.zeros([settings['num_ens'],len(settings['years_total'])],dtype=np.float32)*np.nan
    altimetry['global'][:,np.in1d(settings['years_total'],settings['years'])] = alt_basin['ensembles'][:,6,:].squeeze()
    np.save(settings['fn_alt_ens'],altimetry)
    return

def compute_basin_rsl():
    global mask, altimetry, alt_basin, settings
    # Prepare storage
    alt_basin = {}
    alt_basin['ensembles'] = mp_empty_float([settings['num_ens'],7,len(settings['years'])])

    # Ensemble run
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_basin_rsl_indiv, range(settings['num_ens']))
    return

def compute_basin_rsl_indiv(ens):
    print('   '+str(ens))
    global mask, altimetry, alt_basin, settings, gsl_samples
    rad_gia_ens = read_gia_ens(ens)
    grd_ens = read_grd_ens(ens)
    alt_rsl = altimetry['ssh'] - (settings['years']-settings['years'].mean())[:,np.newaxis,np.newaxis] * rad_gia_ens[np.newaxis,...] - grd_ens['rad']
    # Average over mask
    for basin in range(6):
        mask_lcl = (mask['basin']==basin)
        alt_rsl_basin  = np.nansum((mask['area']*mask_lcl)[np.newaxis,:,:]*alt_rsl,axis=(1,2)) / (mask['area']*mask_lcl).sum()
        alt_rsl_basin += gsl_samples[ens,:]
        alt_rsl_basin -= alt_rsl_basin[-17:].mean()
        alt_basin['ensembles'][ens,basin,:] = alt_rsl_basin
    # Global
    mask_lcl = (np.isfinite(mask['basin']))
    alt_rsl_basin = np.nansum((mask['area'] * mask_lcl)[np.newaxis, :, :] * alt_rsl, axis=(1, 2)) / (mask['area'] * mask_lcl).sum()
    alt_rsl_basin += gsl_samples[ens, :]
    alt_rsl_basin -= alt_rsl_basin[-17:].mean()
    alt_basin['ensembles'][ens, 6, :] = alt_rsl_basin
    return

def read_altimetry():
    global settings, altimetry
    print('Reading altimetry...')
    altimetry = {}
    file_handle = Dataset(settings['fn_altimetry'], 'r')
    file_handle.set_auto_mask(False)
    altimetry['lat']  =  mp_filled_float(file_handle.variables['y'][:])
    altimetry['lon']  =  mp_filled_float(file_handle.variables['x'][:])
    altimetry['time'] =  mp_filled_float(file_handle.variables['t'][:])
    altimetry['ssh'] = mp_filled_float(file_handle.variables['z'][:])
    file_handle.close()
    altimetry['slm'] = mp_filled_bool(np.isfinite(altimetry['ssh'][-1,:,:]))
    return

def read_gia_ens(ens):
    global gia, settings
    file_handle = Dataset(settings['fn_gia_rad'], 'r')
    file_handle.set_auto_mask(False)
    if settings['test_run_ICE6G_D']:
        rad = file_handle.variables['rad'][:]
    else:
        rad = file_handle.variables['rad'][ens,:,:]
    file_handle.close()
    return(rad)

def read_mask():
    global settings, mask
    mask = {}
    mask_raw = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask['lon']   = mp_filled_float(mask_raw['lon'])
    mask['lat']   = mp_filled_float(mask_raw['lat'])
    mask['basin'] = mp_filled_float(mask_raw['basin'])
    mask['area']  = mp_filled_float(gentools.grid_area(mask['lat'],mask['lon']))
    return

def read_grd_ens(ens):
    global settings
    fn = settings['dir_grd']+'grd_'+str(ens)+'.nc'
    time_alt_idx = np.zeros(len(settings['years_total']),dtype=bool)
    time_alt_idx[np.in1d(settings['years_total'],settings['years'])] = True
    grd = {}
    grd['rad'] = Dataset(fn,'r').variables['rad'][time_alt_idx,:,:]._get_data()
    return(grd)

# Altimetry measurement uncertainty
def generate_gsl_samples():
    global settings, gsl_samples
    # Generate random fluctuations
    t_rand_highf  = gen_autocov(1.7,  1.5, 1.2, 2 / 12, settings)
    t_rand_medf   = gen_autocov(1.3,  1.2, 1.0, 1, settings)
    t_rand_wettr  = gen_autocov(1.1,  1.1, 1.1, 5, settings)
    t_rand_lorbit = gen_autocov(1.12, 0.5, 0.5, 10, settings)
    t_rand_intmis = gen_randjump(settings)
    t_rand_dorbit = gen_longterm_drift(0.10,settings)
    t_rand_dtopex = gen_topex_drift(settings)
    gsl_samples = t_rand_highf + t_rand_medf + t_rand_wettr + t_rand_lorbit + t_rand_intmis + t_rand_dorbit + t_rand_dtopex
    gsl_samples -=  gsl_samples[:,:12].mean(axis=1)[:,np.newaxis]
    # trend_ens = np.zeros(5000)
    # for n in range(5000): trend_ens[n] = gentools.lsqtrend(settings['years'],gsl_samples[n,:])
    return

def gen_randjump(settings):
    # Jump indices
    jump_0 = np.argmin(np.abs(settings['years'] - 2002))
    jump_1 = np.argmin(np.abs(settings['years'] - 2008))
    jump_2 = np.argmin(np.abs(settings['years'] - 2016))
    tpx_rnd = np.random.randn(settings['num_ens']) * 2
    j1_rnd = np.random.randn(settings['num_ens']) * 0.5
    j2_rnd = np.random.randn(settings['num_ens']) * 0.5
    t_rand = np.zeros([settings['num_ens'], len(settings['years'])])
    t_rand[:, jump_0:] = tpx_rnd[:, np.newaxis]
    t_rand[:, jump_1:] = t_rand[:, jump_1:] + j1_rnd[:, np.newaxis]
    t_rand[:, jump_2:] = t_rand[:, jump_2:] + j2_rnd[:, np.newaxis]
    return(t_rand)

def gen_longterm_drift(drift,settings):
    t_rand = drift * np.random.randn(settings['num_ens'])[:, np.newaxis] * (settings['years'] - np.mean(settings['years']))[np.newaxis, :]
    return(t_rand)

def gen_topex_drift(settings):
    dstop_A = np.argmin(np.abs(settings['years'] - 1998))
    dstop_B = np.argmin(np.abs(settings['years'] - 2002))
    drift_A = np.zeros([settings['num_ens'], len(settings['years'])])
    drift_B = np.zeros([settings['num_ens'], len(settings['years'])])
    drift_A[:,:dstop_A] = 0.7 * np.random.randn(settings['num_ens'])[:, np.newaxis] * (settings['years'][:dstop_A] - settings['years'][dstop_A])[np.newaxis, :]
    drift_B[:,:dstop_B] = 0.1 * np.random.randn(settings['num_ens'])[:, np.newaxis] * (settings['years'][:dstop_B] - settings['years'][dstop_B])[np.newaxis, :]
    t_rand = drift_A + drift_B
    return(t_rand)

def gen_autocov(sigma_TPX,sig_J1, sig_J23, l_factor, settings):
    jump_0 = np.argmin(np.abs(settings['years'] - 2002))
    jump_1 = np.argmin(np.abs(settings['years'] - 2008))
    sig_array = np.zeros(len(settings['years']))
    sig_array[:jump_0]       = sigma_TPX
    sig_array[jump_0:jump_1] = sig_J1
    sig_array[jump_1:]       = sig_J23
    t_distance = np.abs(settings['years'][:, np.newaxis] - settings['years'][np.newaxis, :])
    covmat = np.exp(-0.5 * (t_distance / (l_factor)) ** 2) + np.eye(len(settings['years'])) * 1.0e-10
    covc = np.linalg.cholesky(covmat)
    t_rand = np.zeros([settings['num_ens'], len(settings['years'])])
    for n in range(settings['num_ens']):
        t_rand[n, :] = sig_array * np.matmul(covc, np.random.randn(len(settings['years'])))
    return(t_rand)

# MP functions
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