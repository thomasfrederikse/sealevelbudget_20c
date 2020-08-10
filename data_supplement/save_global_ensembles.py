# ---------------------------------------------
# Save save GMSL and component ensemble members
# as netcdf files
# ---------------------------------------------
import numpy as np
import os
from netCDF4 import Dataset

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_steric_ensembles'] = settings['dir_budget'] + 'results/steric_basin_global_ens.npy'
    settings['fn_gia_ensembles'] = settings['dir_budget'] + 'results/gia_basin_global_ens.npy'
    settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens.npy'
    settings['fn_obs_ensembles'] = settings['dir_budget'] + 'results/obs_basin_global_ens.npy'
    settings['fn_alt_ensembles'] = settings['dir_budget'] + 'results/alt_basin_global_ens.npy'
    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_save_ens'] = settings['dir_budget'] + 'data_supplement/GMSL_ensembles.nc'
    settings['years'] = np.arange(1900, 2019)
    settings['num_ens'] = 5000
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability'] / settings['probability'].sum()

    obs_ensembles = np.load(settings['fn_obs_ensembles'], allow_pickle=True).all()
    steric_ensembles = np.load(settings['fn_steric_ensembles'], allow_pickle=True).all()
    gia_ensembles = np.load(settings['fn_gia_ensembles'], allow_pickle=True).all()
    grd_ensembles = np.load(settings['fn_grd_ensembles'], allow_pickle=True).all()
    alt_ensembles = np.load(settings['fn_alt_ensembles'], allow_pickle=True).all()

    fh = Dataset(settings['fn_save_ens'], 'w')
    fh.createDimension('likelihood', settings['num_ens'])
    fh.createDimension('time', len(settings['years']))
    fh.copyright   = 'This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/.'
    fh.information = "This work is a data supplement to 'The causes of sea-level rise since 1900' by Thomas Frederikse, Felix Landerer, Lambert Caron, Surendra Adhikari, David Parkes, Vincent W. Humphrey, Soenke Dangendorf, Peter Hogarth, Laure Zanna, Lijing Cheng, Yun-Hao Wu'"
    fh_prob = fh.createVariable('likelihood', 'f4', ('likelihood'), zlib=True,complevel=4)
    fh_prob[:] = settings['probability']
    fh_prob.long_name = 'Normalized likelihood'
    fh_prob.units = '(-)'

    fh_time = fh.createVariable('time', 'i2', ('time'), zlib=True,complevel=4)
    fh_time[:] = settings['years']
    fh_time.long_name = 'Time (value represents year average)'
    fh_time.units = 'Year'

    fh_gmsl = fh.createVariable('GMSL', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*obs_ensembles['global']['tg_full']).astype(int)
    fh_gmsl[:] = save_ens
    fh_gmsl.missing_value = -20000
    fh_gmsl.scale_factor=0.05
    fh_gmsl.units = 'mm'
    fh_gmsl.long_name = 'Observed global-mean sea level'

    fh_bary = fh.createVariable('Barystatic', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*grd_ensembles['global']['total']).astype(int)
    fh_bary[:] = save_ens
    fh_bary.missing_value = -20000
    fh_bary.scale_factor=0.05
    fh_bary.units = 'mm'
    fh_bary.long_name = 'Barystatic sea level (total)'

    fh_bary = fh.createVariable('Glaciers', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*grd_ensembles['global']['glac']).astype(int)
    fh_bary[:] = save_ens
    fh_bary.missing_value = -20000
    fh_bary.scale_factor=0.05
    fh_bary.units = 'mm'
    fh_bary.long_name = 'Glaciers'

    fh_bary = fh.createVariable('GrIS', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*grd_ensembles['global']['GrIS']).astype(int)
    fh_bary[:] = save_ens
    fh_bary.missing_value = -20000
    fh_bary.scale_factor=0.05
    fh_bary.units = 'mm'
    fh_bary.long_name = 'Greenland Ice Sheet'

    fh_bary = fh.createVariable('AIS', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*grd_ensembles['global']['AIS']).astype(int)
    fh_bary[:] = save_ens
    fh_bary.missing_value = -20000
    fh_bary.scale_factor=0.05
    fh_bary.units = 'mm'
    fh_bary.long_name = 'Antarctic Ice Sheet'

    fh_bary = fh.createVariable('TWS', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*grd_ensembles['global']['tws']).astype(int)
    fh_bary[:] = save_ens
    fh_bary.missing_value = -20000
    fh_bary.scale_factor=0.05
    fh_bary.units = 'mm'
    fh_bary.long_name = 'Terrestrial water storage'

    fh_bary = fh.createVariable('Steric', 'i2', ('likelihood', 'time'), zlib=True, complevel=4)
    save_ens = np.rint(20*steric_ensembles['global']).astype(int)
    fh_bary[:] = save_ens
    fh_bary.missing_value = -20000
    fh_bary.scale_factor=0.05
    fh_bary.units = 'mm'
    fh_bary.long_name = 'Thermosteric sea-level anomalies'
    fh.close()