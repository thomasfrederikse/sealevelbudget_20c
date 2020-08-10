# --------------------------------------------------------------------
# Determine whether the trends in 21st-century sea level, steric, mass
# are significantly higher than any moment before
# --------------------------------------------------------------------

import os
import numpy as np
from netCDF4 import Dataset
import mod_hector as hector

def main():
    global settings
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'

    settings['fn_steric_ensembles'] = settings['dir_budget'] + 'results/steric_basin_global_ens.npy'
    settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens.npy'
    settings['fn_obs_ensembles'] = settings['dir_budget'] + 'results/obs_basin_global_ens.npy'

    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_save_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'

    settings['years'] = np.arange(1900, 2019)
    settings['num_ens'] = 5000
    settings['filt_width'] = 30  # Moving window when estimating running trends
    settings['est_serial_corr'] = True  # Estimate serial correlation in linear trends
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability'] / settings['probability'].sum()

    # Load ensemble members
    obs_ensembles = np.load(settings['fn_obs_ensembles'], allow_pickle=True).all()
    steric_ensembles = np.load(settings['fn_steric_ensembles'], allow_pickle=True).all()
    grd_ensembles = np.load(settings['fn_grd_ensembles'], allow_pickle=True).all()

    # Estimate sliding trend of each ensemble member
    obs_ens = estimate_sliding_trend_ensemble(obs_ensembles['global']['tg_full'])
    grd_ens = estimate_sliding_trend_ensemble(grd_ensembles['global']['total'])
    steric_ens = estimate_sliding_trend_ensemble(steric_ensembles['global'])

    # Determine probability that the trend post-2000 is higher trhan the trand at any moment before
    signif_steric = determine_significance(steric_ens)
    signif_grd = determine_significance(grd_ens)
    signif_obs = determine_significance(obs_ens)
    return

def estimate_sliding_trend_ensemble(ensemble):
    # Estimate sliding trend from ensemble
    global settings
    time_acc_idx = np.isfinite(ensemble[0, :])
    ensemble = ensemble[:, time_acc_idx]
    time_lcl = settings['years'][time_acc_idx]
    fhalf = np.int(settings['filt_width'] / 2)
    trend_ens = np.zeros(ensemble.shape) * np.nan
    for idx, yr in enumerate(time_lcl):
        idx_start = np.max([0, idx - fhalf])
        idx_stop = np.min([len(time_lcl), idx + fhalf])
        yr_acc = time_lcl[idx_start:idx_stop]
        amat = np.ones([len(yr_acc), 2])
        amat[:, 1] = yr_acc
        amat_T = amat.T
        amat_sq = np.linalg.inv(np.dot(amat_T, amat))
        for ens in range(len(ensemble)):
            trend_ens[ens, idx] = np.dot(amat_sq, np.dot(amat_T, ensemble[ens, idx_start:idx_stop]))[1]
    trend_ens[:,:np.rint(settings['filt_width'] / 4).astype(int)] = np.nan
    trend_ens[:,-np.rint(settings['filt_width'] / 4).astype(int):] = np.nan
    return (trend_ens)

def determine_significance(trend_ens):
    trend_diff = np.nanmax(trend_ens[:,99:],axis=1) - np.nanmax(trend_ens[:,:85],axis=1)

    trend_diff_bnd = np.zeros(3)
    trend_diff_bnd[0] = (settings['probability'] * trend_diff).sum()

    sort_idx = np.argsort(trend_diff)
    sort_cdf = np.cumsum(settings['probability'][sort_idx])
    trend_diff_bnd[1] = trend_diff[sort_idx][np.argmin(np.abs(sort_cdf - 0.05))]
    trend_diff_bnd[2] = trend_diff[sort_idx][np.argmin(np.abs(sort_cdf - 0.95))]
    return(trend_diff_bnd)

