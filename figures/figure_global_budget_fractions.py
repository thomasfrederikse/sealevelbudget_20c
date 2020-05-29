# ------------------------------------------------------------------
# For GMSL, estimate the fractions of each individual contributor to
# the total sum of contributors, both with and without GMSL
# ------------------------------------------------------------------
import os
import numpy as np
from netCDF4 import Dataset
import mod_gentools as gentools
def main():
    global settings
    set_settings()
    steric_ensembles = np.load(settings['fn_steric_ensembles'], allow_pickle=True).all()
    grd_ensembles = np.load(settings['fn_grd_ensembles'], allow_pickle=True).all()
    budget_total = steric_ensembles['global'][:settings['num_ens'],:] + grd_ensembles['global']['total'][:settings['num_ens'],:]
    budget_notws = steric_ensembles['global'][:settings['num_ens'],:] + grd_ensembles['global']['total'][:settings['num_ens'],:] - grd_ensembles['global']['tws'][:settings['num_ens'],:]

    mass_frac = ['glac', 'GrIS', 'AIS', 'tws','total']
    # Steric fractions
    frac_steric = estimate_fraction_from_ensemble(steric_ensembles['global'][:settings['num_ens'],:],budget_total)
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'frac_steric', settings['years'], frac_steric)
    frac_steric_notws = estimate_fraction_from_ensemble(steric_ensembles['global'][:settings['num_ens'],:],budget_notws)
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'frac_steric_notws', settings['years'], frac_steric_notws)
    # Mass fractions
    for process in mass_frac:
        frac_proc = estimate_fraction_from_ensemble(grd_ensembles['global'][process][:settings['num_ens'],:], budget_total)
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'frac_'+process, settings['years'], frac_proc)
        if process =='total': frac_proc_notws = estimate_fraction_from_ensemble(grd_ensembles['global'][process][:settings['num_ens'],:]-grd_ensembles['global']['tws'][:settings['num_ens'],:], budget_notws)
        else: frac_proc_notws = estimate_fraction_from_ensemble(grd_ensembles['global'][process][:settings['num_ens'],:], budget_notws)
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'frac_'+process+'_notws', settings['years'], frac_proc_notws)
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt'] = os.getenv('HOME') + '/Scripts/GMT/Papers/Budget_20c/Mass_steric_ratio/'
    settings['fn_steric_ensembles'] = settings['dir_budget'] + 'results/steric_basin_global_ens.npy'
    settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens.npy'
    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['years']     = np.arange(1900,2019)
    settings['num_ens']   = 5000
    settings['filt_width'] = 32
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability'] / settings['probability'].sum()
    return

def estimate_fraction_from_ensemble(ensemble_term,ensemble_total):
    # Estimate fraction of the sliding trend in ensemble_total explained by ensemble_term
    global settings
    fhalf = np.int(settings['filt_width'] / 2)
    ensemble_trend_term  = np.zeros(ensemble_term.shape) * np.nan
    ensemble_trend_total = np.zeros(ensemble_total.shape) * np.nan
    for idx, yr in enumerate(settings['years']):
        idx_start = np.max([0, idx - fhalf])
        idx_stop = np.min([len(settings['years']), idx + fhalf])
        yr_acc = settings['years'][idx_start:idx_stop]
        amat = np.ones([len(yr_acc), 2])
        amat[:, 1] = yr_acc
        amat_T = amat.T
        amat_sq = np.linalg.inv(np.dot(amat_T, amat))
        for ens in range(len(ensemble_term)):
            ensemble_trend_term[ens, idx]  = np.dot(amat_sq, np.dot(amat_T, ensemble_term[ens, idx_start:idx_stop]))[1]
            ensemble_trend_total[ens, idx] = np.dot(amat_sq, np.dot(amat_T, ensemble_total[ens, idx_start:idx_stop]))[1]
    ensemble_fraction = ensemble_trend_term/ensemble_trend_total
    sliding_fraction = np.zeros([len(settings['years']), 3]) * np.nan
    #sliding_fraction[:, 1] = np.sum(settings['probability'][:, np.newaxis] * ensemble_fraction, axis=0)
    for t in range(len(settings['years'])):
        sort_idx = np.argsort(ensemble_fraction[:, t])
        sort_cdf = np.cumsum(settings['probability'][sort_idx])
        sliding_fraction[t, 0] = ensemble_fraction[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.05))]
        sliding_fraction[t, 1] = ensemble_fraction[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.50))]
        sliding_fraction[t, 2] = ensemble_fraction[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.95))]
    sliding_fraction[:np.rint(settings['filt_width'] / 4).astype(int), :] = np.nan
    sliding_fraction[-np.rint(settings['filt_width'] / 4).astype(int):, :] = np.nan
    return (sliding_fraction)
