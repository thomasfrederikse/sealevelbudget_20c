# -------------------------------
# Save steric time series for GMT
# Mean and standard deviation of
# final estimate
# as well as individual products
# -------------------------------
import os
import numpy as np
import mod_gentools as gentools
def main():
    set_settings()
    write_data()
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt'] = os.getenv('HOME') + '/Scripts/GMT/Papers/Budget_20c/Steric_tseries/'
    settings['fn_steric_global_indiv_products']   = settings['dir_data'] + 'Budget_20c/results/steric_global_indiv_products.npy'
    settings['fn_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
    settings['years'] = np.arange(1900,2019)
    return

def write_data():
    global settings
    stats = np.load(settings['fn_stats'],allow_pickle=True).all()
    steric_indiv = np.load(settings['fn_steric_global_indiv_products'],allow_pickle=True).all()
    for prod in steric_indiv:
        if prod == 'Zanna':
            steric_tseries = np.zeros([len(settings['years']),3])
            steric_tseries[:,1] = steric_indiv[prod]['global'] - steric_indiv[prod]['global'][-17:].mean() + stats['global']['steric']['tseries'][-17:,1].mean()
            steric_tseries[:,0] = steric_tseries[:,1] - 1.65*steric_indiv[prod]['global_sterr']
            steric_tseries[:,2] = steric_tseries[:,1] + 1.65*steric_indiv[prod]['global_sterr']
            gentools.gmt_save_tseries_ci(settings['dir_gmt'], prod, settings['years'], steric_tseries)
        else:
            steric_tseries = steric_indiv[prod]['global'] - steric_indiv[prod]['global'][-17:].mean()
            save_array = np.vstack([steric_indiv[prod]['years'], steric_tseries]).T
            np.savetxt(settings['dir_gmt'] + prod + '_tseries_m.txt', save_array, fmt='%4.3f;%4.3f')
    gentools.gmt_save_tseries_ci(settings['dir_gmt'],'steric',settings['years'],stats['global']['steric']['tseries'])

    # Zanna deep
    steric_tseries = np.zeros([len(settings['years']), 3])
    steric_tseries[:, 1] = steric_indiv['Zanna']['deep'] - steric_indiv['Zanna']['deep'][-17:].mean()
    steric_tseries[:, 0] = steric_tseries[:, 1] - 1.65 * steric_indiv['Zanna']['deep_sterr']
    steric_tseries[:, 2] = steric_tseries[:, 1] + 1.65 * steric_indiv['Zanna']['deep_sterr']
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'Zanna_deep', settings['years'], steric_tseries)
    return