# Save the individual mass estimates for GMT
import numpy as np
import os

def main():
    settings = {}
    settings['dir_data']    = os.getenv('HOME')+'/Data/'
    settings['dir_GMT']    = os.getenv('HOME')+'/Scripts/GMT/Papers/Budget_20c/indiv_mass/'

    settings['fn_indiv'] = settings['dir_data'] + 'Budget_20c/results/indiv_mass.npy'
    settings['years']  = np.arange(1900,2019)
    indiv = np.load(settings['fn_indiv'],allow_pickle=True).all()

    for src in indiv:
        for model in indiv[src]:
            fn = settings['dir_GMT'] + src +'_'+model+'.txt'
            if indiv[src][model].ndim==1:
                np.savetxt(fn, np.vstack([settings['years'],indiv[src][model]]).T,fmt='%.3f,%.3f')
            else:
                np.savetxt(fn,np.vstack([settings['years'],indiv[src][model][:,[1,0,2]].T]).T,fmt='%.3f,%.3f,%.3f,%.3f')




