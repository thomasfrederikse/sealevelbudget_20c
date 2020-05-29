# --------------------------------------------------------------
# Read all fingerprints, GIA, and compute spatial trend patterns
# 1900 - 2017 GIA + GRD
# 1957 - 2017 GRD + Steric
# 1993 - 2017 GRD + Steric
# --------------------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import mod_gentools as gentools
from scipy.interpolate import griddata

def main():
    global results, settings
    set_settings()
    results = {}
    # GRD
    compute_grd_stats()
    # GIA
    compute_gia_stats()
    # Steric
    compute_steric_stats()

    save_fields()
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt']    = os.getenv('HOME') + '/Scripts/GMT/Presentations/2019_09_CLIVAR/budget_map/'
    settings['fn_basin']   = settings['dir_data']+'Basins/ocean_regions_thompson.grd'
    settings['fn_gia_rsl'] = settings['dir_data']+'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_grd_rsl'] = settings['dir_budget']+'results/grd_stats.nc'

    settings['fn_CZ16']      = settings['dir_data'] + 'Steric/Cheng/cheng_steric_1940_2019.nc'
    settings['fn_I17']       = settings['dir_data'] + 'Steric/I17/I17_1955_2018.nc'
    settings['fn_WOA']   = settings['dir_data'] + 'Steric/Levitus/Levitus_1957_2018.nc'
    settings['fn_save']      = settings['dir_budget'] + 'results/spatial_trends.nc'

    settings['steric_grids'] = ['CZ16','WOA','I17']
    settings['years'] = np.arange(1900,2019)
    settings['years_grid'] = np.arange(1957,2019)
    settings['num_ens'] = 5000
    settings['lon'] = np.arange(0.25,360.25,0.5)
    settings['lat'] = np.arange(-89.75,90.25,0.5)
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability']/settings['probability'].sum()
    return


def compute_steric_stats():
    global results, settings
    trend_ens_1957_2018 = np.zeros([len(settings['steric_grids']),len(settings['lat']),len(settings['lon'])])
    trend_ens_1993_2018 = np.zeros([len(settings['steric_grids']),len(settings['lat']),len(settings['lon'])])

    for idx, pname in enumerate(settings['steric_grids']):
        read_steric_grid(pname)
        trend_ens_1957_2018[idx,:,:],trend_ens_1993_2018[idx,:,:] = read_steric_grid(pname)
    results['steric_1957_2018'] = trend_ens_1957_2018.mean(axis=0)
    results['steric_1993_2018'] = trend_ens_1993_2018.mean(axis=0)
    return

def read_steric_grid(pname):
    global settings
    file_handle = Dataset(settings['fn_' + pname])
    file_handle.set_auto_mask(False)
    time = file_handle.variables['t'][:]
    lon = file_handle.variables['x'][:]
    lat = file_handle.variables['y'][:]
    slm = file_handle.variables['slm'][:]

    time_acc = (time < settings['years_grid'][-1] + 1) & (time >= settings['years_grid'][0])
    time = time[time_acc]
    steric_monthly = file_handle.variables['h_totalsteric'][time_acc, :, :]
    file_handle.close()

    # To annual data
    steric_annual = np.zeros([len(settings['years_grid']), len(lat), len(lon)])
    for idx, yr in enumerate(settings['years_grid']):
        acc_idx = (time >= yr) & (time < yr + 1)
        steric_annual[idx, :, :] = steric_monthly[acc_idx, :, :].mean(axis=0)
    # Compute trends
    acc_1993_2016 = np.in1d(settings['years_grid'],np.arange(1993,2017))
    amat = np.ones([len(settings['years_grid']),2])
    amat[:,1] = settings['years_grid'][:] - settings['years_grid'].mean()
    amat_1957_2016_T  = amat.T
    amat_1957_2016_sq = np.linalg.inv(np.dot(amat_1957_2016_T, amat))
    amat_1993_2016_T  = amat[acc_1993_2016,:].T
    amat_1993_2016_sq = np.linalg.inv(np.dot(amat_1993_2016_T, amat[acc_1993_2016,:]))

    trend_1957_2016 = np.zeros(slm.shape)*np.nan
    trend_1993_2016 = np.zeros(slm.shape)*np.nan
    for i in range(len(lat)):
        for j in range(len(lon)):
            if np.isfinite(steric_annual[0,i,j]):
                trend_1957_2016[i,j] = np.dot(amat_1957_2016_sq, np.dot(amat_1957_2016_T, steric_annual[:,i,j]))[1]
                trend_1993_2016[i,j] = np.dot(amat_1993_2016_sq, np.dot(amat_1993_2016_T, steric_annual[acc_1993_2016,i,j]))[1]

    # Common grid using griddata
    lonmat,latmat = np.meshgrid(lon,lat)
    coord_in_array  = np.vstack([lonmat.ravel(), latmat.ravel()]).T
    lonmat,latmat = np.meshgrid(settings['lon'],settings['lat'])
    coord_out_array  = np.vstack([lonmat.ravel(), latmat.ravel()]).T
    trend_1957_2016_regrid  = griddata(coord_in_array,trend_1957_2016.flatten(),coord_out_array,method='nearest').reshape(360,720)
    trend_1993_2016_regrid  = griddata(coord_in_array,trend_1993_2016.flatten(),coord_out_array,method='nearest').reshape(360,720)
    return(trend_1957_2016_regrid,trend_1993_2016_regrid)

def compute_grd_stats():
    global results, settings
    file_handle = Dataset(settings['fn_grd_rsl'], 'r')
    file_handle.set_auto_mask(False)
    rsl = file_handle.variables['total_rsl_mean'][:]
    file_handle.close()
    gentools.field_trend(settings['years'],rsl)

    results['grd_1900_2018'] = gentools.field_trend(settings['years'],rsl)
    results['grd_1957_2018'] = gentools.field_trend(settings['years'][55:],rsl[55:,...])
    results['grd_1993_2018'] = gentools.field_trend(settings['years'][93:],rsl[93:,...])
    return

def compute_gia_stats():
    global results, settings
    file_handle = Dataset(settings['fn_gia_rsl'], 'r')
    file_handle.set_auto_mask(False)
    rsl = file_handle.variables['rsl'][:settings['num_ens'],:,:]
    file_handle.close()
    results['gia'] = (rsl * settings['probability'][:,np.newaxis,np.newaxis]).sum(axis=0)
    return

def save_fields():
    global results, settings
    file_handle = Dataset(settings['fn_save'], 'w')
    file_handle.createDimension('x', len(settings['lon']))
    file_handle.createDimension('y', len(settings['lat']))
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = settings['lon']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = settings['lat']
    for prod in results.keys():
        file_handle.createVariable(prod, 'f4', ('y', 'x',),zlib=True,complevel=4,least_significant_digit=4)[:] = results[prod]
    file_handle.close()
    return

def save_locations(settings):
    station_list = np.load(settings['fn_station_list'],allow_pickle=True)
    statloc = np.zeros([len(station_list),2])
    for idx, stat in enumerate(station_list):
        statloc[idx,0] = stat['lon']
        statloc[idx,1] = stat['lat']
    fn = settings['dir_gmt'] + 'statloc.txt'
    np.savetxt(fn,statloc,fmt='%4.3f;%4.3f')
    return


