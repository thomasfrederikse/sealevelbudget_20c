# ------------------------------------------
# Read spatial estimates from ocean mass
# and steric fields, generate understandable
# metadata and save
# ------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
from scipy.interpolate import griddata

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_save'] = settings['dir_budget'] + 'data_supplement/'
    settings['fn_mask'] = settings['dir_budget'] +'grd_prep/mask.npy'

    settings['years'] = np.arange(1900, 2019)
    settings['years_steric'] = np.arange(1957,2019)
    settings['steric_products'] = ['CZ16','I17','WOA']
    settings['mass_terms'] = ['glac','GrIS','AIS','tws','total']
    settings['fn_CZ16']      = settings['dir_data'] + 'Steric/Cheng/cheng_steric_1940_2019.nc'
    settings['fn_I17']       = settings['dir_data'] + 'Steric/I17/I17_1955_2018.nc'
    settings['fn_WOA']       = settings['dir_data'] + 'Steric/Levitus/Levitus_1957_2018.nc'
    settings['lon'] = np.arange(0.25,360.25,0.5)
    settings['lat'] = np.arange(-89.75,90.25,0.5)

    process_steric(settings)
    for term in settings['mass_terms']: process_grd(term, settings)

    return

def process_grd(term,settings):
    print('  Processing '+term+'...')
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    fh = Dataset(settings['dir_budget']+'/results/grd_stats.nc','r')
    fh.set_auto_mask(False)
    grd = {}
    grd['rsl_mean']  = fh.variables[term+'_rsl_mean'][:]
    grd['rad_mean']  = fh.variables[term+'_rad_mean'][:]
    grd['rsl_sterr'] = fh.variables[term+'_rsl_sterr'][:]
    grd['rad_sterr'] = fh.variables[term+'_rad_sterr'][:]
    fh.close()

    grd['rsl_mean'][:,mask['land']]  = np.nan
    grd['rsl_sterr'][:,mask['land']] = np.nan
    save_nc(settings['years'], grd, term, settings)
    return

def process_steric(settings):
    steric_indiv = {}
    for prod in settings['steric_products']: steric_indiv[prod] = process_steric_indiv(prod,settings)

    steric ={}
    steric['mean'] = np.zeros([len(settings['years_steric']), len(settings['lat']), len(settings['lon'])])
    for prod in settings['steric_products']: steric['mean']+=(steric_indiv[prod]/len(settings['steric_products']))
    steric['mean']-=steric['mean'][(settings['years_steric']>1999),:,:].mean(axis=0)

    # Standard error
    steric['sterr'] = np.zeros([len(settings['years_steric']), len(settings['lat']), len(settings['lon'])])
    for prod in settings['steric_products']:  steric['sterr'] += (steric_indiv[prod]-steric['mean'])**2
    steric['sterr'] = np.sqrt(steric['sterr']/3)
    save_nc(settings['years_steric'], steric, 'steric', settings)
    return

def process_steric_indiv(prod,settings):
    # Read individual steric field, compute annual means
    # and regrid to the standard coordinates
    print('  Processing '+prod+'...')
    file_handle = Dataset(settings['fn_' + prod])
    file_handle.set_auto_mask(False)
    time = file_handle.variables['t'][:]
    lon = file_handle.variables['x'][:]
    lat = file_handle.variables['y'][:]
    time_acc = (time < settings['years_steric'][-1] + 1) & (time >= settings['years_steric'][0])
    time = time[time_acc]
    steric_monthly = file_handle.variables['h_totalsteric'][time_acc, :, :]
    file_handle.close()
    # To annual data
    steric_annual = np.zeros([len(settings['years_steric']), len(lat), len(lon)])
    for idx, yr in enumerate(settings['years_steric']):
        acc_idx = (time >= yr) & (time < yr + 1)
        steric_annual[idx, :, :] = steric_monthly[acc_idx, :, :].mean(axis=0)
    # Common grid using griddata
    lonmat,latmat = np.meshgrid(lon,lat)
    coord_in_array  = np.vstack([lonmat.ravel(), latmat.ravel()]).T
    lonmat,latmat = np.meshgrid(settings['lon'],settings['lat'])
    coord_out_array  = np.vstack([lonmat.ravel(), latmat.ravel()]).T
    steric_regrid = np.zeros([len(settings['years_steric']),len(settings['lat']),len(settings['lon'])])*np.nan
    print('    Regridding...')
    for idx, yr in enumerate(settings['years_steric']):
        steric_regrid[idx,...]  = griddata(coord_in_array,steric_annual[idx, :, :].flatten(),coord_out_array,method='nearest').reshape(360,720)
    return(steric_regrid)

def save_nc(time,grid,pname,settings):
    name_dict = {}
    name_dict['steric'] = 'Steric sea-level anomaly'
    name_dict['GrIS']   = 'Greenland Ice Sheet'
    name_dict['glac']   = 'Glaciers'
    name_dict['AIS']    = 'Antarctic Ice Sheet'
    name_dict['tws']    = 'Terrestrial water storage'
    name_dict['total']  = 'Glaciers, ice sheets, and terrestrial water storage'

    filename = settings['dir_save'] + pname+'.nc'
    fh = Dataset(filename, 'w')
    fh.copyright   = 'This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/.'
    fh.information = "This work is a data supplement to 'The causes of sea-level rise since 1900' by Thomas Frederikse, Felix Landerer, Lambert Caron, Surendra Adhikari, David Parkes, Vincent W. Humphrey, Soenke Dangendorf, Peter Hogarth, Laure Zanna, Lijing Cheng, Yun-Hao Wu'"
    fh.createDimension('lon', len(settings['lon']))
    fh.createDimension('lat', len(settings['lat']))
    fh.createDimension('time', len(time))
    fh_lon = fh.createVariable('lon', 'f4', ('lon'), zlib=True,complevel=4)
    fh_lon[:] = settings['lon']
    fh_lon.long_name = 'Longitude'
    fh_lon.units = 'Degrees East'

    fh_lat = fh.createVariable('lat', 'f4', ('lat'), zlib=True,complevel=4)
    fh_lat[:] = settings['lat']
    fh_lat.long_name = 'Longitude'
    fh_lat.units = 'Degrees North'

    fh_time = fh.createVariable('time', 'i2', ('time'), zlib=True,complevel=4)
    fh_time[:] = time
    fh_time.long_name = 'Time (value represents year average)'
    fh_time.units = 'Year'

    if pname == 'steric':
        mask = np.isnan(grid['mean'][0,...])
        mean_save = np.rint(20*grid['mean']).astype(int)
        mean_save[:,mask] = -20000
        sterr_save = np.rint(20*grid['sterr']).astype(int)
        sterr_save[:,mask] = -20000

        fh_grid = fh.createVariable(pname+'_mean', 'i2', ('time','lat','lon'), zlib=True,complevel=4)
        fh_grid[:] = mean_save
        fh_grid.missing_value = -20000
        fh_grid.scale_factor=0.05
        fh_grid.units = 'mm'
        fh_grid.long_name = name_dict[pname]+' mean'

        fh_grid = fh.createVariable(pname+'_sterr', 'i2', ('time','lat','lon'), zlib=True,complevel=4)
        fh_grid[:] = sterr_save
        fh_grid.missing_value = -20000
        fh_grid.scale_factor=0.05
        fh_grid.units = 'mm'
        fh_grid.long_name = name_dict[pname]+' standard error'
    else:
        mask = np.isnan(grid['rsl_mean'][0,...])
        rsl_mean_save = np.rint(20*grid['rsl_mean']).astype(int)
        rsl_sterr_save = np.rint(20*grid['rsl_sterr']).astype(int)
        rsl_mean_save[:,mask] = -20000
        rsl_sterr_save[:,mask] = -20000
        rad_mean_save = np.rint(20*grid['rad_mean']).astype(int)
        rad_sterr_save = np.rint(20*grid['rad_sterr']).astype(int)

        fh_grid = fh.createVariable(pname+'_rsl_mean', 'i2', ('time','lat','lon'), zlib=True,complevel=4)
        fh_grid[:] = rsl_mean_save
        fh_grid.missing_value = -20000
        fh_grid.scale_factor=0.05
        fh_grid.units = 'mm'
        fh_grid.long_name = name_dict[pname]+' relative sea level change'

        fh_grid = fh.createVariable(pname+'_rsl_sterr', 'i2', ('time','lat','lon'), zlib=True,complevel=4)
        fh_grid[:] = rsl_sterr_save
        fh_grid.missing_value = -20000
        fh_grid.scale_factor=0.05
        fh_grid.units = 'mm'
        fh_grid.long_name = name_dict[pname]+' relative sea level change standard error'

        fh_grid = fh.createVariable(pname+'_rad_mean', 'i2', ('time','lat','lon'), zlib=True,complevel=4)
        fh_grid[:] = rad_mean_save
        fh_grid.missing_value = -20000
        fh_grid.scale_factor=0.05
        fh_grid.units = 'mm'
        fh_grid.long_name = name_dict[pname]+' Solid-Earth deformation'

        fh_grid = fh.createVariable(pname+'_rad_sterr', 'i2', ('time','lat','lon'), zlib=True,complevel=4)
        fh_grid[:] = rad_sterr_save
        fh_grid.missing_value = -20000
        fh_grid.scale_factor=0.05
        fh_grid.units = 'mm'
        fh_grid.long_name = name_dict[pname]+' solid-Earth deformation standard error'
    fh.close()
    return
