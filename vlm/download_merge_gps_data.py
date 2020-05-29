# ---------------------------------------
# Download the ITRF2014 GPS data from UNR
# Merge with ITRF2008 data for stations
# that have disappeared from the database
# ---------------------------------------
# Save in numpy arrays
# lat lon name tseries jumps
# ---------------------------------------
import urllib.request
import numpy as np
import os
from joblib import Parallel,delayed

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_ITRF2014'] = settings['dir_data'] + 'GPS/UNR/ITRF2014_2020_02_12/'
    settings['dir_ITRF2008'] = settings['dir_data'] + 'GPS/UNR/ITRF2008_2018_06_15/'
    settings['min_GPS_obs'] = 4 * 365
    settings['fn_gps_tseries'] = settings['dir_data'] + 'Budget_20c/vlm/gps_tseries.npy'

    download_ITRF2014(settings)
    ts_ITRF2014 = read_data(2014, settings)
    ts_ITRF2008 = read_data(2008, settings)
    ts_combined = combine_datasets(ts_ITRF2008,ts_ITRF2014,settings)
    np.save(settings['fn_gps_tseries'], ts_combined)
    return

def combine_datasets(ts_ITRF2008,ts_ITRF2014,settings):
    list_2014 = []
    for i in range(len(ts_ITRF2014)):
        list_2014.append(ts_ITRF2014[i]['name'])
    list_2008 = []
    for i in range(len(ts_ITRF2008)):
        list_2008.append(ts_ITRF2008[i]['name'])

    list_combined = []
    for stat in range(len(ts_ITRF2014)):
        if list_2014[stat] not in list_2008: # In 2014 not in 2008, append
            list_combined.append(ts_ITRF2014[stat])
        else: # In both 2014 and 2008
            idx_2008 = list_2008.index(list_2014[stat])
            if ts_ITRF2014[stat]['n_obs'] >= 0.75*ts_ITRF2008[idx_2008]['n_obs']:
                list_combined.append(ts_ITRF2014[stat])
            else:
                print('ITRF2008 longer: '+list_2014[stat])
                list_combined.append(ts_ITRF2008[idx_2008])
    for stat in range(len(ts_ITRF2008)):
        if list_2008[stat] not in list_2014:
            print('Not in ITRF2014: '+list_2008[stat])
            list_combined.append(ts_ITRF2008[stat])

    # Sort list
    list_names = []
    for i in range(len(list_combined)):
        list_names.append(list_combined[i]['name'])
    list_sort = np.argsort(list_names)
    ts_combined = np.zeros(len(list_sort),dtype=object)
    for i in range(len(ts_combined)):
        ts_combined[i] = list_combined[list_sort[i]]
    return(ts_combined)

def download_ITRF2014(settings):
    # 1. Make list of downloads
    print('   Saving download list...')
    save_statlist = settings['dir_ITRF2014']+'statlist.txt'
    null = urllib.request.urlretrieve('http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt',save_statlist)

    # Read filelist
    statlist = np.loadtxt(save_statlist,skiprows=1,usecols=(0,1,2,10),dtype='object')
    statlist[:,1:3] = statlist[:,1:3].astype(float)
    statlist[:,3] = statlist[:,3].astype(int)

    # Remove short records
    acc_len_idx = statlist[:,3]>settings['min_GPS_obs']
    statlist = statlist[acc_len_idx,:]
    station_list = np.zeros(len(statlist),dtype=object)
    for i in range(len(station_list)):
        station_list[i] = {}
        station_list[i]['name'] = statlist[i,0]
        station_list[i]['lat'] = statlist[i,1]
        station_list[i]['lon'] = statlist[i,2]
        station_list[i]['url'] = 'http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/'+station_list[i]['name']+'.tenv3'
    Parallel(n_jobs=8)(delayed(download_gps_indiv)(station_list[i]['url'],settings) for i in range(len(station_list)))

    # steplist list
    save_steplist = settings['dir_ITRF2014']+'steps.txt'
    save_decyear = settings['dir_ITRF2014']+'decyear.txt'
    void = urllib.request.urlretrieve('http://geodesy.unr.edu/NGLStationPages/steps.txt', save_steplist)
    void = urllib.request.urlretrieve('http://geodesy.unr.edu/NGLStationPages/decyr.txt', save_decyear)

    # name lat lon list
    np.savetxt(save_statlist,statlist[:,:3],'%s %.3f %.3f')
    return

def download_gps_indiv(url,settings):
    savename = settings['dir_ITRF2014'] + url[-10:]
    print('   ' + savename)
    null = urllib.request.urlretrieve(url, savename)  # Download
    return

def read_data(run,settings):
    print('   Processing ITRF'+str(run))
    statlist = np.loadtxt(settings['dir_ITRF'+str(run)]+'statlist.txt',dtype=object)
    statlist[:,1:3] = statlist[:,1:3].astype(float)
    ts_data = np.zeros(len(statlist),dtype='object')

    # Prepare step list
    fn_steplist = settings['dir_ITRF'+str(run)]+'steplist.txt'
    fn_decyear  = settings['dir_ITRF'+str(run)]+'decyear.txt'
    steplist     = np.loadtxt(fn_steplist, delimiter='  ', dtype=object, usecols=(0,1,2))
    steplist[:,2] = steplist[:,2].astype(int)
    decyr_lookup = np.loadtxt(fn_decyear, delimiter=' ', dtype=object, usecols=(0,1))
    decyr_dict = dict(zip(decyr_lookup[:,0],decyr_lookup[:,1].astype(float)))
    step_stats = steplist[steplist[:,2]==1,0]
    step_times = steplist[steplist[:,2]==1,1]
    step_decyears = np.zeros(len(step_times))
    for i in range(len(step_times)):
        step_decyears[i] = decyr_dict[step_times[i]]
    step_list = np.vstack([step_stats,step_decyears]).T

    stat_data_raw = Parallel(n_jobs=8)(delayed(read_data_indiv)(stat,step_list,run,statlist,settings) for stat in range(len(statlist)))
    for stat in range(len(statlist)):
        ts_data[stat] = stat_data_raw[stat]
    return(ts_data)

def read_data_indiv(stat,step_list,run,statlist,settings):
    print(str(stat) + '/' + str(len(statlist)))
    stat_indiv = {}
    stat_indiv['name'] = statlist[stat, 0]
    stat_indiv['lat'] = statlist[stat, 1]
    stat_indiv['lon'] = statlist[stat, 2]
    # Read file and store height
    if run==2008: fn = settings['dir_ITRF' + str(run)] + statlist[stat, 0] + '.tenv'
    elif run==2014: fn = settings['dir_ITRF' + str(run)] + statlist[stat, 0] + '.tenv3'
    else: print('ERROR')
    data_raw = np.loadtxt(fn, usecols=(2, 12), skiprows=1)
    stat_indiv['time'] = data_raw[:, 0]
    stat_indiv['height'] = data_raw[:, 1]
    stat_indiv['n_obs'] = len(stat_indiv['height'])
    stat_indiv['ref_frame'] = 'ITRF'+str(run)
    # Jumps
    acc_steps = (step_list[:, 0] == statlist[stat, 0])
    if acc_steps.sum() > 0:  # There are jumps
        stat_indiv['has_jumps'] = True
        stat_indiv['jumps'] = step_list[acc_steps, 1]
    else:
        stat_indiv['has_jumps'] = False
    return(stat_indiv)
