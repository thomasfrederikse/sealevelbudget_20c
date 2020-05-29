# --------------------------------------------------------------------------------------------------------
# Generate loads for each of the indidvidual TWS contributors:
# natural TWS variations, dam retention and groundwater depletion data
# Further processed to get grd fingerprints
# For each of the above:
# 1. Read data
# 2. Compute annual averages
# 3. Mascon representation
# 4. Save as netcdf
# --------------------------------------------------------------------------------------------------------

import numpy as np
import mod_gentools as gentools
from netCDF4 import Dataset
import os
import pandas as pd
import Levenshtein

def main():
    print('TWS loads...')
    settings = {}
    settings['time']    = np.arange(1900,2004)
    settings['dir_data']    = os.getenv('HOME') + '/Data/'

    # Directories
    settings['dir_gwd_20c']      = settings['dir_data'] + 'Hydrology/Groundwater/Past/'
    settings['dir_TWS_humphrey'] = settings['dir_data'] + 'Hydrology/Humphrey/02_monthly_grids_ensemble_model2/'
    settings['dir_gwd_prep']     = settings['dir_data'] + 'Budget_20c/grd_prep/'

    # Filenames
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    settings['fn_GRACE']   = settings['dir_data'] + 'GRACE/JPL_mascon/JPL_mascon_RL06v02_CRI_noGIA_noEQ_noseas.nc'
    settings['fn_mask_sl'] = settings['dir_data'] + 'GRACE/JPL_mascon/LAND_MASK.CRIv01.nc'
    settings['fn_mask_hy'] = settings['dir_data'] + 'GRACE/JPL_mascon/CLM4.SCALE_FACTOR.JPL.MSCNv01CRIv01.nc'
    settings['fn_mascon_coords'] = settings['dir_data'] + 'GRACE/JPL_mascon/mascon_coords.npy'
    settings['fn_gwd_wada'] = settings['dir_data'] + 'Hydrology/Groundwater/waterdemand_30min_groundwaterdepletion_month_1960-2010.nc'
    settings['fn_gwd_doll_irr'] = settings['dir_data'] + 'Hydrology/Groundwater/Doll/Model_output/IRR70_S_TOTAL_WATER_STORAGES_mm_1960_2009.12.nc'
    settings['fn_gwd_doll_nouse'] = settings['dir_data'] + 'Hydrology/Groundwater/Doll/Model_output/NOUSE_S_TOTAL_WATER_STORAGES_mm_1960_2009.12.nc'

    settings['fn_chao_list'] = settings['dir_data'] + 'Hydrology/Dams/wu_res.xlsx'
    settings['fn_chao_loc']  = settings['dir_data'] + 'Hydrology/Dams/wu_locs.xlsx'
    settings['fn_lehner']    = settings['dir_data'] + 'Hydrology/Dams/dam_edit.xlsx'
    settings['fn_save'] = settings['dir_gwd_prep'] + 'load_tws_dam_gwd.nc'

    GRACE = GRACE_data(settings)
    humphrey_annual(settings)
    depletion_mscn = prep_GWD(GRACE,settings)
    depletion_doll_mscn = prep_GWD_Doll(GRACE, depletion_mscn, settings)
    dam_load_mscn = prep_dams(GRACE, settings)
    save_loads(GRACE, depletion_mscn, depletion_doll_mscn, dam_load_mscn, settings)
    return

def humphrey_annual(settings):
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    # Read monthly data and save annual-means 1900-2004
    for ens in range(100):
        print(ens)
        fn = settings['dir_TWS_humphrey'] + 'GRACE_REC_v03beta_JPL_GSWP3_monthly_ens'+str(ens+1).zfill(3)+'.nc'
        file_handle = Dataset(fn)
        file_handle.set_auto_mask(False)
        lat = file_handle.variables['lat'][:]
        lon = file_handle.variables['lon'][:]
        time = gentools.monthly_time(1901,2004)
        tws = file_handle.variables['rec_ensemble_member'][:len(time),:,:]
        tws[tws<-10000]=0
        file_handle.close()
        h_year = np.arange(1900,2004)

        # Annual means
        humphrey = np.zeros([len(h_year),len(lat),len(lon)])
        for idx,yr in enumerate(h_year[1:]):
            t_idx = (time > yr) & (time < yr+1)
            humphrey[idx+1,:,:] = np.mean(tws[t_idx,:,:],axis=0)

        # Flip array to longitudes 0-360
        humphrey = np.dstack([humphrey[:,:,360:],humphrey[:,:,:360]])

        # Mask out non-tws regions
        humphrey = humphrey * mask['tws'][np.newaxis,:,:]

        fn = settings['dir_gwd_prep'] + 'load_tws_natural_' + str(ens) + '.nc'
        file_handle = Dataset(fn, 'w')
        file_handle.createDimension('x', len(lon))
        file_handle.createDimension('y', len(lat))
        file_handle.createDimension('t', len(h_year))
        file_handle.createVariable('x', 'f4', ('x',), zlib=True)[:] = np.arange(0.25,360.25,0.5)
        file_handle.createVariable('y', 'f4', ('y',), zlib=True)[:] = lat
        file_handle.createVariable('t', 'i4', ('t',), zlib=True)[:] = h_year
        file_handle.createVariable('z', 'i2', ('t', 'y', 'x',), zlib=True, complevel=4)[:] = humphrey
        file_handle.close()
    return

def save_loads(GRACE, depletion_mscn, depletion_doll_mscn, dam_load_mscn, settings):
    # Save dam retention and GWD data in netCDF file
    file_handle = Dataset(settings['fn_save'], 'w')
    file_handle.createDimension('x', len(GRACE['lon']))
    file_handle.createDimension('y', len(GRACE['lat']))
    file_handle.createDimension('t', 104)
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = GRACE['lon']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = GRACE['lat']
    file_handle.createVariable('t', 'i4', ('t',),zlib=True)[:] = np.arange(1900,2004)
    file_handle.createVariable('GWD_wada', 'i2', ('t', 'y', 'x',),zlib=True)[:] = depletion_mscn
    file_handle.createVariable('GWD_doll', 'i2', ('t', 'y', 'x',),zlib=True)[:] = depletion_doll_mscn
    file_handle.createVariable('Dams', 'i2', ('t', 'y', 'x',),zlib=True)[:] = dam_load_mscn
    file_handle.close()
    return

def prep_dams(GRACE,settings):
    # --------------------------------------------------------------------------------------------
    # Prepare dam water storage data
    # 1. Read dam lists from the GRanD list and the data from Ben Chao
    # 2. Remove duplicates as good as it gets
    # 3. Fill global grid with dam retention data and seepage estimates based on Chao et al., 2008
    # --------------------------------------------------------------------------------------------
    print('   Reading spreadsheets...')
    chao_list_raw = pd.read_excel(settings['fn_chao_list'], header=None).values  # Read Chao list
    chao_loc_raw  = pd.read_excel(settings['fn_chao_loc'], header=None).values   # Read Chao locs
    lehner_raw    = pd.read_excel(settings['fn_lehner'], header=None).values     # Read Lehner locs

    # name lat lon year cap loc_avail
    dam_list_grand = np.zeros([len(lehner_raw), 6], dtype=object)
    print('   Processing GRanD list...')
    lehner_raw[:, 2][lehner_raw[:, 2]<0] = lehner_raw[:, 2][lehner_raw[:, 2]<0] + 360
    for i in range(len(lehner_raw)):
        dam_list_grand[i, 0] = lehner_raw[i, 1].split('(')[0].lower().replace(" ", "")
        dam_list_grand[i, 1] = lehner_raw[i, 3]
        dam_list_grand[i, 2] = lehner_raw[i, 2]
        dam_list_grand[i, 3] = lehner_raw[i, 7]
        dam_list_grand[i, 4] = lehner_raw[i, 6] * 1000
        dam_list_grand[i, 5] = True
    dam_list_grand = dam_list_grand[dam_list_grand[:, 4] > 0, :]

    print('   Processing Ben Chaos location list...')
    # --------------------------------------------------------------
    # GRanD list and Benjamin Chao's location list have overlap
    # Challenge: find dams that are in Chao's list, but not in GRanD
    # These dams have a known lat/lon location
    # --------------------------------------------------------------
    dam_list_chao_loc = []
    chao_loc_raw[:, 6][chao_loc_raw[:, 6]<0] = chao_loc_raw[:, 6][chao_loc_raw[:, 6]<0] + 360
    for i in range(len(chao_loc_raw)):
        dam_list_lcl = np.zeros(6, dtype=object)
        dam_list_lcl[0] = chao_loc_raw[i, 2].split('(')[0].lower().replace(" ", "")
        diff_score = np.zeros(len(dam_list_grand))
        for lcl in range(len(diff_score)):
            diff_score[lcl] = Levenshtein.ratio(dam_list_lcl[0], dam_list_grand[lcl, 0])
            vol_ratio = chao_loc_raw[i, 3] / (dam_list_grand[lcl, 4])
            if (vol_ratio < 0.5) | (vol_ratio > 2): diff_score[lcl] = -1                     # If capacity is vastly different, probably different dam
            if np.abs(chao_loc_raw[i, 4] - dam_list_grand[lcl, 3]) > 2: diff_score[lcl] = -1 # If year of construction is vastly different, probably different dam
        if np.max(diff_score) < 0.9:
            dam_list_lcl[1] = chao_loc_raw[i, 5]
            dam_list_lcl[2] = chao_loc_raw[i, 6]
            dam_list_lcl[3] = chao_loc_raw[i, 4]
            dam_list_lcl[4] = chao_loc_raw[i, 3]
            dam_list_lcl[5] = True
            dam_list_chao_loc.append(dam_list_lcl)
    dam_list_chao_loc = np.array(dam_list_chao_loc)
    dam_list_loc = np.vstack([dam_list_grand, dam_list_chao_loc])

    # --------------------------------------------------------------
    # GRanD list and Benjamin Chao's location list have overlap
    # Challenge: find dams that are in Chao's list, but not in GRanD
    # These dams don't have a known lat/lon location
    # --------------------------------------------------------------
    print('   Processing Ben Chaos full list...')
    dam_list_chao_full = []
    for i in range(len(chao_list_raw)):
        if chao_list_raw[i, 3] > 1000:
            dam_list_lcl = np.zeros(6, dtype=object)
            dam_list_lcl[0] = chao_list_raw[i, 2].split('(')[0].lower().replace(" ", "")
            diff_score = np.zeros(len(dam_list_loc))
            for lcl in range(len(dam_list_loc)):
                diff_score[lcl] = Levenshtein.ratio(dam_list_lcl[0], dam_list_loc[lcl, 0])
            if np.max(diff_score) < 0.7:
                dam_list_lcl[3] = chao_list_raw[i, 4]
                dam_list_lcl[4] = chao_list_raw[i, 3]
                dam_list_lcl[5] = False
                dam_list_chao_full.append(dam_list_lcl)
    dam_list_chao_full = np.array(dam_list_chao_full)
    dam_list = np.vstack([dam_list_grand, dam_list_chao_loc,dam_list_chao_full])

    # From 1000 m3 to kg
    dam_list[:, 4] = dam_list[:, 4] * 1e6
    dam_years      = np.arange(1800, 2004, 1)

    # Index for dams with known locations
    nlocs = np.sum(dam_list[:, 5])
    dam_load_list = np.zeros([nlocs, len(dam_years)])
    no_seepage = ['manicouagan', 'jenpeg', 'smallwoodreservoir', 'missifallscontrol', 'earfalls', 'whitesandrapids', 'pipmuacan', 'keenleyside', 'sanhezha', 'tainionkoski', 'irkutsk', 'verkhnetulomskaya', 'ondakumskaya', 'verkhnesvirskaya', 'structure308']

    # Total volumes
    print('   Computing storage and seepage...')
    total_volume  = np.zeros(len(dam_years)) # Total water in dam assuming 85% full
    total_scaled  = np.zeros(len(dam_years)) # Total water in dam assuming 85% full
    total_seepage = np.zeros(len(dam_years)) # Total TWS due to seepage after dam construction
    total_storage = np.zeros(len(dam_years)) # Total TWS due to water in dam and seepage
    for i in range(len(dam_list)):
        local_volume  = np.zeros(len(dam_years))
        local_seepage = np.zeros(len(dam_years))
        if (dam_list[i, 3] > 1799) & (dam_list[i, 3] < 2004):
            startindex = int(dam_list[i, 3] - 1800)
        else:
            startindex = 0
        local_volume[startindex:] = dam_list[i, 4] * 0.85
        seepage_growth = np.minimum(np.cumsum(1 / np.sqrt(dam_years[startindex + 1:] - (dam_years[startindex]))), 20)
        local_seepage[startindex + 1:] = dam_list[i, 4] * 0.05 * seepage_growth
        total_scaled = total_scaled + local_volume
        total_seepage = total_seepage + local_seepage
        if dam_list[i, 0] in no_seepage:
            local_storage = local_volume
        else:
            local_storage = local_seepage + local_volume
        total_storage = total_storage + local_storage
        total_volume = total_volume + local_volume/0.85
        # If location known, add as new entry, otherwise spread over all other locations
        if i < nlocs:
            dam_load_list[i, :] = local_storage
        else:
            dam_load_list = dam_load_list + (local_storage / nlocs)[np.newaxis, :]

    # List with dam indices
    print('   Adding all dams to grid...')
    area = gentools.grid_area(GRACE['lat'],GRACE['lon'])
    dam_load = np.zeros([len(dam_years),len(GRACE['lat']),len(GRACE['lon'])])
    for i in range(nlocs):
        idx_lat = np.argmin(np.abs(GRACE['lat'] - dam_list[i, 1]))
        idx_lon = np.argmin(np.abs(GRACE['lon'] - dam_list[i, 2]))
        dam_load[:, idx_lat, idx_lon] = dam_load[:, idx_lat, idx_lon] + dam_load_list[i, :] / area[idx_lat, idx_lon]

    # Only 1900-2004
    acc_idx = np.in1d(dam_years, settings['time'])
    dam_load = dam_load[acc_idx, :, :]
    total_scaled = total_scaled[acc_idx]
    total_seepage = total_seepage[acc_idx]
    total_storage = total_storage[acc_idx]
    total_volume = total_volume[acc_idx]

    # smooth_trend_s = np.zeros(len(settings['time']))*np.nan
    #
    # flen=30
    # fhalf = int(flen/2)
    # for idx,yr in enumerate(settings['time']):
    #     if (idx>(fhalf-1)) & (idx<len(settings['time'])-fhalf-1):
    #         amat = np.ones([flen,2])
    #         amat[:,1] = settings['time'][idx-fhalf:idx+fhalf]
    #         smooth_trend_s[idx] = np.linalg.lstsq(amat,total_storage[idx-fhalf:idx+fhalf],rcond=None)[0][1]


    # To mascon
    print('   Grid to mascon...')
    dam_load_mscn = np.zeros(dam_load.shape)
    # Into mascons
    for k in range(len(GRACE['mscn_coords'])):
        lat_acc = np.where((GRACE['lat'] >= GRACE['mscn_coords'][k, 0]) & (GRACE['lat'] < GRACE['mscn_coords'][k, 1]))[0]
        lon_acc = np.where((GRACE['lon'] >= GRACE['mscn_coords'][k, 2]) & (GRACE['lon'] < GRACE['mscn_coords'][k, 3]))[0]
        lsm_acc = GRACE['land'][lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1]
        dp_mean = np.mean(dam_load[:,lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1]*lsm_acc[np.newaxis,:,:],axis=(1,2))
        dam_load_mscn[:,lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1] = dp_mean[:,np.newaxis,np.newaxis]
    dam_load_mscn = dam_load_mscn * GRACE['land'][np.newaxis,:,:]
    dam_load_global_mscn = np.nansum(dam_load_mscn * area, axis=(1, 2))
    dam_load_mscn = dam_load_mscn * (total_storage/dam_load_global_mscn)[:,np.newaxis,np.newaxis]
    return(dam_load_mscn)

def prep_GWD(GRACE,settings):
    # ------------------------------------------------------------------
    # Groundwater depletion
    #  Use estimates from Wada et al., 2010, 2014, which cover 1900-2010
    #  Merge the 2 data sets
    #  Sum depletion rate to get total depletion
    #  Wada 2016: not all depleted water ends up in ocean. We apply a
    #  correction factor of 0.6 to account for this effect.
    # ------------------------------------------------------------------
    wada_2016_scale = 0.6 # Scale factor (Wada 2016)

    area = gentools.grid_area(GRACE['lat'], GRACE['lon'])
    time_tot = np.arange(1900,2004)
    depletion = np.zeros([len(time_tot),len(GRACE['lat']),len(GRACE['lon'])])
    # Old Wada 2010 data
    time_past = np.arange(1900,1960)
    for idx,yr in enumerate(time_past):
        data_raw = np.loadtxt(settings['dir_gwd_20c'] + 'gwd0' + str(yr) + '.asc', skiprows=6)
        data_raw[data_raw == -9999] = 0
        data_raw = np.flipud(data_raw)
        data_raw = np.hstack([data_raw[:, 360:], data_raw[:, 0:360]])
        depletion[idx, :, :] = data_raw
    # Present (Wada 2014)
    file_handle = Dataset(settings['fn_gwd_wada'], 'r')
    file_handle.set_auto_mask(False)
    depletion_pd_monthly = file_handle.variables["anrg"][:]
    depletion_pd_monthly[depletion_pd_monthly < 0] = 0
    file_handle.close()
    time_pd_m = gentools.monthly_time(1960, 2010)
    time_present = np.arange(1960, 2004)
    # To normal grid
    depletion_rate_monthly = np.fliplr(depletion_pd_monthly)
    depletion_rate_monthly = np.dstack([depletion_rate_monthly[:, :, 360:], depletion_rate_monthly[:, :, 0:360]])
    for idx,yr in enumerate(time_present):
        idx_m = (time_pd_m>yr) & (time_pd_m<yr+1)
        depletion[idx+60, :, :] = np.sum(depletion_rate_monthly[idx_m, :, :], axis=0)
    depletion = depletion/ area  * -1e9
    depletion = np.cumsum(depletion, axis=0)* wada_2016_scale

    # Test global
    depletion_global = np.nansum(depletion * area, axis=(1, 2))
    depletion_mscn = np.zeros(depletion.shape)
    # Into mascons
    for k in range(len(GRACE['mscn_coords'])):
        lat_acc = np.where((GRACE['lat'] >= GRACE['mscn_coords'][k, 0]) & (GRACE['lat'] < GRACE['mscn_coords'][k, 1]))[0]
        lon_acc = np.where((GRACE['lon'] >= GRACE['mscn_coords'][k, 2]) & (GRACE['lon'] < GRACE['mscn_coords'][k, 3]))[0]
        lsm_acc = GRACE['land'][lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1]
        dp_mean = np.mean(depletion[:,lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1]*lsm_acc[np.newaxis,:,:],axis=(1,2))
        depletion_mscn[:,lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1] = dp_mean[:,np.newaxis,np.newaxis]
    depletion_mscn = depletion_mscn * GRACE['land'][np.newaxis,:,:]
    depletion_global_mscn = np.nansum(depletion_mscn * area, axis=(1, 2))
    depletion_mscn = depletion_mscn * (depletion_global/depletion_global_mscn)[:,np.newaxis,np.newaxis]
    return(depletion_mscn)

def prep_GWD_Doll(GRACE, depletion_mscn, settings):
    # ---------------------------------------------
    # Ground water depletion from Doll et al (2014)
    # - Transform from monthly to annual
    # - 0-360
    # - Mask out Greenland
    # ---------------------------------------------
    file_handle = Dataset(settings['fn_gwd_doll_irr'], 'r')
    file_handle.set_auto_mask(False)
    time = 1960+file_handle.variables['time'][:]/12 + 1/24
    lat = file_handle.variables['lat'][:]
    lon = file_handle.variables['lon'][:]
    tws_irr = file_handle.variables['TWS_mm'][:]
    tws_irr[tws_irr < -9998] = np.nan
    file_handle.close()

    tws_nouse = Dataset(settings['fn_gwd_doll_nouse'], 'r').variables['TWS_mm'][:]._get_data()
    tws_nouse[tws_nouse < -9998] = np.nan

    # To 0-360
    tws_irr   = np.dstack([tws_irr[:,:,360:],tws_irr[:,:,:360]])
    tws_nouse = np.dstack([tws_nouse[:,:,360:],tws_nouse[:,:,:360]])

    # Depletion in irrigated minus natural
    tws = tws_irr - tws_nouse

    area = gentools.grid_area(GRACE['lat'], GRACE['lon'])

    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask_tws = (mask['land']) & (~mask['AIS']) & (~mask['GrIS'])

    # To annual
    tws[:,~mask_tws] = 0
    tws[np.isnan(tws)] = 0
    time_ann_doll = np.arange(1960,2004,1)
    tws_ann_doll = np.zeros([len(time_ann_doll),len(lat),len(lon)])
    for idx, t in enumerate(time_ann_doll):
        acc_t = (np.floor(time).astype(int)==t)
        tws_ann_doll[idx,:,:] = tws[acc_t,:,:].mean(axis=0)

    # To mascon
    tws_ann_doll_global = np.nansum(tws_ann_doll * area, axis=(1, 2))
    tws_ann_doll_mscn = np.zeros(tws_ann_doll.shape)
    # Into mascons
    for k in range(len(GRACE['mscn_coords'])):
        lat_acc = np.where((GRACE['lat'] >= GRACE['mscn_coords'][k, 0]) & (GRACE['lat'] < GRACE['mscn_coords'][k, 1]))[0]
        lon_acc = np.where((GRACE['lon'] >= GRACE['mscn_coords'][k, 2]) & (GRACE['lon'] < GRACE['mscn_coords'][k, 3]))[0]
        lsm_acc = GRACE['land'][lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1]
        dp_mean = np.mean(tws_ann_doll[:,lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1]*lsm_acc[np.newaxis,:,:],axis=(1,2))
        tws_ann_doll_mscn[:,lat_acc[0]:lat_acc[-1] + 1, lon_acc[0]:lon_acc[-1] + 1] = dp_mean[:,np.newaxis,np.newaxis]
    tws_ann_doll_mscn = tws_ann_doll_mscn * GRACE['land'][np.newaxis,:,:]
    tws_ann_doll_global_mscn = np.nansum(tws_ann_doll_mscn * area, axis=(1, 2))
    tws_ann_doll_mscn *= (tws_ann_doll_global/tws_ann_doll_global_mscn)[:,np.newaxis,np.newaxis]

    depletion_doll_mscn = np.zeros(depletion_mscn.shape)
    depletion_doll_mscn[:61,...] = depletion_mscn[:61,...]
    depletion_doll_mscn[60:,...] = tws_ann_doll_mscn - tws_ann_doll_mscn[0,...][np.newaxis,...] + depletion_doll_mscn[60,...][np.newaxis,...]
    # tst = np.nansum(tws_annual * area, axis=(1, 2))
    # tst2 = np.nansum(depletion_mscn * area, axis=(1, 2))
    # tst3 = np.nansum(depletion_doll_mscn * area, axis=(1, 2))
    return(depletion_doll_mscn)

def GRACE_data(settings):
    GRACE = {}
    file_handle = Dataset(settings['fn_GRACE'])
    file_handle.set_auto_mask(False)
    GRACE['lat'] = file_handle.variables['lat'][:]
    GRACE['lon'] = file_handle.variables['lon'][:]
    GRACE['land']  = Dataset(settings['fn_mask_sl'], 'r').variables["land_mask"][:]._get_data().astype(bool)
    GRACE['mscn_coords'] = np.load(settings['fn_mascon_coords'],allow_pickle=True) # Coordinates of mascons
    file_handle.close()
    return(GRACE)
