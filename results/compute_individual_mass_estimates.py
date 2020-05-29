# -----------------------------------------------------
# Save the individual mass estimates, together with the
# averaged mass estimates used in this study for each
# process.
# -----------------------------------------------------
import numpy as np
import os
from netCDF4 import Dataset
import mod_gentools as gentools
import pandas as pd
import Levenshtein

def main():
    settings = {}
    if os.uname().nodename == 'MT-110180': settings['nproc'] = 4
    else: settings['nproc'] = 128
    settings['dir_data']    = os.getenv('HOME')+'/Data/'

    # Directories
    settings['dir_gwd_prep'] = settings['dir_data'] + 'Budget_20c/grd_prep/'
    settings['dir_glacier_zemp'] =  settings['dir_data'] + 'Glaciers/Zemp_2019/'
    settings['dir_gwd_20c']      = settings['dir_data'] + 'Hydrology/Groundwater/Past/'
    # Files
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    settings['fn_love'] = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    settings['fn_grace']         = settings['dir_data'] + 'Budget_20c/grd_prep/ewh_GRACE_annual_noGIA.nc'
    settings['fn_mascon_coords'] = settings['dir_data'] + 'GRACE/JPL_mascon/mascon_coords.npy'
    settings['fn_gia_ewh']       = settings['dir_data'] + 'GIA/Caron/Ensemble/ewh_ens_05.nc'
    settings['fn_marzeion_t'] = settings['dir_data']+'Glaciers/Marzeion_2015/data_for_thomas.txt'
    settings['fn_marzeion_p'] = settings['dir_data']+'Glaciers/Marzeion_2015/data_marzeion_etal_update_2015_regional.txt'
    settings['fn_parkes']     = settings['dir_data']+'Glaciers/Parkes2018/annual.csv'
    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'

    settings['fn_grd_glacier_regions'] = settings['dir_data'] + 'Budget_20c/grd_prep/grd_glacier_regions.nc'
    settings['fn_kjeldsen'] = settings['dir_data']+'IceSheets/Kjeldsen/Kjeldsen.csv'
    settings['fn_bamber']   = settings['dir_data']+'IceSheets/Bamber_2018/Bamber-etal_2018.tab'
    settings['fn_imbie']    = settings['dir_data']+'IceSheets/IMBIE/IMBIE2.txt'
    settings['fn_mouginot'] = settings['dir_data']+'IceSheets/Mouginot_2019/GrIS_total.csv'

    settings['fn_gwd_wada'] = settings['dir_data'] + 'Hydrology/Groundwater/waterdemand_30min_groundwaterdepletion_month_1960-2010.nc'
    settings['fn_gwd_doll_irr'] = settings['dir_data'] + 'Hydrology/Groundwater/Doll/Model_output/IRR70_S_TOTAL_WATER_STORAGES_mm_1960_2009.12.nc'
    settings['fn_gwd_doll_nouse'] = settings['dir_data'] + 'Hydrology/Groundwater/Doll/Model_output/NOUSE_S_TOTAL_WATER_STORAGES_mm_1960_2009.12.nc'

    settings['fn_chao_list'] = settings['dir_data'] + 'Hydrology/Dams/wu_res.xlsx'
    settings['fn_chao_loc']  = settings['dir_data'] + 'Hydrology/Dams/wu_locs.xlsx'
    settings['fn_lehner']    = settings['dir_data'] + 'Hydrology/Dams/dam_edit.xlsx'

    settings['fn_indiv'] = settings['dir_data'] + 'Budget_20c/results/indiv_mass.npy'

    settings['years']         = np.arange(1900,2019)
    settings['years_grace']   = np.arange(2003,2019)

    settings['num_ens'] = 100
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability'] / settings['probability'].sum()

    indiv_glaciers = compute_indiv_glaciers(settings)
    indiv_GrIS, indiv_AIS = compute_indiv_ice(indiv_glaciers, settings)
    indiv_tws = compute_indiv_tws(settings)

    indiv = {}
    indiv['tws'] = indiv_tws
    indiv['GrIS'] = indiv_GrIS
    indiv['AIS'] = indiv_AIS
    indiv['glac'] = indiv_glaciers
    np.save(settings['fn_indiv'],indiv)
    return

def compute_indiv_glaciers(settings):
    indiv_glaciers = {}
    marzeion_mean  = np.zeros([18, len(settings['years'])],dtype=np.float32)*np.nan
    # Process Marzeion
    data_p = np.loadtxt(settings['fn_marzeion_p'])
    glac_mn = -(data_p[:, 1:]).T
    for reg in range(18):
        mean_lcl = glac_mn[reg,:]
        marzeion_mean[reg,:114] = np.hstack((mean_lcl[0], mean_lcl[0], mean_lcl))
    indiv_glaciers['marzeion_GRL']  = np.cumsum(marzeion_mean[4,:])
    indiv_glaciers['marzeion_glac'] = np.cumsum(np.delete(marzeion_mean,4,axis=0).sum(axis=0))
    region_split = np.nanmean(marzeion_mean,axis=1) / np.nanmean(marzeion_mean,axis=1).sum()

    # Process Zemp
    zemp_mean  = np.zeros(marzeion_mean.shape) * np.nan
    zemp_flist = os.listdir(settings['dir_glacier_zemp'])
    for reg in range(18):
        fname = settings['dir_glacier_zemp']+[i for i in zemp_flist if '_'+str(reg+1)+'_' in i][0]
        raw_data = np.loadtxt(fname,skiprows=28,delimiter=',',usecols=(0,10,18))
        years      = raw_data[:,0]
        acc_a = np.in1d(settings['years'],years)
        acc_b = np.in1d(years,settings['years'])
        zemp_mean[reg,acc_a]  = -raw_data[acc_b,1]/362 # Rate in gigatons
    indiv_glaciers['zemp_GRL'] = np.zeros(len(settings['years']))*np.nan
    indiv_glaciers['zemp_glac'] = np.zeros(len(settings['years']))*np.nan
    indiv_glaciers['zemp_GRL'][62:-2] = np.cumsum(zemp_mean[4,62:-2])
    indiv_glaciers['zemp_glac'][62:-2]       = np.cumsum((np.delete(zemp_mean,4,axis=0).sum(axis=0))[62:-2])

    # Process Parkes
    data_parkes = np.loadtxt(settings['fn_parkes'], delimiter=',', skiprows=1)
    rate_parkes_global = np.zeros(len(settings['years']),dtype=np.float32)
    rate_parkes_global[2:-3] = -data_parkes[:,3]
    rate_parkes_global[:2] = -data_parkes[:10, 3].mean()
    parkes_global = (1.29/2)*np.cumsum(rate_parkes_global)
    indiv_glaciers['parkes_GRL'] = region_split[4] * parkes_global
    indiv_glaciers['parkes_glac'] = (1-region_split[4]) * parkes_global

    indiv_glaciers['marzeion_parkes_glac'] = indiv_glaciers['marzeion_glac'] + indiv_glaciers['parkes_glac']
    indiv_glaciers['zemp_parkes_glac']     = indiv_glaciers['zemp_glac'] + indiv_glaciers['parkes_glac']

    indiv_glaciers['marzeion_parkes_GRL'] = indiv_glaciers['marzeion_GRL'] + indiv_glaciers['parkes_GRL']
    indiv_glaciers['zemp_parkes_GRL']     = indiv_glaciers['zemp_GRL'] + indiv_glaciers['parkes_GRL']

    # Remove 2004
    for model in indiv_glaciers: indiv_glaciers[model]-=indiv_glaciers[model][103:108].mean()

    # Read GRACE and full
    full_ens = np.zeros([settings['num_ens'],len(settings['years'])])
    for ens in range(settings['num_ens']):
        full_ens[ens,:] = Dataset(settings['dir_data']+'Budget_20c/grd/grd_glac_'+str(ens)+'.nc','r').variables['barystatic'][:]._get_data()
    ts_full = ensemble_ts(full_ens, settings)
    indiv_glaciers['GRACE_glac'] = np.zeros(len(settings['years'])) * np.nan
    indiv_glaciers['GRACE_glac'][103:] = ts_full[103:,1]
    indiv_glaciers['full_glac'] = ts_full
    return(indiv_glaciers)

def compute_indiv_ice(indiv_glaciers,settings):
    # Read time series:
    # - AIS:  Bamber/IMBIE2/Adhikari
    # - GrIS: Bamber/Mouginot/Kjeldsen
    indiv_GrIS = {}
    indiv_AIS = {}

    imbie_raw    = np.loadtxt(settings['fn_imbie'], delimiter=';',usecols=(0,3,4))
    mouginot_raw = np.loadtxt(settings['fn_mouginot'], delimiter=';')
    kjeldsen_raw = np.loadtxt(settings['fn_kjeldsen'],delimiter=',')
    bamber_raw   = np.loadtxt(settings['fn_bamber'],delimiter='	',skiprows=21,usecols=(0,1,2,3,4,5,6))

    ### AIS
    # IMBIE
    imbie_yrs =np.arange(1992,2017)
    imbie_ann = np.zeros(len(imbie_yrs))
    for idx,yr in enumerate(imbie_yrs):
        yr_idx = np.floor(imbie_raw[:,0])==yr
        imbie_ann[idx] = imbie_raw[yr_idx,1].mean()
    indiv_AIS['IMBIE'] = np.zeros(len(settings['years'])) * np.nan
    indiv_AIS['IMBIE'][np.in1d(settings['years'],imbie_yrs)] = imbie_ann
    # Bamber
    indiv_AIS['Bamber'] = np.zeros(len(settings['years'])) * np.nan
    indiv_AIS['Bamber'][np.in1d(settings['years'],bamber_raw[:, 0])] = np.cumsum(bamber_raw[:, 3]+bamber_raw[:, 5])/-362
    for model in indiv_AIS: indiv_AIS[model]-=indiv_AIS[model][103:108].mean()

    ### GrIS
    pmg = indiv_glaciers['marzeion_parkes_GRL'].copy()
    pmg[-6:] = pmg[-6]

    indiv_GrIS['Bamber'] = np.zeros(len(settings['years'])) * np.nan
    indiv_GrIS['Bamber'][np.in1d(settings['years'],bamber_raw[:, 0])] = np.cumsum(bamber_raw[:, 1])/-362

    indiv_GrIS['Kjeldsen'] = np.zeros(len(settings['years'])) * np.nan
    indiv_GrIS['Kjeldsen'][:113] = -np.cumsum(kjeldsen_raw[60:,1])/362

    indiv_GrIS['Kjeldsen_P'] = np.zeros(len(settings['years'])) * np.nan
    indiv_GrIS['Kjeldsen_P'][:113] = -np.cumsum(kjeldsen_raw[60:,1])/362 + indiv_glaciers['marzeion_parkes_GRL'][:113]


    indiv_GrIS['Mouginot'] = np.zeros(len(settings['years'])) * np.nan
    indiv_GrIS['Mouginot'][72:] = -np.cumsum(mouginot_raw[:,1])/362

    indiv_GrIS['Mouginot_P'] = np.zeros(len(settings['years'])) * np.nan
    indiv_GrIS['Mouginot_P'][72:] = -np.cumsum(mouginot_raw[:,1])/362 + pmg[72:]

    for model in indiv_GrIS: indiv_GrIS[model]-=indiv_GrIS[model][103:108].mean()
    # for model in indiv_GrIS:
    #     plt.plot(settings['years'],indiv_GrIS[model])

    # Read GRACE and full
    full_ens_GrIS = np.zeros([settings['num_ens'],len(settings['years'])])
    full_ens_AIS = np.zeros([settings['num_ens'],len(settings['years'])])

    for ens in range(settings['num_ens']):
        full_ens_GrIS[ens,:] = Dataset(settings['dir_data']+'Budget_20c/grd/grd_GrIS_'+str(ens)+'.nc','r').variables['barystatic'][:]._get_data()
        full_ens_AIS[ens,:] = Dataset(settings['dir_data']+'Budget_20c/grd/grd_AIS_'+str(ens)+'.nc','r').variables['barystatic'][:]._get_data()

    ts_full_GrIS = ensemble_ts(full_ens_GrIS, settings)
    full_ens_AIS = ensemble_ts(full_ens_AIS, settings)

    indiv_GrIS['GRACE'] = np.zeros(len(settings['years'])) * np.nan
    indiv_AIS['GRACE'] = np.zeros(len(settings['years'])) * np.nan
    indiv_GrIS['GRACE'][103:] = ts_full_GrIS[103:,1]
    indiv_AIS['GRACE'][103:] = full_ens_AIS[103:,1]
    indiv_GrIS['full'] = ts_full_GrIS
    indiv_AIS['full'] = full_ens_AIS
    return(indiv_GrIS, indiv_AIS)

def compute_indiv_tws(settings):
    indiv_tws = {}

    # Natural: Read Humphrey average
    fn = settings['dir_data'] + 'Hydrology/Humphrey/01_monthly_grids_ensemble_means_allmodels/GRACE_REC_v03_JPL_GSWP3_monthly_ensemble_mean.nc'
    file_handle = Dataset(fn,'r')
    file_handle.set_auto_mask(False)
    lat = file_handle.variables['lat'][:]
    lon = file_handle.variables['lon'][:]
    time = file_handle.variables['time'][:]/365.25+1901
    load = file_handle.variables['rec_ensemble_mean'][:]
    file_handle.close()
    mask = np.load(settings['fn_mask'], allow_pickle=True).all()

    load = np.dstack([load[:, :, 360:], load[:, :, :360]])
    load[load<-32e4] = 0

    bary_annual = np.zeros(len(settings['years']))*np.nan

    area = gentools.grid_area(lat,lon)
    for idx, year in enumerate(settings['years']):
        acc_idx = np.floor(time).astype(int) == year
        if acc_idx.sum()>0:
            bary_annual[idx] = np.nansum(area*mask['tws']*load[acc_idx,:,:].mean(axis=0))
    indiv_tws['Humphrey'] = bary_annual/-362e12

    # GWD Wada/Doll
    wada_2016_scale = 0.6 # Scale factor (Wada 2016)
    wada_depletion = np.zeros(len(settings['years']))*np.nan
    # Old Wada 2010 data
    time_past = np.arange(1900,1960)
    for idx,yr in enumerate(time_past):
        data_raw = np.loadtxt(settings['dir_gwd_20c'] + 'gwd0' + str(yr) + '.asc', skiprows=6)
        data_raw[data_raw == -9999] = 0
        data_raw = np.flipud(data_raw)
        data_raw = np.hstack([data_raw[:, 360:], data_raw[:, 0:360]])
        wada_depletion[idx] = data_raw.sum()*1e9
    # Present (Wada 2014)
    file_handle = Dataset(settings['fn_gwd_wada'], 'r')
    file_handle.set_auto_mask(False)
    depletion_pd_monthly = file_handle.variables["anrg"][:]
    depletion_pd_monthly[depletion_pd_monthly < 0] = 0
    file_handle.close()
    time_pd_m = gentools.monthly_time(1960, 2010)
    time_present = np.arange(1960, 2011)
    # To normal grid
    depletion_rate_monthly = np.fliplr(depletion_pd_monthly)
    depletion_rate_monthly = np.dstack([depletion_rate_monthly[:, :, 360:], depletion_rate_monthly[:, :, 0:360]])
    for idx,yr in enumerate(time_present):
        idx_m = (time_pd_m>yr) & (time_pd_m<yr+1)
        wada_depletion[idx+60] = (np.sum(depletion_rate_monthly[idx_m, :, :], axis=0)).sum()*1e9
    indiv_tws['Wada'] = wada_2016_scale*np.cumsum(wada_depletion)/362e12


    ### DOLL
    file_handle = Dataset(settings['fn_gwd_doll_irr'], 'r')
    file_handle.set_auto_mask(False)
    time = 1960+file_handle.variables['time'][:]/12 + 1/24
    tws_irr = file_handle.variables['TWS_mm'][:]
    tws_irr[tws_irr < -9998] = np.nan
    file_handle.close()

    tws_nouse = Dataset(settings['fn_gwd_doll_nouse'], 'r').variables['TWS_mm'][:]._get_data()
    tws_nouse[tws_nouse < -9998] = np.nan

    # To 0-360
    tws_irr   = np.dstack([tws_irr[:,:,360:],tws_irr[:,:,:360]])
    tws_nouse = np.dstack([tws_nouse[:,:,360:],tws_nouse[:,:,:360]])

    # Depletion in irrigated minus natural
    tws_doll = tws_irr - tws_nouse

    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask_tws = (mask['land']) & (~mask['AIS']) & (~mask['GrIS'])

    doll_depletion = np.zeros(len(settings['years']))*np.nan

    # To annual
    tws_doll[:,~mask_tws] = 0
    tws_doll[np.isnan(tws_doll)] = 0
    time_ann_doll = np.arange(1960,2010,1)
    for idx, t in enumerate(time_ann_doll):
        acc_t = (np.floor(time).astype(int)==t)
        doll_depletion[idx+60] = (area*(tws_doll[acc_t,:,:].mean(axis=0))).sum()
    indiv_tws['Doll'] = doll_depletion/-362e12

    #################### Dams
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
    dam_years      = np.arange(1800, 2019, 1)

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

    indiv_tws['Chao'] = -total_storage[100:]/362e12
    for model in indiv_tws: indiv_tws[model] -= indiv_tws[model][103:108].mean()

    full_ens = np.zeros([settings['num_ens'],len(settings['years'])])
    for ens in range(settings['num_ens']):
        full_ens[ens,:] = Dataset(settings['dir_data']+'Budget_20c/grd/grd_tws_'+str(ens)+'.nc','r').variables['barystatic'][:]._get_data()
    ts_full = ensemble_ts(full_ens, settings)
    indiv_tws['GRACE'] = np.zeros(len(settings['years'])) * np.nan
    indiv_tws['GRACE'][103:] = ts_full[103:,1]
    indiv_tws['full'] = ts_full
    return(indiv_tws)



def ensemble_ts(ensemble,settings):
    tseries = np.zeros([len(settings['years']), 3])
    tseries[:, 1] = np.sum(settings['probability'][:, np.newaxis] * ensemble, axis=0)
    for t in range(len(settings['years'])):
        sort_idx = np.argsort(ensemble[:, t])
        sort_cdf = np.cumsum(settings['probability'][sort_idx])
        tseries[t, 0] = ensemble[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.05))]
        tseries[t, 2] = ensemble[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.95))]
    tseries-=tseries[103:108,1].mean()
    return(tseries)
