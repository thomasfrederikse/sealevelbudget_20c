# -----------------------------------------
# Toolbox for QC'ing all tide-gauge regions
# Also contains list of final station
# selection per ocean basin
# -----------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import os
import mod_gentools as gentools
from netCDF4 import Dataset
import datetime as dt


def main():
    date_now = dt.datetime.now()
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_regions_for_selection'] = settings['dir_budget'] + 'tg/regions_for_selection.npy'
    settings['fn_alttg_for_region_selection'] = settings['dir_budget'] + 'vlm/alttg_for_region_selection.npy'
    settings['fn_gps_data'] = settings['dir_budget'] + 'vlm/gps_data.npy'
    settings['fn_region_list'] = settings['dir_budget'] + 'region_data/region_list_' + str(date_now.year) + '_' + str(date_now.month).zfill(2) + '_' + str(date_now.day).zfill(2) + '_' + str(date_now.hour).zfill(2) + '_' + str(date_now.minute).zfill(2) + '.npy'
    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['years'] = np.arange(1900, 2019)
    settings['num_ens'] = 100
    settings['min_years_reg'] = 10

    # Write region list file from selected tide gauge regions
    region_list = np.zeros(6, dtype=object)
    for reg in range(len(region_list)): region_list[reg] = {}
    region_list[0]['name'] = 'Subpolar North Atlantic'
    region_list[0]['list'] = create_SPNA_list()
    region_list[1]['name'] = 'Indian Ocean - South Pacific'
    region_list[1]['list'] = create_IOSP_list()
    region_list[2]['name'] = 'Subtropical North Atlantic'
    region_list[2]['list'] = create_STNA_list()
    region_list[3]['name'] = 'East Pacific'
    region_list[3]['list'] = create_EAPA_list()
    region_list[4]['name'] = 'South Atlantic'
    region_list[4]['list'] = create_SATL_list()
    region_list[5]['name'] = 'Northwest Pacific'
    region_list[5]['list'] = create_NWPA_list()
    np.save(settings['fn_region_list'], region_list)
    return

def sel_stats(settings):
    # Loop over all stations to perform QC and select stations
    regions_for_selection = np.load(settings['fn_regions_for_selection'], allow_pickle=True).all()
    gps_data = np.load(settings['fn_gps_data'], allow_pickle=True)
    alttg_data = np.load(settings['fn_alttg_for_region_selection'], allow_pickle=True).all()
    probability = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    probability = probability / probability.sum()
    gps_dict = {}
    for region in range(len(gps_data)):
        if gps_data[region]['has_gps']:
            for station in range(len(gps_data[region]['id'])):
                gps_dict[gps_data[region]['id'][station]] = region
    alttg_dict = {}
    for region in range(len(alttg_data['id'])):
        for station in range(len(alttg_data['id'][region])):
            alttg_dict[alttg_data['id'][region][station]] = region

    basin_num = 3  # Set number for which region selection must be watched
    for region in range(len(regions_for_selection['id'])):
        if (regions_for_selection['basin_num'][region] == basin_num) & (np.isfinite(regions_for_selection['height_corr'][region, :]).sum() > settings['min_years_reg']):
            vlm_name = []
            vlm_est = []
            vlm_err = []
            if regions_for_selection['id'][region][0] in gps_dict:
                vlm_name = vlm_name + gps_data[gps_dict[regions_for_selection['id'][region][0]]]['gps_name']
                vlm_est = vlm_est + gps_data[gps_dict[regions_for_selection['id'][region][0]]]['resvlm_trend_mean'][:,0].tolist()
                vlm_err = vlm_err + gps_data[gps_dict[regions_for_selection['id'][region][0]]]['vlm_trend_mean'][:,1].tolist()
            if regions_for_selection['id'][region][0] in alttg_dict:
                vlm_name.append(alttg_data['code'][alttg_dict[regions_for_selection['id'][region][0]]])
                vlm_est.append(alttg_data['resvlm_trend_mean'][alttg_dict[regions_for_selection['id'][region][0]], 0])
                vlm_err.append(alttg_data['resvlm_sterr_AR1'][alttg_dict[regions_for_selection['id'][region][0]]])
            vlm_acc = np.array(vlm_err) < 1.02
            vlm_est = np.array(vlm_est)[vlm_acc]
            vlm_err = np.array(vlm_err)[vlm_acc]
            vlm_name = np.array(vlm_name)[vlm_acc]
            plt.plot(settings['years'], regions_for_selection['height_corr'][region, :] - regions_for_selection['rsl_grd_dev'][region] - (settings['years'] * regions_for_selection['rsl_gia_dev'][region]),'o',color='black')
            plt.plot(settings['years'], regions_for_selection['height_corr'][region, :] - regions_for_selection['rsl_grd_dev'][region] - (settings['years'] * regions_for_selection['rsl_gia_dev'][region]),color='C0')
            plt.grid()
            plt.draw()
            plt.pause(0.001)
            print('TG id:  ' + str(regions_for_selection['name'][region]))
            print('TG id:  ' + str(regions_for_selection['id'][region]))
            print('TG time: ' + str(settings['years'][np.where(np.isfinite(regions_for_selection['height_corr'][region, :]))[0][0]]) + ' - ' + str(settings['years'][np.where(np.isfinite(regions_for_selection['height_corr'][region, :]))[0][-1]]))
            trend_basin = gentools.lsqtrend(settings['years'], regions_for_selection['height_corr'][region, :] - regions_for_selection['rsl_grd_dev'][region] - (settings['years'] * regions_for_selection['rsl_gia_dev'][region]))
            print('Basin trend: ' + str(trend_basin)[:4])
            print('-----------------')
            if len(vlm_name) == 0:
                print('No VLM estimate available')
            else:
                for vlm in range(len(vlm_name)):
                    print(vlm_name[vlm].rjust(5, ' ') + ' ' + str(vlm_est[vlm])[:5] + ' ' + str(vlm_err[vlm])[:5] + '   ' + str(vlm_est[vlm] + trend_basin)[:5])
            print('-----------------')
            print("    region_list.append({'id': " + str(regions_for_selection['id'][region]) + ", 'vlm_id':" + str(vlm_name.tolist()) + "})")
            print('-----------------')
            input("Press Enter to continue...")
            plt.close()
    return

def create_NWPA_list():
    region_list = []
    region_list.append({'id': [130], 'vlm_id': []})  # JAPAN
    region_list.append({'id': [132], 'vlm_id': ['G111', 'J053', 'J574', 'J971', 'J972']})  # JAPAN
    region_list.append({'id': [133], 'vlm_id': ['J094', 'ALTTG']})  # JAPAN
    region_list.append({'id': [134], 'vlm_id':['G012', 'G155', 'G208', 'J070', 'J377', 'J653', 'J826', 'SMST', 'Z113', 'Z203']})  # JAPAN
    region_list.append({'id': [152, 1356], 'vlm_id': []})  # TAIWAN
    region_list.append({'id': [153, 545], 'vlm_id':[]}) # KEELUNG
    region_list.append({'id': [155, 823], 'vlm_id':['HNLC', 'ZHN1', 'ALTTG']})
    region_list.append({'id': [300], 'vlm_id': ['HILO', 'ALTTG']})
    region_list.append({'id': [513], 'vlm_id': []})
    region_list.append({'id': [521], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [523], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [528], 'vlm_id': []})
    region_list.append({'id': [540, 2130], 'vlm_id':['GUAM', 'GUUG']})
    region_list.append({'id': [595], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [598], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [614], 'vlm_id': []})
    region_list.append({'id': [723, 1513], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [756], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [934], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [955], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [956], 'vlm_id': ['SKM1', 'ALTTG']})
    region_list.append({'id': [979], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [997], 'vlm_id': []})
    region_list.append({'id': [1007], 'vlm_id': []})
    region_list.append({'id': [1066,1627], 'vlm_id':['JEJU', 'ALTTG']})
    region_list.append({'id': [1108], 'vlm_id': []})
    region_list.append({'id': [1155], 'vlm_id':['ALTTG']})
    region_list.append({'id': [1164], 'vlm_id': []})
    region_list.append({'id': [1217, 1838], 'vlm_id':['MAJB', 'ALTTG']})
    region_list.append({'id': [1365], 'vlm_id': ['SKCH', 'ALTTG']})
    region_list.append({'id': [1370, 1925], 'vlm_id':['POHN']})
    region_list.append({'id': [1372], 'vlm_id': []})
    region_list.append({'id': [1445], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1446], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1473], 'vlm_id': []})
    region_list.append({'id': [1474], 'vlm_id':['CNMR']}) # Spikes
    region_list.append({'id': [1489], 'vlm_id': []})
    region_list.append({'id': [1490], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1527, 2042], 'vlm_id':['ALTTG']})
    region_list.append({'id': [1546], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1568], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1588], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1628], 'vlm_id':['ALTTG']})
    region_list.append({'id': [1675], 'vlm_id': []})
    region_list.append({'id': [1699], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1739, 1804], 'vlm_id': ['KIRI']})
    region_list.append({'id': [1844], 'vlm_id':['NAUR']})
    region_list.append({'id': [1860], 'vlm_id':['PNGM']})
    region_list.append({'id': [1861], 'vlm_id':['SOLO']})
    region_list.append({'id': [2128], 'vlm_id':['NIHO', 'UPO5', 'UPO6']})
    return (region_list)

def create_SATL_list():
    region_list = []
    region_list.append({'id': [157, 832], 'vlm_id': ['IGM1', 'BUE1']})
    region_list.append({'id': [223], 'vlm_id': []})
    region_list.append({'id': [284], 'vlm_id': ['DRBA']})
    region_list.append({'id': [431], 'vlm_id': ['UYMO']})
    region_list.append({'id': [433], 'vlm_id':[]})
    region_list.append({'id': [764], 'vlm_id':['UYLP']})
    region_list.append({'id': [820], 'vlm_id': ['PELB', 'PLBA']})
    region_list.append({'id': [826], 'vlm_id': ['SIMO']})
    region_list.append({'id': [910], 'vlm_id':[]})
    region_list.append({'id': [911], 'vlm_id': []})
    region_list.append({'id': [914], 'vlm_id': ['FG08']})
    region_list.append({'id': [950], 'vlm_id':['PLET']})
    region_list.append({'id': [1443], 'vlm_id':['RBAY']})
    return (region_list)

def create_EAPA_list():
    region_list = []
    region_list.append({'id': [10, 437], 'vlm_id':['EBMD','TIBB','PBL1','P224','UCSF','ALTTG']})
    region_list.append({'id': [127], 'vlm_id':['SEAI', 'SEAT', 'SMAI', 'ALTTG']})
    region_list.append({'id': [158, 256], 'vlm_id':['NSSS','SIO3','ALTTG']})
    region_list.append({'id': [163, 581], 'vlm_id':['ACP1','IGN1']})
    region_list.append({'id': [165], 'vlm_id': []}) ##### NEW Tofino
    region_list.append({'id': [166], 'vlm_id':['ALBH', 'BCNS', 'FRID', 'PGC5', 'ALTTG']})
    region_list.append({'id': [167], 'vlm_id': ['BCPR']})
    region_list.append({'id': [175, 193, 1245, 1255], 'vlm_id': ['BCSF', 'BRNB', 'NWE2', 'ALTTG']})
    region_list.append({'id': [225], 'vlm_id': []})
    region_list.append({'id': [245, 717], 'vlm_id':['LBC2','LBCH','ALTTG']})
    region_list.append({'id': [265], 'vlm_id': ['FTS1', 'FTS5', 'FTS6']})
    region_list.append({'id': [266], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [377], 'vlm_id':['BRAN', 'CAL5', 'CASM', 'CBHS', 'CRHS', 'CSDH', 'CSN1', 'DSHS', 'ECCO', 'ELSC', 'FXHS', 'LAPC', 'LASC', 'LEEP', 'LFRS', 'MHMS', 'MTA1', 'NOPK', 'P799', 'P800', 'PKRD', 'PVHS', 'SILK', 'SPK1', 'TORP', 'UCLP', 'USC1', 'VIMT', 'WMAP', 'WRHS', 'ALTTG']})
    region_list.append({'id': [378], 'vlm_id':['CACC','PTSG','ALTTG']})
    region_list.append({'id': [384], 'vlm_id':['FRID','P439', 'SC02']})
    region_list.append({'id': [385], 'vlm_id': ['MKAH', 'NEAH','ALTTG']})
    region_list.append({'id': [405], 'vlm_id':['AB50', 'JNU1', 'ALTTG']})
    region_list.append({'id': [426], 'vlm_id': ['BIS1', 'BIS5', 'BIS6']})
    region_list.append({'id': [458], 'vlm_id': []}) #### NEW MATARANI
    region_list.append({'id': [464], 'vlm_id': []}) #### NEW PUNTARENAS
    region_list.append({'id': [475], 'vlm_id': []}) #### NEW TALARA
    region_list.append({'id': [495], 'vlm_id': ['AB44', 'ALTTG']})
    region_list.append({'id': [508], 'vlm_id':['CHOR', 'DAPK', 'DCAN', 'HIMT','P523', 'P525', 'P528', 'USLO', 'ALTTG']})
    region_list.append({'id': [510, 511], 'vlm_id': ['BN02', 'UCNF']})
    region_list.append({'id': [527], 'vlm_id': ['PTAL', 'QUL2']})
    region_list.append({'id': [554], 'vlm_id': ['BCOV']})
    region_list.append({'id': [566], 'vlm_id': ['EYAC', 'ALTTG']})
    region_list.append({'id': [567], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [619], 'vlm_id': []}) #### NEW CALDERA
    region_list.append({'id': [688, 1152], 'vlm_id':['NCOW', 'PGC5', 'SC04']})
    region_list.append({'id': [689], 'vlm_id': []})
    region_list.append({'id': [695], 'vlm_id': ['IPAZ', 'LPAZ']}) #### NEW LA PAZ
    region_list.append({'id': [766], 'vlm_id':['CALK', 'CCCS','OEOC', 'SBCC', 'THMS', 'TRAK', 'WHYT']})
    region_list.append({'id': [795], 'vlm_id': ['CIC1', 'IMIE']}) #### NEW ENSENADA
    region_list.append({'id': [829], 'vlm_id': []}) #### NEW QUEEN CHARLOTTE CITY
    region_list.append({'id': [984], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1011], 'vlm_id': []}) #### NEW ACAJUTLA
    region_list.append({'id': [1013], 'vlm_id':['CSST', 'HVYS', 'OVLS', 'P548', 'P549']})
    region_list.append({'id': [1067], 'vlm_id':['AC44', 'ANC2', 'TBON', 'TSEA', 'ZAN1', 'ALTTG']})
    region_list.append({'id': [1071], 'vlm_id': ['BCPH', 'ALTTG']})
    region_list.append({'id': [1196], 'vlm_id':['NEWP', 'ONAB', 'P367', 'ALTTG']})
    region_list.append({'id': [1242], 'vlm_id': ['BAMF', 'ALTTG']})
    region_list.append({'id': [1269], 'vlm_id':['P364', 'P365', 'ALTTG']})
    region_list.append({'id': [1274], 'vlm_id': ['CAL3']})
    region_list.append({'id': [1323], 'vlm_id': ['QUAD']})
    region_list.append({'id': [1325], 'vlm_id':['CHCM', 'COUP', 'P436', 'P437', 'SQIM', 'WHD1', 'WHD5', 'WHD6', 'ALTTG']})
    region_list.append({'id': [1329], 'vlm_id':['ALTTG']})
    region_list.append({'id': [1352], 'vlm_id':['CAMT', 'P171', 'P210', 'P231']})
    region_list.append({'id': [1354], 'vlm_id': ['P398', 'P415', 'RYMD', 'ALTTG']})
    region_list.append({'id': [1371], 'vlm_id': []})
    region_list.append({'id': [1394], 'vlm_id':['P193', 'P194', 'PTRB', 'ALTTG']})
    region_list.append({'id': [1450], 'vlm_id': []})
    region_list.append({'id': [1633], 'vlm_id': ['BELI', 'P439', 'ALTTG']})
    region_list.append({'id': [1634], 'vlm_id': ['AB07', 'ALTTG']})
    region_list.append({'id': [1639], 'vlm_id':['P058', 'P159', 'P160', 'P161', 'P162', 'P169', 'ALTTG']}) #### NEW N. SPIT, HUMBOLDT BAY
    region_list.append({'id': [1640], 'vlm_id': ['CABL', 'ALTTG']})
    region_list.append({'id': [1645], 'vlm_id': []})
    region_list.append({'id': [1799], 'vlm_id': ['HOLB', 'ALTTG']})
    region_list.append({'id': [2125], 'vlm_id': ['P059', 'P184', 'ALTTG']})
    region_list.append({'id': [2126], 'vlm_id':['CASN', 'COPR', 'CSST', 'P519', 'P520', 'P548', 'RCA2', 'UCSB', 'ALTTG']})
    region_list.append({'id': [2127], 'vlm_id':['ALBH', 'P064', 'P435', 'P436', 'PTAA', 'SQIM', 'ALTTG']})
    region_list.append({'id': [2214], 'vlm_id': ['CHZZ', 'P396', 'P405', 'TILL']})
    region_list.append({'id': [2298], 'vlm_id': ['P401', 'P402', 'P815', 'P816']}) # LA PUSH
    region_list.append({'id': [2301], 'vlm_id': ['AC02']}) # ALITAK
    region_list.append({'id': [2302], 'vlm_id': ['AC25', 'BAY1', 'BAY2', 'BAY5', 'BAY6', 'CDB8']}) ### KING COVE, DEER PASSAGE
    region_list.append({'id': [2329], 'vlm_id':['CAHB', 'CAP1', 'CHAB', 'HSIB', 'JRSC', 'MSHP', 'OXMT', 'P176', 'P177', 'P178', 'P219', 'P220', 'P221', 'P222', 'P223', 'P225', 'SCCP', 'SLAC', 'STFU', 'SWEP', 'WIN2', 'WINT', 'ZOA1', 'ZOA2', 'ALTTG']})
    region_list.append({'id': [2330], 'vlm_id':['BRIB', 'CAFA', 'DIAB', 'LRA3', 'OHLN', 'P224', 'P248', 'P261', 'P262', 'P266', 'PTRO', 'ROCP', 'SRB1']})
    return (region_list)

def create_STNA_list():
    region_list = []
    region_list.append({'id': [112], 'vlm_id': []})
    region_list.append({'id': [161, 828], 'vlm_id':['TXGA', 'TXLM', 'UHC1', 'UHC2', 'UHC3','ALTTG']})
    region_list.append({'id': [169], 'vlm_id': ['ACP1', 'ALTTG']})
    region_list.append({'id': [188], 'vlm_id':['FLKW', 'KWST', 'KYW1', 'KYW2', 'KYW5', 'KYW6', 'ALTTG']})
    region_list.append({'id': [234], 'vlm_id':['SCCC', 'SCHA']})
    region_list.append({'id': [246], 'vlm_id': ['PCLA', 'ALTTG']})
    region_list.append({'id': [270], 'vlm_id': ['ORMD']})
    region_list.append({'id': [316, 716], 'vlm_id': []})
    region_list.append({'id': [363], 'vlm_id': ['MIA3', 'ZMA1']})
    region_list.append({'id': [368], 'vlm_id':['BMPD', 'BRMU']})
    region_list.append({'id': [395], 'vlm_id':['SAVA', 'ALTTG']})
    region_list.append({'id': [396], 'vlm_id':['NCCH', 'NCFF']})
    region_list.append({'id': [418], 'vlm_id': []})
    region_list.append({'id': [428], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [497, 919], 'vlm_id':['AMBI', 'TXLN', 'ALTTG']})
    region_list.append({'id': [520], 'vlm_id': ['MCD5', 'MCD6', 'ALTTG']})
    region_list.append({'id': [526], 'vlm_id': ['GRIS', 'ALTTG']})
    region_list.append({'id': [538], 'vlm_id':['TXPO', 'ALTTG']})
    region_list.append({'id': [563], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [759], 'vlm_id':['PRLT', 'PRMI', 'ALTTG']})
    region_list.append({'id': [768], 'vlm_id': []})
    region_list.append({'id': [1001], 'vlm_id':['BYSP', 'ZSU4', 'ALTTG']})
    region_list.append({'id': [1038], 'vlm_id': []})
    region_list.append({'id': [1106], 'vlm_id':['FMYR', 'ALTTG']})
    region_list.append({'id': [1107], 'vlm_id': ['NAPL', 'ALTTG']})
    region_list.append({'id': [1156], 'vlm_id':['ALDI', 'MOB1', 'MOB2', 'MOB5', 'MOB6', 'ALTTG']})
    region_list.append({'id': [1193], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1393], 'vlm_id':['STVI', 'VITH', 'ALTTG']})
    region_list.append({'id': [1444], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1447, 2118], 'vlm_id': ['CRO1', 'VIKH']})
    region_list.append({'id': [1638], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1641], 'vlm_id': ['PNCY', 'ALTTG']})
    region_list.append({'id': [1701], 'vlm_id':['ALTTG']})
    region_list.append({'id': [1835], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1858], 'vlm_id': ['MIA3', 'RMND', 'ZMA1', 'ALTTG']})
    region_list.append({'id': [1903], 'vlm_id': ['TXCC', 'ALTTG']})
    region_list.append({'id': [1910], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [2021], 'vlm_id': []})
    region_list.append({'id': [2123], 'vlm_id': ['CCV3', 'CCV5', 'CCV6', 'ALTTG']})
    region_list.append({'id': [2294], 'vlm_id':['NCBX']})
    region_list.append({'id': [2295], 'vlm_id': ['NCBE', 'ALTTG']})
    return (region_list)

def create_IOSP_list():
    region_list = []
    region_list.append({'id': [43], 'vlm_id': ['ALTTG']})  # MUMBAY
    region_list.append({'id': [44], 'vlm_id': []})  # ADEN
    region_list.append({'id': [65, 196, 549], 'vlm_id': ['BANK', 'FTDN', 'SYDN', 'UNSW','ALTTG']})  # SYDNEY
    region_list.append({'id': [93], 'vlm_id': ['MOBS', 'SG36']})  # ADD WILLIAMSTOWN
    region_list.append({'id': [111], 'vlm_id': []})  # FREMANTLE
    region_list.append({'id': [136, 252, 1643], 'vlm_id':['DUND', 'DUNT', 'OUS2', 'OUSD', 'ALTTG']})  # DUNEDIN
    region_list.append({'id': [150, 217], 'vlm_id': ['AUCK', 'AUKT','ALTTG']})  # AUCKLAND
    region_list.append({'id': [189], 'vlm_id':['PTHL']}) #### NEW PORT HEDLAND
    region_list.append({'id': [204], 'vlm_id': ['ALTTG']})  # KARACHI
    region_list.append({'id': [205], 'vlm_id': []})  # CHENNAI
    region_list.append({'id': [216], 'vlm_id':['CBRK']})  # PORT PIRIE
    region_list.append({'id': [221, 500], 'vlm_id':['AVLN', 'WEL1', 'WGTN', 'WGTT']}) # WELLINGTON
    region_list.append({'id': [230], 'vlm_id':[]}) #### NEW PORT LINCOLN
    region_list.append({'id': [247, 259], 'vlm_id': ['LYTT', 'MQZG', 'ALTTG']})  # PORT LYTTELTON
    region_list.append({'id': [248, 1183, 1248, 1275, 1351, 1746, 1894, 1895, 1896], 'vlm_id': ['NTUS','ALTTG']})  # SINGAPORE
    region_list.append({'id': [260], 'vlm_id': []})  # JOLO, SULU
    region_list.append({'id': [267, 320, 837, 1335], 'vlm_id': ['NEWE']})  # NEWCASTLE
    region_list.append({'id': [310], 'vlm_id':['ALTTG']})  # YAMBA
    region_list.append({'id': [333, 1034, 1674, 1698, 1891], 'vlm_id':['HKFN', 'HKKT', 'HKLT', 'HKMW', 'HKOH', 'HKPC', 'HKSC', 'HKSL', 'HKSS', 'HKST', 'HKWS', 'T430']}) ### HONGKONG
    region_list.append({'id': [386], 'vlm_id':['ALTTG']})  # THEVENARD
    region_list.append({'id': [394], 'vlm_id': []})  # CEBU
    region_list.append({'id': [438], 'vlm_id':[]})   ####  COCHIN
    region_list.append({'id': [448], 'vlm_id':['ADE1','ADE2','PTSV','S021', 'ALTTG']})  # PORT ADELAIDE
    region_list.append({'id': [449], 'vlm_id':[]}) #### NEW KO SICHANG
    region_list.append({'id': [539], 'vlm_id': ['ASPA']}) #####  PAGO PAGO
    region_list.append({'id': [564], 'vlm_id': ['SLAD', 'ALTTG']})  # MACKAY
    region_list.append({'id': [637], 'vlm_id': ['TOW2', 'ALTTG']})  # TOWNSVILLE
    region_list.append({'id': [683], 'vlm_id':['BUR1', 'BUR2', 'RHPT', 'ALTTG']})  # BURNIE
    region_list.append({'id': [822], 'vlm_id':['CLEV', 'WOOL']})  # BRISBANE
    region_list.append({'id': [825], 'vlm_id': ['ALTTG']})  # GLADSTONE
    region_list.append({'id': [831], 'vlm_id':['PTKL', 'ALTTG']})  # PORT KEMBLA
    region_list.append({'id': [834], 'vlm_id': ['ALTTG']})  # BUNBURY (CHECKEN)
    region_list.append({'id': [835], 'vlm_id': []}) #### NEW KURUMBA
    region_list.append({'id': [838], 'vlm_id':['HOB2', 'ALTTG']})  # HOBART
    region_list.append({'id': [935], 'vlm_id': ['00NA', '01NA', 'DARM', 'LKYA']}) #### NEW DARWIN
    region_list.append({'id': [953], 'vlm_id': []})  # CAIRNS
    region_list.append({'id': [957], 'vlm_id': ['ALBY', 'ALTTG']})  # ALBANY
    region_list.append({'id': [1033], 'vlm_id':['MRNT', 'STNY', 'ALTTG']})  # STONY POINT
    region_list.append({'id': [1065], 'vlm_id': ['WHNG']})  # WHANGAREI HARBOUR (MARSDEN POINT)
    region_list.append({'id': [1069], 'vlm_id':['MCKN', 'ALTTG']})  # VICTOR HARBOUR
    region_list.append({'id': [1114], 'vlm_id': ['ESPA', 'ALTTG']})  #### NEW ESPERANCE
    region_list.append({'id': [1154], 'vlm_id': ['BNDY', 'ALTTG']})  # BUNDABERG, BURNETT HEADS
    region_list.append({'id': [1159], 'vlm_id': ['BRO1']}) #### NEW BROOME
    region_list.append({'id': [1160], 'vlm_id': []}) ##### NEW MILNER BAY
    region_list.append({'id': [1216], 'vlm_id': ['SPBY', 'ALTTG']})  # SPRING BAY
    region_list.append({'id': [1246], 'vlm_id':['SLAD']})  # HAY POINT
    region_list.append({'id': [1249], 'vlm_id': []}) #### NEW MORMUGAO
    region_list.append({'id': [1252], 'vlm_id': ['PALA']})  # MALAKAL-B
    region_list.append({'id': [1270], 'vlm_id': ['ALTTG']}) ##### NEW HALDIA
    region_list.append({'id': [1369], 'vlm_id': []})  # GANGRA
    region_list.append({'id': [1395], 'vlm_id': []})  # OKHA
    region_list.append({'id': [1420], 'vlm_id': ['ALTTG']})  # WALLAROO II
    region_list.append({'id': [1423], 'vlm_id': []})  # MANGALORE
    region_list.append({'id': [1471], 'vlm_id': []})  # PORT DOUGLAS 2
    region_list.append({'id': [1475], 'vlm_id': ['ALTTG']}) #### NEW DANANG
    region_list.append({'id': [1492], 'vlm_id': ['TOW2', 'ALTTG']})  # CAPE FERGUSON
    region_list.append({'id': [1493], 'vlm_id': ['ALTTG']})  # MOOLOOLABA
    region_list.append({'id': [1501], 'vlm_id': ['SLEU']}) #### NEW POINTE DES GALETS
    region_list.append({'id': [1547], 'vlm_id': ['PTLD', 'ALTTG']})  # PORTLAND
    region_list.append({'id': [1550], 'vlm_id': ['ALTTG']})  # PORT GILES
    region_list.append({'id': [1569], 'vlm_id': ['ALTTG']})  # SHUTE HARBOUR
    region_list.append({'id': [1589], 'vlm_id': ['ALTTG']})  # TANJUNG GELANG
    region_list.append({'id': [1590], 'vlm_id':['TRNG', 'ALTTG']}) #### NEW TAURANGA
    region_list.append({'id': [1591], 'vlm_id': ['ALTTG']}) #### NEW PELABUHAN KELANG
    region_list.append({'id': [1592], 'vlm_id':['KUAL']})  # CENDERING
    region_list.append({'id': [1593], 'vlm_id': ['ALTTG']}) #### NEW TANJUNG KELING
    region_list.append({'id': [1594], 'vlm_id':['ALTTG']})  # LUMUT
    region_list.append({'id': [1595], 'vlm_id':['ALTTG']})  # PULAU PINANG
    region_list.append({'id': [1600], 'vlm_id': ['ALTTG']})  # ZANZIBAR
    region_list.append({'id': [1630], 'vlm_id': ['ALTTG']})  # LUCINDA
    region_list.append({'id': [1673], 'vlm_id': ['VACS']}) #### NEW PORT LOUIS II
    region_list.append({'id': [1676], 'vlm_id': ['ALTTG']}) #### NEW PULAU LANGKAWI
    region_list.append({'id': [1678], 'vlm_id': ['ALTTG']})  # PULAU TIOMAN
    region_list.append({'id': [1702], 'vlm_id': []})  # TANJUNG SEDILI
    region_list.append({'id': [1703], 'vlm_id': ['GETI', 'ALTTG']}) #### NEW GETING
    region_list.append({'id': [1707], 'vlm_id': []})  # GAN
    region_list.append({'id': [1733], 'vlm_id': []})  # KOTA KINABALU
    region_list.append({'id': [1734], 'vlm_id': []})  # TAWAU
    region_list.append({'id': [1753], 'vlm_id': ['ALTTG']})  # MALE-B, HULULE
    region_list.append({'id': [1760], 'vlm_id': ['ALTTG']})  # ROSSLYN BAY
    region_list.append({'id': [1761], 'vlm_id': ['HIL1','ALTTG']}) #### NEW HILLARYS
    region_list.append({'id': [1805], 'vlm_id': ['LAUT', 'ALTTG']}) #### NEW LAUTOKA
    region_list.append({'id': [1834], 'vlm_id': []}) #### NEW SANDAKAN
    region_list.append({'id': [1836], 'vlm_id': ['ALTTG']})  # LORNE
    region_list.append({'id': [1876], 'vlm_id': []}) #### NEW KUDAT
    region_list.append({'id': [1877], 'vlm_id': []}) #### NEW LAHAT DATU
    region_list.append({'id': [1879], 'vlm_id': ['ALTTG']}) #### NEW LABUAN 2
    region_list.append({'id': [2072], 'vlm_id': ['ALTTG']})  # PORT ALMA
    region_list.append({'id': [2073], 'vlm_id': ['ALTTG']}) #### NEW URANGAN II
    region_list.append({'id': [2074], 'vlm_id': ['ALTTG']})  # BOWEN II
    region_list.append({'id': [2310], 'vlm_id': ['ALTTG']}) #### NEW BRUNSWICK HEADS
    region_list.append({'id': [2312], 'vlm_id': ['ALTTG']}) #### NEW JERVIS BAY II
    region_list.append({'id': [2313], 'vlm_id':['ALTTG']})  # PORT MACQUARIE
    region_list.append({'id': [2316], 'vlm_id':['ALTTG']})  # PORT STEPHENS (TOMAREE)
    return (region_list)

def create_SPNA_list():
    region_list = []
    region_list.append({'id': [1, 1294], 'vlm_id': ['BRST', 'GUIP', 'LPPZ', 'PLAB', 'ALTTG']})
    region_list.append({'id': [3, 334], 'vlm_id':['MAID', 'SHEE', 'SHOE', 'ALTTG']})
    region_list.append({'id': [5], 'vlm_id': ['HOLY', 'LYN1']})
    region_list.append({'id': [7], 'vlm_id': ['CADE', 'TGBU', 'TGCU', 'ALTTG']})
    region_list.append({'id': [8, 54], 'vlm_id': ['D777', 'D779']})
    region_list.append({'id': [9, 19, 22], 'vlm_id': ['DELF', 'DLF1', 'DLF5', 'DLFT', 'OVRN', 'RDAM']})
    region_list.append({'id': [11], 'vlm_id': ['D773', 'WARN']})
    region_list.append({'id': [12, 1637], 'vlm_id':['NJHT', 'NJI2', 'NYBK', 'NYBP', 'NYBR', 'NYJM', 'NYOB', 'NYPR', 'NYQN', 'SA22', 'ALTTG']})
    region_list.append({'id': [13], 'vlm_id': []})
    region_list.append({'id': [15, 765], 'vlm_id': ['DARE', 'LIVE']})
    region_list.append({'id': [20], 'vlm_id':['EEKL', 'OOS1', 'SASG', 'VLIS', 'ZEEB', 'ALTTG']})
    region_list.append({'id': [21, 361], 'vlm_id': ['ABER', 'KINT']})
    region_list.append({'id': [23], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [24], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [25], 'vlm_id': ['MAKK', 'ALTTG']})
    region_list.append({'id': [32], 'vlm_id': ['ALK2', 'HEER', 'IJMU', 'WIJK', 'ALTTG']})
    region_list.append({'id': [33], 'vlm_id': ['OSLS', 'ALTTG']})
    region_list.append({'id': [34, 1748], 'vlm_id': ['TRDS', 'ALTTG']})
    region_list.append({'id': [35, 36, 1551], 'vlm_id': ['SDYK', 'SVNS', 'ALTTG']})
    region_list.append({'id': [39, 87, 168, 2100], 'vlm_id': ['CAVA', 'CGIA', 'SDNA', 'SFEL', 'VEN1', 'VEN2', 'ALTTG']})
    region_list.append({'id': [45], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [47], 'vlm_id': ['STAS', 'ALTTG']})
    region_list.append({'id': [52], 'vlm_id': ['CASC', 'FCUL', 'IGP0']})
    region_list.append({'id': [55], 'vlm_id': ['WES1']})
    region_list.append({'id': [56], 'vlm_id': ['D799', 'SASS']})
    region_list.append({'id': [58], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [59, 2090], 'vlm_id': ['GENO', 'GENU']})
    region_list.append({'id': [61], 'vlm_id': ['AXPV', 'MARS', 'PRIE']})
    region_list.append({'id': [62], 'vlm_id': ['OSLS']})
    region_list.append({'id': [76], 'vlm_id': []})
    region_list.append({'id': [80], 'vlm_id': ['ESBC', 'ESBH']})
    region_list.append({'id': [81], 'vlm_id': ['SMID']})
    region_list.append({'id': [82], 'vlm_id':['0LOD', '1MAL', 'BUDP', 'ALTTG']})
    region_list.append({'id': [89], 'vlm_id': ['HIRS']})
    region_list.append({'id': [91], 'vlm_id': []})
    region_list.append({'id': [95], 'vlm_id':['MOR7', 'MORP', 'NCAS', 'NSLG', 'NSTG', 'ALTTG']})
    region_list.append({'id': [96, 2352], 'vlm_id':['DART', 'DRT2', 'HFXA', 'HLFX', 'TANT', 'ALTTG']})
    region_list.append({'id': [98], 'vlm_id': []})
    region_list.append({'id': [104], 'vlm_id':['CA05', 'CAG1', 'CAGL', 'CAGZ', 'UCAG']})
    region_list.append({'id': [113], 'vlm_id': []})
    region_list.append({'id': [119], 'vlm_id': []})
    region_list.append({'id': [120], 'vlm_id': ['GESR']})
    region_list.append({'id': [135], 'vlm_id':['NJGC', 'NJTW', 'PAPH', 'ALTTG']})
    region_list.append({'id': [148], 'vlm_id':['ANP5', 'ANP6', 'BACO', 'LOYK', 'SA15', 'UMBC', 'ALTTG']})
    region_list.append({'id': [154, 1009, 1817, 2099], 'vlm_id': ['KOPE', 'TRIE']})
    region_list.append({'id': [162], 'vlm_id': ['LAGO']})
    region_list.append({'id': [180], 'vlm_id':['NJGT', 'ALTTG']})
    region_list.append({'id': [183], 'vlm_id':['MEGO', 'MEYA', 'ALTTG']})
    region_list.append({'id': [190], 'vlm_id':['SABB', 'SABS']})
    region_list.append({'id': [195], 'vlm_id':['SJNB', 'SJPA']})
    region_list.append({'id': [202], 'vlm_id':['CAM1', 'LIZ1', 'NEWL', 'ALTTG']})
    region_list.append({'id': [208, 960], 'vlm_id': ['ALAC']})
    region_list.append({'id': [219, 861], 'vlm_id': ['BELF']})
    region_list.append({'id': [224], 'vlm_id':['DED2', 'DEMI', 'ALTTG']})
    region_list.append({'id': [235], 'vlm_id':['MAMI', 'ALTTG']})
    region_list.append({'id': [236], 'vlm_id': ['TERS', 'ALTTG']})
    region_list.append({'id': [255, 2277], 'vlm_id':['DVTG']})
    region_list.append({'id': [257, 1832], 'vlm_id':['CARF']})
    region_list.append({'id': [258], 'vlm_id':['PDEL', 'VFDC', 'ALTTG']})
    region_list.append({'id': [263], 'vlm_id':['CHIO', 'OSHQ', 'PMTG', 'SOTN']})
    region_list.append({'id': [286], 'vlm_id': ['EASN', 'SWAN']})
    region_list.append({'id': [288], 'vlm_id': ['MESB', 'NHUN']})
    region_list.append({'id': [299, 399, 1635], 'vlm_id':['DRV1', 'DRV5', 'DRV6', 'LOY2', 'LOYZ', 'ALTTG']})
    region_list.append({'id': [302], 'vlm_id':['TGDE', 'ALTTG']})
    region_list.append({'id': [303], 'vlm_id':['GRAF', 'LLAG', 'TN01']})
    region_list.append({'id': [311], 'vlm_id':['ANP5', 'ANP6', 'LOYF', 'USNA', 'ALTTG']})
    region_list.append({'id': [312], 'vlm_id': []})
    region_list.append({'id': [313], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [314], 'vlm_id':['STOR', 'SWTG', 'ALTTG']})
    region_list.append({'id': [332], 'vlm_id': ['EPRT', 'ALTTG']})
    region_list.append({'id': [335, 336], 'vlm_id':['BAR3', 'MAID', 'SHEE', 'STR3', 'TBSB', 'UEL0']})
    region_list.append({'id': [350], 'vlm_id': ['OSHQ', 'PMTG', 'SANO', 'SOTN', 'ALTTG']})
    region_list.append({'id': [351], 'vlm_id':['MAD4', 'NPRI', 'URIL', 'ALTTG']})
    region_list.append({'id': [352, 685], 'vlm_id': ['SPLT']})
    region_list.append({'id': [353], 'vlm_id': []})
    region_list.append({'id': [360], 'vlm_id':['GODE', 'GODN', 'GODZ', 'LOYB', 'NRL1', 'S071', 'SA02', 'USN3', 'USNO', 'ALTTG']})
    region_list.append({'id': [362, 856], 'vlm_id':['LAMT', 'NJHT', 'NYBK', 'NYBP', 'NYBR', 'NYJM', 'NYQN']})
    region_list.append({'id': [366], 'vlm_id': ['SHK1', 'SHK5', 'SHK6', 'ALTTG']})
    region_list.append({'id': [367], 'vlm_id': ['AMTS', 'ALTTG']})
    region_list.append({'id': [392], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [393], 'vlm_id':['STJ2', 'STJ5', 'STJO']})
    region_list.append({'id': [397], 'vlm_id': ['D799', 'SASS']})
    region_list.append({'id': [398], 'vlm_id': []})
    region_list.append({'id': [401], 'vlm_id': []})
    region_list.append({'id': [412], 'vlm_id':['MDLT', 'SOL1', 'ALTTG']})
    region_list.append({'id': [413, 489], 'vlm_id': ['OOST', 'ZEEB', 'ALTTG']})
    region_list.append({'id': [425], 'vlm_id':['AND1', 'ANDO']})
    region_list.append({'id': [427], 'vlm_id':['HUN3', 'PEGS', 'PETI', 'ALTTG']})
    region_list.append({'id': [429], 'vlm_id': ['CTGR', 'ALTTG']})
    region_list.append({'id': [430], 'vlm_id':['RIRU', 'ALTTG']})
    region_list.append({'id': [432], 'vlm_id': ['SWR1', 'SWRD', 'TLLG']})
    region_list.append({'id': [435], 'vlm_id': ['IOMS']})
    region_list.append({'id': [453], 'vlm_id':['CT3F', 'CT3I', 'HAVR', 'TANC', 'ALTTG']})
    region_list.append({'id': [454], 'vlm_id': ['DIPL', 'GOUE', 'SMTG', 'ALTTG']})
    region_list.append({'id': [455], 'vlm_id': ['BLER']})
    region_list.append({'id': [457, 1078], 'vlm_id':['CROI', 'STNA']})
    region_list.append({'id': [462], 'vlm_id': ['LOY3', 'RIC1']})
    region_list.append({'id': [466], 'vlm_id':['AUNI', 'ILDX', 'LROC', 'POR8', 'ALTTG']})
    region_list.append({'id': [468], 'vlm_id': ['COUD', 'DGLG', 'ALTTG']})
    region_list.append({'id': [470], 'vlm_id':['BRGG', 'OOS1', 'OOST', 'VLIS', 'ZEEB']})
    region_list.append({'id': [474], 'vlm_id': ['AMBL', 'SMEC', 'ALTTG']})
    region_list.append({'id': [483, 1898], 'vlm_id': ['BAIO', 'VIGO']})
    region_list.append({'id': [484, 763, 1808], 'vlm_id': ['ACOR']})
    region_list.append({'id': [485, 1051, 1807], 'vlm_id':['CANT', 'LARE', 'TRLV', 'TRVG', 'ALTTG']})
    region_list.append({'id': [486], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [488, 490, 2117], 'vlm_id':['ALGC', 'CEU1', 'ALTTG']})
    region_list.append({'id': [496, 1810], 'vlm_id': ['MALA', 'MLGA', 'ALTTG']})
    region_list.append({'id': [498], 'vlm_id':['ALGC', 'CEU1', 'ALTTG']})
    region_list.append({'id': [509], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [519], 'vlm_id':['ALTTG']})
    region_list.append({'id': [524], 'vlm_id': ['VARD', 'VARS']})
    region_list.append({'id': [525], 'vlm_id': ['BARH', 'ALTTG']})
    region_list.append({'id': [531], 'vlm_id': []})
    region_list.append({'id': [562], 'vlm_id': []})
    region_list.append({'id': [565, 1802], 'vlm_id':['AGUI', 'PLUZ', 'TERR', 'ALTTG']})
    region_list.append({'id': [568, 2064], 'vlm_id': ['LPAL', 'MAZO', 'ALTTG']})
    region_list.append({'id': [593], 'vlm_id': ['HRIA', 'TIAS', 'YAIZ', 'ALTTG']})
    region_list.append({'id': [597], 'vlm_id': ['GLPT', 'LOYX', 'VAGP']})
    region_list.append({'id': [593], 'vlm_id': ['HRIA', 'TIAS', 'YAIZ', 'ALTTG']})
    region_list.append({'id': [597], 'vlm_id':['GLPT', 'LOYX', 'VAGP']})
    region_list.append({'id': [636], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [638], 'vlm_id': ['REYK', 'ALTTG']})
    region_list.append({'id': [680], 'vlm_id': ['TRO1', 'TROM']})
    region_list.append({'id': [681], 'vlm_id': []})
    region_list.append({'id': [682], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [703], 'vlm_id': []})
    region_list.append({'id': [742], 'vlm_id':['ARDL', 'COLC', 'COLH']})
    region_list.append({'id': [754], 'vlm_id': ['GORS', 'LOWE', 'ALTTG']})
    region_list.append({'id': [755], 'vlm_id':['ALTTG']})
    region_list.append({'id': [758], 'vlm_id': []})
    region_list.append({'id': [760], 'vlm_id': ['DUB2', 'DUBR', 'ALTTG']})
    region_list.append({'id': [786], 'vlm_id': ['DENE', 'RED1', 'RED5', 'RED6', 'ALTTG']})
    region_list.append({'id': [788, 1468], 'vlm_id': ['EZEV', 'NICA', 'NICE', 'ALTTG']})
    region_list.append({'id': [789], 'vlm_id': ['HOL2']})
    region_list.append({'id': [791], 'vlm_id': ['GAIA']})
    region_list.append({'id': [802, 1074], 'vlm_id': ['EDIN']})
    region_list.append({'id': [830], 'vlm_id': ['LERI', 'LERW', 'LWTG']})
    region_list.append({'id': [839], 'vlm_id': ['ARGI']})
    region_list.append({'id': [848], 'vlm_id': ['NYCI', 'SG06', 'ZNY1']})
    region_list.append({'id': [916], 'vlm_id': []})
    region_list.append({'id': [936], 'vlm_id':['BLAP']})
    region_list.append({'id': [939, 1731], 'vlm_id': []})
    region_list.append({'id': [958], 'vlm_id':['AGDE', 'AGDS', 'CAPA', 'MTP2', 'MTPL', 'PZNA', 'SETE']})
    region_list.append({'id': [980], 'vlm_id':['PQRL', 'ALTTG']})
    region_list.append({'id': [981], 'vlm_id':['ALGC', 'CEU1']})
    region_list.append({'id': [982], 'vlm_id': ['PMTH']})
    region_list.append({'id': [1030, 2024], 'vlm_id': ['FUNC']})
    region_list.append({'id': [1036], 'vlm_id': ['HOE2', 'TGDA']})
    region_list.append({'id': [1037], 'vlm_id': ['BORJ', 'BORK', 'KRUM', 'TGBF', 'TGD2', 'ALTTG']})
    region_list.append({'id': [1068], 'vlm_id': ['CTDA', 'CTNH', 'SG06', 'ALTTG']})
    region_list.append({'id': [1075], 'vlm_id': []})
    region_list.append({'id': [1109], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1111], 'vlm_id': ['IMTS', 'ALTTG']})
    region_list.append({'id': [1112], 'vlm_id': []})
    region_list.append({'id': [1113], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1121], 'vlm_id': []})
    region_list.append({'id': [1153], 'vlm_id': ['NJCM', 'ALTTG']})
    region_list.append({'id': [1158], 'vlm_id': ['TSKT', 'ALTTG']})
    region_list.append({'id': [1197], 'vlm_id': ['FYHA']})
    region_list.append({'id': [1213], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1214], 'vlm_id': ['APPL', 'ALTTG']})
    region_list.append({'id': [1215], 'vlm_id': ['STRN']})
    region_list.append({'id': [1241], 'vlm_id': ['VIKC', 'ALTTG']})
    region_list.append({'id': [1247], 'vlm_id': ['GROI', 'PLOE']})
    region_list.append({'id': [1267], 'vlm_id': ['HONS']})
    region_list.append({'id': [1294], 'vlm_id': ['BRST', 'GUIP', 'LPPZ', 'ALTTG']})
    region_list.append({'id': [1295], 'vlm_id': ['HNPT', 'ALTTG']})
    region_list.append({'id': [1299], 'vlm_id': ['CBRM', 'SDNY', 'ALTTG']})
    region_list.append({'id': [1321], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1347], 'vlm_id': ['ROTG', 'ALTTG']})
    region_list.append({'id': [1349], 'vlm_id': ['ESCU', 'ALTTG']})
    region_list.append({'id': [1421], 'vlm_id': ['NYA1', 'NYAC', 'NYAL', 'ALTTG']})
    region_list.append({'id': [1448], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [1456], 'vlm_id': ['SCAC', 'ALTTG']})
    region_list.append({'id': [1469], 'vlm_id': ['BEAR', 'CREU', 'PERP']})
    region_list.append({'id': [1491], 'vlm_id': []})
    region_list.append({'id': [1505], 'vlm_id': ['LOFT']})
    region_list.append({'id': [1524], 'vlm_id': ['MEMA']})
    region_list.append({'id': [1526], 'vlm_id': ['EDIN']})
    region_list.append({'id': [1548], 'vlm_id': ['HERO', 'HERS', 'HERT', 'ALTTG']})
    region_list.append({'id': [1636], 'vlm_id': ['NCDU', 'ALTTG']})
    region_list.append({'id': [1747], 'vlm_id': ['SABL', 'SAGI', 'ALTTG']})
    region_list.append({'id': [1758], 'vlm_id': ['TAUT']})
    region_list.append({'id': [1759], 'vlm_id': ['0STR', 'ALTTG']})
    region_list.append({'id': [1764], 'vlm_id':['GIRO', 'ALTTG']})
    region_list.append({'id': [1771], 'vlm_id':['MACY']})
    region_list.append({'id': [1773], 'vlm_id': ['PBIL', 'ALTTG']})
    region_list.append({'id': [1774], 'vlm_id': ['DARE', 'LIVE']})
    region_list.append({'id': [1775], 'vlm_id':['KINL', 'ALTTG']})
    region_list.append({'id': [1794], 'vlm_id': ['BEE1', 'STBE']})
    region_list.append({'id': [1795], 'vlm_id': []})
    region_list.append({'id': [1797], 'vlm_id':['CSAR', 'ALTTG']})
    region_list.append({'id': [1803], 'vlm_id':['GRAF', 'LLAG', 'TN01', 'ALTTG']})
    region_list.append({'id': [1806], 'vlm_id':['GERN', 'SOPU', 'ALTTG']})
    region_list.append({'id': [1809], 'vlm_id':['LEBR', 'ALTTG']})
    region_list.append({'id': [1811], 'vlm_id':['BCLN', 'GAR1', 'PLAN']})
    region_list.append({'id': [1812], 'vlm_id':['TEJH']})
    region_list.append({'id': [1813], 'vlm_id': ['VALE', 'VCIA']})
    region_list.append({'id': [1854], 'vlm_id':['ASAP']})
    region_list.append({'id': [1855], 'vlm_id':['SCIL']})
    region_list.append({'id': [1859], 'vlm_id':['ZADA', 'ALTTG']})
    region_list.append({'id': [1867], 'vlm_id':['KLRE', 'ALTTG']})
    region_list.append({'id': [1871], 'vlm_id':['AVLS', 'ALTTG']})
    region_list.append({'id': [1915], 'vlm_id':['ARVE', 'ROYA']})
    region_list.append({'id': [2323], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [2324], 'vlm_id': ['ALTTG']})
    region_list.append({'id': [2325], 'vlm_id': ['NCPI', 'ALTTG']})
    return (region_list)
