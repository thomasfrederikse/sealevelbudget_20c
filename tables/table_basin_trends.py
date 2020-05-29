# ---------------------------------------------------------------
# Write basin-mean trends to text string for easy import in LaTeX
# ---------------------------------------------------------------
import numpy as np
import os
import mod_gentools as gentools

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'

    print_table = []
    stats = np.load(settings['fn_stats'],allow_pickle=True).all()
    # Rij 1 GMSL
    region_names = ['Subpolar North Atlantic', 'Indian Ocean-South Pacific', 'Subtropical North Atlantic', 'East Pacific', 'South Atlantic', 'Northwest Pacific']
    for basin in range(6):
        print_table.append('\\textbf{'+region_names[basin]+ '}& \\multicolumn{2}{l}{1900-2018} & \\multicolumn{2}{l}{1957-2018} & \\multicolumn{2}{l}{1993-2018} \\\\')
        print_table.append('\\hline')
        print_table.append('Glaciers &' + print_sl_trend('grd_glac',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('grd_glac',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('grd_glac',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Greenland Ice Sheet &' + print_sl_trend('grd_GrIS',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('grd_GrIS',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('grd_GrIS',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Antarctic Ice Sheet &' + print_sl_trend('grd_AIS',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('grd_AIS',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('grd_AIS',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Terrestrial Water Storage &' + print_sl_trend('grd_tws',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('grd_tws',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('grd_tws',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Barystatic &' + print_sl_trend('grd_total',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('grd_total',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('grd_total',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Glacial Isostatic Adjustment &' + print_sl_trend('gia',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('gia',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('gia',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Steric & - & &'+print_sl_trend('steric',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('steric',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Sum of Contributors & - & &'+print_sl_trend('budget',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('budget',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Observed sea level &' + print_sl_trend('obs',basin,stats,'trend_1900_2018')+'&'+print_sl_trend('obs',basin,stats,'trend_1957_2018')+'&'+print_sl_trend('obs',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('Altimetry & - & & - & &' +print_sl_trend('altimetry',basin,stats,'trend_1993_2018')+'\\\\')
        print_table.append('& & & & & & \\\\')
    for i in print_table:
        print(i)


def print_sl_trend(varname,basin,stats,period):
    trend_mn = "{:.2f}".format(stats['basin'][basin][varname][period][0])
    trend_lo = "{:.2f}".format(stats['basin'][basin][varname][period][0]+stats['basin'][basin][varname][period][1])
    trend_hi = "{:.2f}".format(stats['basin'][basin][varname][period][0]+stats['basin'][basin][varname][period][2])
    print_str = trend_mn + ' &['+trend_lo+' '+trend_hi+']'
    return(print_str)
