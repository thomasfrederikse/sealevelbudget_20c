# -----------------------------------------------
# Plot 4 problematic stations for reponse letter
# - Fort Phrachula Chomklao (subsidence)
# - Port Taranaki (Variance)
# - Ilha Fiscal (Jump)
# - Nagaevo (Wobbles)
# -----------------------------------------------
import os
import numpy as np
import matplotlib.pyplot as plt
settings = {}
settings['dir_data'] = os.getenv('HOME') + '/Data/'
settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
settings['fn_regions_for_selection'] = settings['dir_budget'] + 'tg/regions_for_selection.npy'
settings['fn_alttg_for_region_selection'] = settings['dir_budget'] + 'vlm/alttg_for_region_selection.npy'
settings['fn_gps_data'] = settings['dir_budget'] + 'vlm/gps_data.npy'
settings['years'] = np.arange(1900, 2019)
regions_for_selection = np.load(settings['fn_regions_for_selection'], allow_pickle=True).all()
id_stat = [444,827,996,1032]
plt.figure(figsize=(8,6))
for idx,id in enumerate(id_stat):
    for stat in range(len(regions_for_selection['id'])):
        if id in regions_for_selection['id'][stat]:
            plt.subplot(2,2,idx+1)
            plt.plot(settings['years'],regions_for_selection['height_corr'][stat],color='black')
            plt.plot(settings['years'],regions_for_selection['height_corr'][stat],'x',color='black')

            plt.grid()
            plt.title(regions_for_selection['name'][stat],fontsize=8)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.ylabel('Height (mm)',fontsize=8)
plt.tight_layout()
plt.savefig('problematic_stations.png')






