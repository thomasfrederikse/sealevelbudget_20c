# Save Love numbers in numpy format
import numpy as np
import os

def main():
    print('Love numbers...')
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['fn_love_CM'] = settings['dir_data'] + 'Love_Wang/PREM-CM.dat'
    settings['fn_love_save'] = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    # Love numbers
    khl = np.loadtxt(settings['fn_love_CM'], skiprows=1)
    love = {}
    love['degree'] = khl[:360,0].astype(int)
    love['k'] = khl[:360+1,3]
    love['h'] = khl[:360+1,1]
    love['l'] = khl[:360+1,2]
    np.save(settings['fn_love_save'],love)
    return

if __name__ == '__main__':
    main()
