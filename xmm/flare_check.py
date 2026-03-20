
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys

rate_limit_check = float(sys.argv[1])
SN_name = sys.argv[2]
obsid = sys.argv[3]
instrument = sys.argv[4] #pn m1 m2

project = 'projekt_3'

def flare_check_func(rate_limit_check, SN_name, obsid, instrument):
    file = fits.open(f'/Users/juliaahlvind/Documents/{project}/data/xmm/{obsid}/{obsid}_{instrument}_bkgrate.fits')
    hdu= file[1].header
    data = file[1].data

    time = np.arange(hdu['TSTART'],hdu['TSTOP'],hdu['TIMEDEL'])

    # Estimate fallout time
    # –––––––––––––––––––––––
    removed_time = []
    removed_flux = []
    count = 0
    for i in range(len(data)):
        if data['RATE'][i] >= rate_limit_check:
            removed_time.append(time[i])
            removed_flux.append(data['RATE'][i])
        else:
            count += hdu['TIMEDEL']

    obsid_s = obsid
    plt.figure(1)
    plt.title(str(SN_name)+", "+str(obsid_s)+"effective time:"+str(count/1000)+"ks")
    plt.plot(time/1000, data['RATE'], label='Ratio of BKG to SRC?')
    plt.plot(np.divide(removed_time, 1000), removed_flux, '.', label='points above limit')
    plt.axhline(y=rate_limit_check, linestyle='dashed', label='RATE='+str(rate_limit_check), color='black')
    plt.xlabel('ks')
    plt.ylabel('RATE')
    plt.legend()

    
    print('Total time: ', (hdu['TSTOP']-hdu['TSTART'])/1000, 'kS')
    print('Effective time: ', count/1000, 'kS')
    
    pn_path = f"/Users/juliaahlvind/Documents/{project}/data/xmm/{obsid}/"
    plt.savefig(f"/Users/juliaahlvind/Documents/{project}/data/xmm/lightcurves/LightCurve_bkg_flare_{SN_name}_{obsid}_{instrument}.png")
    
    # save effective time in separate txt file
    f = open(f"{pn_path}/exp_time_{SN_name}_{obsid_s}.txt", '+w')
    f.write("total_time(ks) \t effective_time(ks) \n")
    f.write(str((hdu['TSTOP']-hdu['TSTART'])/1000)+ "\t")
    f.write(str(count/1000)+ "\n")
    #plt.show()
    f.close()
    #return plt.show(), count/1000, f
    return count/1000, f

flare_check_func(rate_limit_check,SN_name,obsid,instrument)





