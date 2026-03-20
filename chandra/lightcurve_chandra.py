
import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from lightcurves import *


obsid = sys.argv[1]


# check background flares: bkgrates_tot.fits
# check final source-bakcground lightcurve: src_lc.fits   bkgrates_tot.fits



file = fits.open(f'/Users/juliaahlvind/Documents/projekt_3/data/chandra/{obsid}/repro/lightcurve_0310.fits')
hdu= file[1].header
data = file[1].data

time = np.arange(hdu['TSTART'],hdu['TSTOP'],hdu['TIMEDEL'])
time = np.arange(0,hdu['LIVETIME'],hdu['LIVETIME']/hdu['NAXIS2'])

rate_limit_check = 0.5
# Estimate fallout time
# –––––––––––––––––––––––

saved_time = []
saved_flux = []
count = 0
for i in range(len(data)):
    if data['COUNT_RATE'][i] != 0: #rate_limit_check:
        saved_time.append(time[i])
        saved_flux.append(data['COUNT_RATE'][i])
        count += hdu['TIMEDEL']
mean = np.mean(saved_flux)

removed_time = []
removed_flux = []
for i in range(len(data)):
    if data['COUNT_RATE'][i] > 1.2*mean or data['COUNT_RATE'][i]==0: #rate_limit_check:
        removed_time.append(time[i])
        removed_flux.append(data['COUNT_RATE'][i])
        count += hdu['TIMEDEL']

plt.figure(1)
plt.plot(time/1000, data['COUNT_RATE'], label='Ratio of BKG to SRC?')
plt.plot(np.divide(removed_time, 1000), removed_flux, '.', label='removed points')
plt.axhline(y=mean, linestyle='dashed', label='mean='+str(mean), color='black')
plt.xlabel('ks')
plt.ylabel('COUNT_RATE')
plt.title('Total time: '+ str((hdu['TSTOP']-hdu['TSTART'])/1000) + 'kS, Effective time: '+ str(count/1000)+ 'kS. with scale 1.2')
plt.legend()

plt.show()
obsid_s = obsid.split('/')[-1]
plt.savefig(f'/Users/juliaahlvind/Documents/projekt_3/data/chandra/lightcurves/lightcurve_{obsid_s}.png')
plt.clf()

print('Total time: ', (hdu['TSTOP']-hdu['TSTART'])/1000, 'kS')
print('Effective time: ', count/1000, 'kS')









