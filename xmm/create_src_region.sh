#! /bin/bash

obsid=$1
sn_name=$2





def coords2pix(image, ra, de):
    raref = fits.getval(image,'TCRVL11', 1)
    radel = fits.getval(image,'TCDLT11', 1)
    rapix = fits.getval(image,'TCRPX11', 1)
    deref = fits.getval(image,'TCRVL12', 1)
    dedel = fits.getval(image,'TCDLT12', 1)
    depix = fits.getval(image,'TCRPX12', 1)
    xx = (ra-raref)*np.cos(np.deg2rad(de))/radel+rapix-1
    yy = (de-deref)/dedel+depix-1
    return xx, yy

fits_path = glob('{obsid}_pn_filt_gti.fits')[0]
cc = SkyCoord('${ra}', '${de}', frame='icrs')
xx, yy = coords2pix(fits_path, cc.ra.degree, cc.dec.degree)
xx = xx+1
yy = yy+1

ff = open('par.sh', 'w+')
ff.write('src_reg=\"circle(' + str(xx) + ',' + str(yy) + ',16)\"\\\\n')
ff.write('bkg_reg=\"annulus(' + str(xx) + ',' + str(yy) + ',24,48)\"\\\\n')
ff.write('xx=' + str(xx) + '\\\\n')
ff.write('yy=' + str(yy) + '\\\\n')
ff.write('xlo=' + str(xx-20) + '\\\\n')
ff.write('xhi=' + str(xx+20) + '\\\\n')
ff.write('ylo=' + str(yy-20) + '\\\\n')
ff.write('yhi=' + str(yy+20) + '\\\\n')
ff.write('ra_deg=' + str(cc.fk5.ra.degree) + '\\\\n')
ff.write('de_deg=' + str(cc.fk5.dec.degree) + '\\\\n')
ff.close()" > get_par.py
}