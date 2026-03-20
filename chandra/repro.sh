#! /bin/bash

obsid=$1
full_obsid=$2

input_dir="/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}"
output_dir="${input_dir}/repro"
chandra_repro "${input_dir}" outdir="${output_dir}" clobber=yes


cd "/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}/repro"

# Apply a 0.3 to 10 keV filter
dmcopy infile="acisf${full_obsid}_repro_evt2.fits[energy=300:10000]" outfile=event_file_0310.fits  clobber=yes

# Create background light curve
dmextract infile="event_file_0310.fits[bin time=::200]" outfile='lightcurve_0310.fits' opt=ltc1  clobber=yes

# Find GTI
deflare lightcurve_0310.fits lightcurve_0310.gti method=clean save=deflare.png

# plott background light curve
python "/Users/juliaahlvind/Documents/projekt_3/pipelines/lightcurve_chandra.py" $obsid

# Apply GTI on fits
dmcopy infile="event_file_0310.fits[@lightcurve_0310.gti]" outfile="acisf${obsid}_repro_evt2_03_10_clean.fits" clobber=yes

# Create images
dmcopy "acisf${obsid}_repro_evt2_03_10_clean.fits[bin X=1,Y=1][energy=300:10000]" image_0310_bin1.fits clobber=yes
dmcopy "acisf${obsid}_repro_evt2_03_10_clean.fits[bin X=1,Y=1][energy=500:8000]" image_058_bin1.fits clobber=yes


punlearn ardlib


echo "remove folders primary and secondary"
rm -rf "${input_dir}/primary"
rm -rf "${input_dir}/secondary"
