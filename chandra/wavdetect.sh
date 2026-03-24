#! /bin/bash

sn_name=$1
obsid=$2

cd "/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}/repro"

# create psf map
echo "create psfmap if it does not exsitst"
FILE="mypsfmap_058_${sn_name}_${obsid}.fits"
if [ -e "$FILE" ]
then
    echo "Exposuremap already exists"
else
    echo "making exposuremap"
    mkpsfmap infile="small_image_058_bin1_${sn_name}.fits" outfile="mypsfmap_058_${sn_name}_${obsid}.fits" energy=1.49 ecf=0.90 clobber=yes
fi

fits_file="image_058_bin1.fits"
livetime_full=$(fitsheader "$fits_file" | grep "LIVETIME" | awk -F '=' '{print $2}' | tr -d ' ')
livetime=$(echo $livetime_full | sed 's/[^0-9.E+-]//g')
echo $livetime

# run wavdetect for given interval
set clobber=yes
wavdetect \
    infile="small_image_058_bin1_${sn_name}.fits" \
    outfile="sources_058_wavdetect_${sn_name}.fits" \
    imagefile="image_058_wavdetect_${sn_name}.fits" \
    defnbkgfile="background058_${sn_name}.fits" \
    scellfile="new_sources_058_${sn_name}.fits" \
    scales="1.0 2.0 3.0" \
    exptime=$livetime \
    psffile="mypsfmap_058_${sn_name}_${obsid}.fits" \
    regfile="srcs_wavdetect_058_${sn_name}.reg" \
    sigthresh=1e-5 \
    clobber=yes


