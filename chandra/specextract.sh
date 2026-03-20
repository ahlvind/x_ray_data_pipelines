#! /bin/bash

sn_name=$1
obsid=$2

root_path="/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}/repro"
cd $root_path

src="psfsize_src_${sn_name}_${obsid}.reg"
bkg="bkg_annulus_${sn_name}_${obsid}.reg"


# creates a spectrum, rmf, arf, bkg
specextract \
  infile="acisf${obsid}_repro_evt2_03_10_clean.fits[sky=region(${src})]" \
  bkgfile="acisf${obsid}_repro_evt2_03_10_clean.fits[sky=region(${bkg})]" \
  grouptype=NUM_CTS \
  binspec=1 \
  outroot="spectra_${sn_name}_${obsid}" \
  clobber=yes \
  weight=no \
  weight_rmf=no \
  correctpsf=yes \
  &> specextract.log

if grep -q "has zero counts" specextract.log; then
    echo "Zero counts detected"
    exit 42
fi
# output: spectra_snname_obsid_grp.pi


