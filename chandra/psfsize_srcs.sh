#! /bin/bash

sn_name=$1
obsid=$2
ra=$3
dec=$4

cd "/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}/repro"

psfsize_srcs "acisf${obsid}_repro_evt2_03_10_clean.fits" "${ra} ${dec}" off_axis_out.fits energy=broad ecf=0.9 clobber=yes

regphystocel off_axis_out.fits "psfsize_src_${sn_name}_${obsid}.reg" wcsfile="acisf${obsid}_repro_evt2_03_10_clean.fits" clobber=yes
    
    

