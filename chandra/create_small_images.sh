#! /bin/bash

sn_name=$1
obsid=$2
xx=$3
yy=$4
rad_phys=$5

cd "/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}/repro"


dmcopy "image_058_bin1.fits[sky=circle(${xx},${yy},${rad_phys})]" "small_image_058_bin1_${sn_name}.fits" clobber=yes   
    
cd "/Users/juliaahlvind/Documents/projekt_3/pipelines/chandra"

