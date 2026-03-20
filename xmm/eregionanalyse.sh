#! /bin/bash

obsid=$1
sn_name=$2
x_physical=$3
y_physical=$4

export SAS_ODF="/Users/juliaahlvind/Documents/projekt_3/data/xmm/${obsid}"
export SAS_CCFPATH="/Users/juliaahlvind/XMM_SAS/ccf"
export SAS_CCF="/Users/juliaahlvind/Documents/projekt_3/data/xmm/${obsid}/ccf.cif"

pipeline_dir="/Users/juliaahlvind/Documents/projekt_3/pipelines/xmm"

cd "/Users/juliaahlvind/Documents/projekt_3/data/xmm/${obsid}"


eregionanalyse imageset=${obsid}_pn_img_filt_gti_05_10.fits srcexp="(X,Y) in CIRCLE(${x_physical},${y_physical},100)" psfmodel=ELLBETA ulsig=0.997 > "${sn_name}_${obsid}_eregion_pn_05_10.txt"

cd $pipeline_dir
