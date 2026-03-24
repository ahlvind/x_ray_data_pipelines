#! /bin/bash

# OBS, måste initiate SAS först med . SAS3.sh !!!!

obs_id=$1 # 0882400101
name=$2 # 2018bsz

folder="projekt_3" 

home_dir="/Users/juliaahlvind/Documents/projekt_3/data/xmm"
pipeline_dir="/Users/juliaahlvind/Documents/projekt_3/pipelines/xmm"
    

function setup_sas {
    echo "setting up SAS"
    export SAS_CCFPATH="/Users/juliaahlvind/SAS_22.1.0/ccf"
    export SAS_CCF="$home_dir/$obs_id/ccf.cif"
    export SAS_ODF=$(ls "$home_dir/$obs_id/"*SUM.SAS)
}


function _ewavelet {
    
    
    cd "/Users/juliaahlvind/Documents/${folder}/data/xmm/${obs_id}"

    att_file=$(find ./chain -type f -name '*ATT*' -print -quit)

    # Skriv ut första t.ex.
    echo "ATT file: ${att_file}"
    
    
    FILE="expmap_05_10_pn.fits"
    if [ -e "$FILE" ]
    then
        echo "Exposuremap already exists"
    else
        # create exposure map
        echo "exposuremap needs to be created"
        eexpmap imageset=${obs_id}_pn_img_filt_gti_05_10.fits expimageset=expmap_05_10_pn.fits attitudeset=$att_file eventset=${obs_id}_pn_filt_gti.fits withvignetting=true

    fi

    # extract detected sources
    echo "running ewavelet"
    ewavelet imageset=${obs_id}_pn_img_filt_gti_05_10.fits expmapset=expmap_05_10_pn.fits srclistset=source_list05_10_pn_3σ.fits minscale=1 maxscale=10 threshold=3 makerecon=no edgethreshold=2.0

}


setup_sas
_ewavelet

cd $pipeline_dir