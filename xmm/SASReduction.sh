#! /bin/bash

# OBS, måste initiate SAS först med . SAS3.sh !!!!

obs_id=$1 # 0882400101
name=$2 # 2018bsz
preprocess=$3 # true or false

folder="projekt_3" 

home_dir="/Users/juliaahlvind/Documents/projekt_3/data/xmm"
pipeline_dir="/Users/juliaahlvind/Documents/projekt_3/pipelines/xmm"
    
pn_flare_limit=0.4
min_cts=1
oversample=3
soft=300
hard=12000


function prep {
    sas_first_steps
    check_flares
    filter_event_list
    make_image
    filter_event_list_05_10
    make_image_05_10
}


function spec_setup {
    setup_sas
    make_EPIC_pn_spectra_bin1
}

function setup_sas {
    echo "setting up SAS"
    export SAS_CCFPATH="/Users/juliaahlvind/SAS_22.1.0/ccf"
    export SAS_CCF="$home_dir/$obs_id/ccf.cif"
    export SAS_ODF=$(ls "$home_dir/$obs_id/"*SUM.SAS)
}

function sas_first_steps {

    export SAS_ODF="$home_dir/$obs_id"
    export SAS_CCFPATH="/Users/juliaahlvind/SAS_22.1.0/ccf"

    cd "/Users/juliaahlvind/Documents/${folder}/data/xmm"
    cd $obs_id

    cifbuild fullpath=yes

    export SAS_CCF="$home_dir/$obs_id/ccf.cif"

    odfingest
    export SAS_ODF=$(ls "$home_dir/$obs_id/"*SUM.SAS)

    echo "SAS_ODF:${SAS_ODF}"
    echo "SAS_CCFPATH:${SAS_CCFPATH}"
    echo "SAS_CCF:${SAS_CCF}"

    epchain
    #emchain
    
    mv *PIEVLI* raw_pn.fits
    #mv *M1*MIEVLI* raw_m1.fits
    #mv *M2*MIEVLI* raw_m2.fits
    
    mkdir chain/
    mv *.FIT chain/
    
}

run_flare_check() {
    local instrument=$1
    local default_rate=$2

    echo ""
    echo "Checking flares for instrument: $instrument"

    # Initial run with default rate
    python "$pipeline_dir/flare_check.py" "$default_rate" "$name" "$obs_id" "$instrument"

    # Prompt for new rate
    #read -u 1 -p "Enter new threshold rate for $instrument [default: $default_rate]: " user_rate
    #local rate="${user_rate:-$default_rate}"

    # Re-run if changed
    #if [[ "$rate" != "$default_rate" ]]; then
    #    python "$pipeline_dir/flare_check.py" "$rate" "$name" "$obs_id" "$instrument"
    #fi
}

function check_flares {
    cd "/Users/juliaahlvind/Documents/${folder}/data/xmm/${obs_id}"
    
    evselect table=raw_pn.fits withrateset=yes rateset=${obs_id}_pn_bkgrate.fits timecolumn=TIME timebinsize=50 makeratecolumn=yes maketimecolumn=yes expression="(PATTERN == 0)&&(#XMMEA_EP)&&(PI IN [10000:12000])"
    
    run_flare_check "pn" $pn_flare_limit
    
    tabgtigen table=${obs_id}_pn_bkgrate.fits gtiset=${obs_id}_pn_gti.fits expression="RATE<$pn_flare_limit" #typical 0.4
}

function filter_event_list {
    cd "/Users/juliaahlvind/Documents/${folder}/data/xmm/${obs_id}"

    evselect table=raw_pn.fits withfilteredset=true filteredset=${obs_id}_pn_filt_gti.fits keepfilteroutput=true destruct=true expression="(gti(${obs_id}_pn_gti.fits,TIME)&&(PI IN [$soft:$hard])&&(PATTERN<=4)&&(FLAG==0))"

    evselect table=raw_pn.fits withfilteredset=true filteredset=${obs_id}_pn_filt.fits keepfilteroutput=true destruct=true expression="((PI IN [$soft:$hard])&&(PATTERN<=4)&&(FLAG==0))"
    }

function make_image {
    evselect table=${obs_id}_pn_filt_gti.fits withimageset=true imageset=${obs_id}_pn_img_filt_gti.fits xcolumn=X ycolumn=Y imagebinning=binSize ximagebinsize=80 yimagebinsize=80
  }

function filter_event_list_05_10 {
    evselect table=raw_pn.fits withfilteredset=true filteredset=${obs_id}_pn_filt_gti_05_10.fits keepfilteroutput=true destruct=true expression="(gti(${obs_id}_pn_gti.fits,TIME)&&(PI IN [500:10000])&&(PATTERN<=4)&&(FLAG==0))"
  
}

function make_image_05_10 {
    evselect table=${obs_id}_pn_filt_gti_05_10.fits withimageset=true imageset=${obs_id}_pn_img_filt_gti_05_10.fits xcolumn=X ycolumn=Y imagebinning=binSize ximagebinsize=80 yimagebinsize=80
   }



function make_EPIC_pn_spectra_bin1 {
    cd "/Users/juliaahlvind/Documents/${folder}/data/xmm/${obs_id}"

    
    if [ -f "src_${name}_${obs_id}.reg" ]; then
        echo "Creating EPIC-pn spectra"
        
        # source and background region files created in python and manually from ${obs_id}_pn_img_filt_gti_05_10.fits in ds9 and saved in physical coordiantes
        src_file="/Users/juliaahlvind/Documents/${folder}/data/xmm/${obs_id}/src_${name}_${obs_id}.reg" 
        pn_src=$(grep -v '^\s*$' "$src_file" | tail -n 1)
        echo "src reg: ${src_file} och line is: ${pn_src}"
        
        bkg_file="/Users/juliaahlvind/Documents/${folder}/data/xmm/${obs_id}/bkg_annulus_${name}_${obs_id}.reg"
        pn_bkg=$(grep -v '^\s*$' "$bkg_file" | tail -n 1)
        echo "bkg reg: ${bkg_file} och line is: ${pn_bkg}"
            
        #Source spectra
        evselect table=${obs_id}_pn_filt_gti_05_10.fits updateexposure=yes withspecranges=yes withspectrumset=yes energycolumn=PI expression="(X,Y) IN $pn_src" specchannelmax=20479 specchannelmin=0 spectralbinsize=1 spectrumset=${name}_${obs_id}_pn_spec_src_bin1.fits
        
        backscale spectrumset=${name}_${obs_id}_pn_spec_src_bin1.fits badpixlocation=${obs_id}_pn_filt_gti_05_10.fits withbadpixcorr=yes useodfatt=no

        #Background spectra
        evselect table=${obs_id}_pn_filt_gti_05_10.fits updateexposure=yes withspecranges=yes withspectrumset=yes energycolumn=PI expression="(X,Y) IN $pn_bkg" specchannelmax=20479 specchannelmin=0 spectralbinsize=1 spectrumset=${name}_${obs_id}_pn_spec_bkg_bin1.fits

        backscale spectrumset=${name}_${obs_id}_pn_spec_bkg_bin1.fits badpixlocation=${obs_id}_pn_filt_gti_05_10.fits withbadpixcorr=yes useodfatt=no
  
        #Create response files
        rmfgen rmfset=${name}_${obs_id}_pn_rmf_bin1.fits spectrumset=${name}_${obs_id}_pn_spec_src_bin1.fits

        arfgen arfset=${name}_${obs_id}_pn_arf_bin1.fits spectrumset=${name}_${obs_id}_pn_spec_src_bin1.fits withrmfset=yes rmfset=${name}_${obs_id}_pn_rmf_bin1.fits badpixlocation=${obs_id}_pn_filt_gti_05_10.fits

        #Group spectrum and associate files
        specgroup spectrumset=${name}_${obs_id}_pn_spec_src_bin1.fits groupedset=${name}_${obs_id}_pn_spec_grp_bin1.fits mincounts=$min_cts oversample=$oversample rmfset=${name}_${obs_id}_pn_rmf_bin1.fits arfset=${name}_${obs_id}_pn_arf_bin1.fits backgndset=${name}_${obs_id}_pn_spec_bkg_bin1.fits

    fi
}



if [[ "$preprocess" == "true" || "$preprocess" == "True" ]]; then
        echo "Preprocessing the data"
        prep
else
        echo "No preprocess, spectra are generated"
        spec_setup

fi


     



