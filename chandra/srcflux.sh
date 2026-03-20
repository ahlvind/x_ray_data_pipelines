#! /bin/bash
sn_name=$1
obsid=$2


punlearn srcflux

cd "/Users/juliaahlvind/Documents/projekt_3/data/chandra/${obsid}/repro"


src_reg_file="psfsize_src_${sn_name}_${obsid}.reg"
bkg_reg_file="bkg_annulus_${sn_name}_${obsid}.reg"


# Läs in raden med ellipse och extrahera RA och DEC
line=$(grep '^circle' ${src_reg_file})
# Ta ut RA och DEC med sed
if [[ $line =~ circle\(([0-9:.+-]+),([0-9:.+-]+), ]]; then
    RA="${BASH_REMATCH[1]}"
    DEC="${BASH_REMATCH[2]}"
fi
    echo "Problem with identifying RA, DEC"

echo "RA,DEC: ${RA} ${DEC}"

prefix='./'
fovfile_path=$(find . -name "*repro_fov1.fits*")
fov_file=${fovfile_path/#$prefix}
asolfile_path=$(find . -name '*asol1.fits*')
asol_file=${asolfile_path/#$prefix}
mskfile_path=$(find . -name '*_msk1.fits*')
msk_file=${mskfile_path/#$prefix}
bpixfile_path=$(find . -name '*repro_bpix1.fits*')
bpix_file=${bpixfile_path/#$prefix}
event_file="acisf${obsid}_repro_evt2_03_10_clean.fits"

if [[ -f "spectra_${sn_name}_${obsid}.corr.arf" ]]; then
    echo "arf exists -> calc effective energy"
    energ_lo=0.5
    energ_hi=8.0
    arf="spectra_${sn_name}_${obsid}.corr.arf"
    dmtcalc \
        $arf \
        arf_weights \
        expression="mid_energy=(${energ_lo}+${energ_hi})/2.0;weights=(mid_energy*specresp)" clob+

    dmstat \
        "arf_weights[mid_energy=${energ_lo}:${energ_hi}][cols weights,specresp]" verbose=0

    res=$(pget dmstat out_sum)
    a=${res%,*}
    b=${res#*,}
    eff_en=$(python -c "print($a/$b)")
    echo "effektiv energi är:${eff_en}"

else
    echo "arf does not exist -> assume 4"
    eff_en=3
fi

# include rmf and arf if spectrum was generated
if [[ -f "spectra_${sn_name}_${obsid}.corr.arf" ]]; then

    arf="spectra_${sn_name}_${obsid}.corr.arf"
    rmf="spectra_${sn_name}_${obsid}.rmf"

    echo "${src_reg_file} är källan"

    srcflux \
        $event_file \
        "${RA},${DEC}" \
        outroot=rhooph05_8_${sn_name} \
        srcreg="region(${src_reg_file})" \
        bkgreg="region(${bkg_reg_file})" \
        conf=0.997 \
        bands=${energ_lo}:${energ_hi}:${eff_en} \
        binsize=1 \
        rmffile="$rmf" \
        arffile="$arf" \
        fovfile=$fov_file \
        asolfile=$asol_file \
        mskfile=$msk_file \
        bpixfile=$bpix_file \
        clobber=yes \
     
else

    srcflux \
        $event_file \
        "${RA},${DEC}" \
        outroot=rhooph05_8_${sn_name} \
        srcreg="region(${src_reg_file})" \
        bkgreg="region(${bkg_reg_file})" \
        conf=0.997 \
        bands=${energ_lo}:${energ_hi}:${eff_en} \
        binsize=1 \
        fovfile=$fov_file \
        asolfile=$asol_file \
        mskfile=$msk_file \
        bpixfile=$bpix_file \
        clobber=yes  \
      

fi

mv rhooph05_8_${sn_name}_summary.txt "srcflux_output_05_8_${sn_name}_${obsid}.txt"

rm *rhooph*