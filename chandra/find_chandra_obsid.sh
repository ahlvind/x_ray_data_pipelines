#! /bin/bash



cross_match_chandra() {
    local input_file="/Users/juliaahlvind/Documents/projekt_1/OSNC/OSNC_60MPc.txt" 
    local output_file="/Users/juliaahlvind/Documents/projekt_3/sample/raw_cross_matches_optical_chandra.txt"

    # Check if the input file excists
    if [[ ! -f "$input_file" ]]; then
        echo "Error: Inputfile '$input_file' does not excist." >&2
        return 1
    fi

    # Clear previous output file before printing
    : > "$output_file"

    # Read row by row 
    while read -r Name Type dist_Mpc RA DEC Date_of_discovery; do
        
        # Korrigera variabelnamn (RA/DEC verkar komma från $ra1/$dec1 i ditt original)
        local ra_value="$RA"
        local dec_value="$DEC"

        # If Name starts with 'SN'+number, remove 'SN'
        # E.g. SN1987A, SN2021abc, SN3, etc.
        if [[ "$Name" =~ ^SN[0-9] ]]; then
            clean_name="${Name#SN}"
        else
            clean_name="$Name"
        fi

        # Print Name to output file
        echo "$clean_name" >> "$output_file"

        # Run Chandra process find_chandra_obsid
        # instrument=acis includes all ACIS chips
        find_chandra_obsid "$ra_value" "$dec_value" instrument=acis >> "$output_file"

    done < "$input_file"
}


cross_match_chandra 