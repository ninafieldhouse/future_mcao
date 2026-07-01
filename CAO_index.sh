#!/bin/bash

cd /home/nfieldhouse/ma/run/cesm  # temp run directory

module load dyn_tools  # load tools

# Base output directory
member_id=$1  # Get from command line
output_base="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/clim"
mkdir -p "$output_base"  # Ensure base directory exists

# Define regions and corresponding longitude-latitude boxes
declare -A lonlat_boxes
lonlat_boxes[NA]="-sellonlatbox,-75,-65,30,45"
lonlat_boxes[JS]="-sellonlatbox,130,150,30,60" #-15, +30, -10, 0

# Define timeframes, corresponding year ranges, and data sources
declare -A year_ranges
declare -A data_sources

year_ranges[present]="1990 2000"
year_ranges[future]="2090 2100"

data_sources[present]="/net/meso/climphys/cesm212/b.e212.BHISTcmip6.f09_g17.${member_id}/archive/atm/hist"
data_sources[future]="/net/meso/climphys/cesm212/b.e212.BSSP370cmip6.f09_g17.${member_id}/archive/atm/hist"

# Define regions and timeframes
regions=("JS" "NA")
timeframes=("present" "future")

# Loop through regions and timeframes
for region in "${regions[@]}"; do
    echo "Processing region: $region"
    # Define output directory (same for both present and future)
    output_dir="${output_base}/${region}"
    mkdir -p "$output_dir"  # Ensure directory exists

    for timeframe in "${timeframes[@]}"; do
        echo "  Timeframe: $timeframe"
        # Get longitude-latitude box, year range, and data source
        lonlat_box=${lonlat_boxes[$region]}
        read start_year end_year <<< ${year_ranges[$timeframe]}
        data_path=${data_sources[$timeframe]}

        # Set correct h3 and h1 filenames based on timeframe
        if [[ "$timeframe" == "present" ]]; then
            h3="b.e212.BHISTcmip6.f09_g17.${member_id}.cam.h3."
            h1="b.e212.BHISTcmip6.f09_g17.${member_id}.cice.h1."
        else
            h3="b.e212.BSSP370cmip6.f09_g17.${member_id}.cam.h3."
            h1="b.e212.BSSP370cmip6.f09_g17.${member_id}.cice.h1."
        fi

        # Process files for the given year range
        for ((yy = start_year; yy <= end_year; yy++)); do
            echo "    Processing year: $yy"
            # Calculate variables
            cdo $lonlat_box -select,name=T,levidx=26 ${data_path}/${h3}${yy}-01-01-00000.nc T${yy}
            cdo expr,'TH=T*(1000/(85953.5/100))^0.286;' T${yy} TH${yy}

            # Calculate potential SST
            cdo $lonlat_box -selname,SST,PS ${data_path}/${h3}${yy}-01-01-00000.nc SST${yy}
            cdo setrtomiss,0.1,1 SST${yy} SST${yy}_masked
            cdo expr,'SSTH=SST*(1000/(PS/100))^0.286;' SST${yy}_masked SSTH${yy}

            # Mask out areas with sea ice
            cdo setrtomiss,0,273 SST${yy}_masked CI${yy}_masked
            cdo setrtoc,273,350,1 CI${yy}_masked CI${yy}_mask

            # Calculate differences
            cdo merge TH${yy} SSTH${yy} CI${yy}_mask M${yy}
            cdo expr,'CAOI=(TH-SSTH)*SST;' M${yy} C${yy}

            # Move final file to output directory
            mv C${yy} $output_dir/C${yy}
            mv SST${yy} $output_dir/SST${yy}
            mv CI${yy}_mask $output_dir/CI${yy}_mask

            # Remove temporary files
            rm CI${yy}* SSTH${yy}* TH${yy}* M${yy}* T${yy}*
        done
    done
done