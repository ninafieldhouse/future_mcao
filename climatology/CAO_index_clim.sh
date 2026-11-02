#!/bin/bash

source "$(dirname "$0")/../config/config.sh"

# Define base directory where C-files are stored
member_id=$1  # Get from command line
output_base="${CAO_DATA_ROOT}/${member_id}/cao/clim"
mkdir -p "$output_base"  # Ensure base directory exists

# Define regions and their corresponding directories
regions=("JS" "NA")

# Define year ranges for present and future
declare -A year_ranges
year_ranges[present]="1990 2000"
year_ranges[future]="2090 2100"

# Loop through each region
for region in "${regions[@]}"; do
    output_dir="${output_base}/${region}"
    
    # Process present and future timeframes
    for timeframe in "${!year_ranges[@]}"; do
        read start_year end_year <<< ${year_ranges[$timeframe]}

        # Change to the correct output directory
        cd "$output_dir" || exit 1  # Exit if directory doesn't exist

        # Ensure at least one file exists before proceeding
        if [ -z "$(ls C* 2>/dev/null)" ]; then
            echo "No C-files found in $output_dir. Skipping..."
            continue
        fi

        # Set timestamps for each file
        for ((yy = start_year; yy <= end_year; yy++)); do
            if [ -f "C${yy}" ]; then
                cdo settaxis,${yy}-01-01,00:00:00,6h C${yy} C${yy}t
            else
                echo "Warning: File C${yy} not found in $output_dir"
            fi
        done

        # Calculate hourly mean climatology
        cdo -yhourmean -mergetime C*t "yearly_hour_mean_${timeframe}"

        # Merge all files to calculate anomaly
        cdo mergetime C*t "all_times_${timeframe}"

        # Calculate anomalies
        cdo -yhoursub "all_times_${timeframe}" "yearly_hour_mean_${timeframe}" "yearly_hour_anomalies_${timeframe}"

        # Remove temporary C*t files
        rm C*t
    done
done
