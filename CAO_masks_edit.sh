#!/bin/bash
module unload dyn_tools
module load dyn_tools

cd /home/nfieldhouse/ma/run/          # change to your directory

# Define regions and timeframes
regions=("JS" "NA")
# regions=("NA")
timeframes=("present" "future")
# timeframes=("future")  # Uncomment to test with only the present timeframe
member_id=$1  # Get from command line

# Loop through regions and timeframes
for region in "${regions[@]}"
do
        # region="NA"  # Set region to NA for testing
        dircao="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d"
        mkdir -p ${dircao}  # Ensure CAO data directory exists
        dirout="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/lagranto/cao_input/${region}"
        mkdir -p ${dirout}  # Ensure output directory exists
        #     file_events_present="${dircao}/summary_events_low_CAOI_DJF10_overstep_5_present_${region}"
        #     file_events_future="${dircao}/summary_events_low_CAOI_DJF10_overstep_5_future_${region}"

        for timeframe in "${timeframes[@]}"
        do
                # Define paths and filenames based on region and timeframe
                file_events="${dircao}/summary_events_DJF90_${timeframe}_overstep_24_CAOI_${region}_filtered"

                # Create output directory if it doesn't exist
                mkdir -p ${dirout}

                nlines=$(wc -l $file_events | awk '{print $1}')
                echo $nlines

                for ((i=1; i<=$nlines; i++)); do            # loops through all dates of your events file
                        # read date and value from datelist
                        date=$(awk -v line=$i 'NR==line{print $1}' $file_events)
                        value=$(awk -v line=$i 'NR==line{print $4}' $file_events)
                        echo $date
                        echo $value

                        year=$(echo -e ${date} | cut -c1-4)
                        month=$(echo -e ${date} | cut -c5-6)

                        mkdir -p ${dircao}/${year}  # create year directory in CAO data directory
                        mkdir -p ${dirout}/${year}  # create year directory in output directory

                        echo "Region: ${region}, Timeframe: ${timeframe}, Output Directory: ${dirout}"

                        python /home/nfieldhouse/ma/scripts/CAO_masks_2.py ${date} ${dirout}/${year} ${region} ${timeframe} ${member_id} # run python script for input date (makes the CAO mask)
                        # run first part of python script to get first coa mask
                done
        done     
done
