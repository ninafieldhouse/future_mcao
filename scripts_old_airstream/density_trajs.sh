#!/bin/bash

module load dyn_tools

regions=("JS" "NA")
# regions=("JS")
timeframes=("present" "future")
# timeframes=("present")
airstreams=("ALL" "BCAR" "LL" "DI" "nonBCAR")
# airstreams=("ALL")

base_input_path="/net/litho/atmosdyn2/nfieldhouse/data/lagranto"
for region in "${regions[@]}"; do
    echo "Processing region: $region"
    for timeframe in "${timeframes[@]}"; do
        echo "Processing timeframe: $timeframe"
        for airstream in "${airstreams[@]}"; do
            TRAJECTORY_FILE="CAO_${airstream}_trajs_${region}_${timeframe}_ALL_with_dTH.txt"
            OUTPUT_FILE_1="density_all_CAO_${airstream}_trajs_${region}_${timeframe}_ALL_with_dTH"
            OUTPUT_FILE_2="density_timesteps_CAO_${airstream}_trajs_${region}_${timeframe}_ALL_with_dTH"
            
            echo "Processing $airstream..."
            # Define the fields to process
            FIELDS=("dTH" "dTH_24" "p" "TH" )

            # Define specific time steps (in hours)
            SPECIFIC_TIMES=(0.00 -24.00 -48.00 -72.00 -96.00 -120.00 -144.00 -168.00)

            echo "Starting density gridding process..."

            # Step 1: Create file with all fields at all trajectory times
            # Create file with dTH field at all times
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_1} -field dTH

            # Add dTH_24 field at all times
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_1} -field dTH_24

            # Add p field at all times
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_1} -field p

            # Add TH field at all times
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_1} -field TH

            
            # Step 2: Create file with specific time steps
            # Create file with dTH field at each specific time step
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time 0.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -24.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -48.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -72.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -96.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -120.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -144.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH -time -168.00

            # Add dTH_24 field at each specific time step
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time 0.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -24.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -48.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -72.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -96.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -120.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -144.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field dTH_24 -time -168.00
            
            # Add p field at each specific time step
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time 0.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -24.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -48.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -72.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -96.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -120.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -144.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field p -time -168.00

            # Add TH field at each specific time step
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time 0.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -24.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -48.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -72.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -96.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -120.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -144.00
            density ${TRAJECTORY_FILE} ${OUTPUT_FILE_2} -field TH -time -168.00


            echo "Density gridding complete!"
            echo "Output file: $OUTPUT_FILE_1 and $OUTPUT_FILE_2"
            echo "Fields included: ${FIELDS[*]}"
            echo "Available at all trajectory times AND at specific times: ${SPECIFIC_TIMES[*]}"
        done
    done
done