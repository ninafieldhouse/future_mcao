#!/bin/bash
member_id=$1  # Get from command line

REGIONS=("JS" "NA")
TIMEFRAMES=("present" "future")
BASE_DIR="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}"
echo "-------------------------------------"
echo "for ensemble member ${member_id}"
echo "-------------------------------------"
for REGION in "${REGIONS[@]}"; do
    for TIMEFRAME in "${TIMEFRAMES[@]}"; do
        echo "Processing ${REGION} for ${TIMEFRAME} time frame"
        DIR="${BASE_DIR}/cao/events_d/"
        TRAJ_DIR="${BASE_DIR}/lagranto/traj/${REGION}"
        DATES_FILE="${DIR}/summary_events_DJF90_${TIMEFRAME}_overstep_24_CAOI_${REGION}_filtered"

        # Count B_ files in the folder
        num_b_files=$(find "${TRAJ_DIR}" -maxdepth 1 -type f -name 'B_*' | wc -l)
        # Count lines in the date file
        num_lines=$(wc -l < "$DATES_FILE")

        echo "${REGION}: number of B_ files in folder = $num_b_files"
        echo "${REGION}: number of lines in CAO_max_intens file = $num_lines"

        missing=()
        present=()
        while read -r line; do
            date_col=$(echo "$line" | awk '{print $1}')
            ymd=$(echo "$date_col" | cut -d'-' -f1-3 | tr -d '-')
            hh=$(echo "$date_col" | awk -F'-' '{print $4}')
            fname="B_${ymd}_${hh}"
            if [ -f "${TRAJ_DIR}/${fname}" ]; then
                present+=("$date_col")
                # echo "${REGION}: present $date_col"
            else
                missing+=("$date_col")
                echo "${REGION}: missing $date_col"
            fi
        done < "$DATES_FILE"
        echo "${REGION}: total present = ${#present[@]}"
        echo "${REGION}: total missing = ${#missing[@]}"
    done
done

REGIONS=("JS" "NA")
TIMEFRAMES=("present" "future")
BASE_DIR="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}"

echo "-------------------------------------"
echo "Checking for B_ files not in CAO events files"
echo "-------------------------------------"
for REGION in "${REGIONS[@]}"; do
    TRAJ_DIR="${BASE_DIR}/lagranto/traj/${REGION}"
    REMOVED_DIR="${BASE_DIR}/removed"
    mkdir -p "$REMOVED_DIR"

    # Collect all valid B_ filenames for both timeframes
    valid_files=()
    for TIMEFRAME in "${TIMEFRAMES[@]}"; do
        DIR="${BASE_DIR}/cao/events_d/"
        DATES_FILE="${DIR}/summary_events_DJF90_${TIMEFRAME}_overstep_24_CAOI_${REGION}_filtered"
        while read -r line; do
            date_col=$(echo "$line" | awk '{print $1}')
            ymd=$(echo "$date_col" | cut -d'-' -f1-3 | tr -d '-')
            hh=$(echo "$date_col" | awk -F'-' '{print $4}')
            fname="B_${ymd}_${hh}"
            valid_files+=("$fname")
        done < "$DATES_FILE"
    done

    # Convert valid_files to a string for grep
    valid_files_str=$(printf "%s\n" "${valid_files[@]}")

    # Now move any B_ file not in the valid list
    for file in "${TRAJ_DIR}"/B_*; do
        fname=$(basename "$file")
        if ! grep -qx "$fname" <<< "$valid_files_str"; then
            echo "Moving $fname to $REMOVED_DIR (not in events files)"
            echo "Would move $file"
        fi
    done
done