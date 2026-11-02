#!/bin/bash
source "$(dirname "$0")/../config/config.sh"

member_id=$1  # Get from command line, e.g., "member_001"
regions=("JS" "NA")

for region in "${regions[@]}"; do
INPUT_DIR="${CAO_DATA_ROOT}/${member_id}/lagranto/traj/${region}"
for INPUT_FILE in "$INPUT_DIR"/B_*; do
    OUTPUT_FILE="${INPUT_FILE}_with_TH.txt"

    awk '
    BEGIN {
        th_col_name = "TH"
    }
    NR == 3 {
        # Add "TH" to the header line
        print $0, "\t", th_col_name > "'"$OUTPUT_FILE"'"
        next
    }
    NR < 5 {
        # Copy the metadata lines to the output file
        print $0 > "'"$OUTPUT_FILE"'"
        next
    }
    $1 ~ /^-?[0-9]/ {
        # Match rows where the first column (time) is numeric
        if (NF >= 6) {  # Ensure there are enough columns
            p = $4
            T = $6
            if (p > 0 && T != "") {
                # Perform the TH calculation (corrected formula)
                TH = T * (1000 / p) ^ 0.286
                # Append the calculated TH value to the line
                printf "%s %.2f\n", $0, TH >> "'"$OUTPUT_FILE"'"
            } else {
                # If p or T is invalid, append NaN
                print $0 " NaN" >> "'"$OUTPUT_FILE"'"
            }
        } else {
            # If the line does not have enough columns, append NaN
            print $0 " NaN" >> "'"$OUTPUT_FILE"'"
        }
        next
    }
    {
        # Copy any remaining lines to the output file
        print $0 >> "'"$OUTPUT_FILE"'"
    }' "$INPUT_FILE"
done
done