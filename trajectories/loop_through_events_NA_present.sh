#!/bin/bash
source "$(dirname "$0")/../config/config.sh"

# Wrapper script to call _cao for all events
member_id=$1  # Get from command line

# Path to the main script (_cao)
main_script="$(dirname "$0")/backward_forward_trajectories_cao_NA_present.sh"

# Path to the file_events file
file_events="${CAO_DATA_ROOT}/${member_id}/cao/events_d/summary_events_DJF90_present_overstep_24_CAOI_NA_filtered"

# Read the number of lines in the file_events
nlines=$(wc -l $file_events | awk '{print $1}')
echo "Number of events: $nlines"

# Loop through all events
for ((i=1; i<=$nlines; i++)); do
    echo "Processing event $i..."
    bash $main_script $i $file_events $member_id 
    # Pass the event number and file_events to the main script
    if [ $? -ne 0 ]; then
        echo "Error: Script failed for event $i. Skipping to the next event."
        continue
    fi
    echo "Event $i processed successfully."
done

echo "All events processed."