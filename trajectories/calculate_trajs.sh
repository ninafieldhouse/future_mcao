#!/bin/bash
source "$(dirname "$0")/../config/config.sh"

member_id=$1  # Get from command line

# Paths to your loop_through_events scripts
SCRIPT_DIR="$(dirname "$0")"
SCRIPT1="${SCRIPT_DIR}/loop_through_events_JS_present.sh"
SCRIPT2="${SCRIPT_DIR}/loop_through_events_JS_future.sh"
SCRIPT3="${SCRIPT_DIR}/loop_through_events_NA_present.sh"
SCRIPT4="${SCRIPT_DIR}/loop_through_events_NA_future.sh"

# Run all scripts in parallel
bash "$SCRIPT1" $member_id &
bash "$SCRIPT2" $member_id &
bash "$SCRIPT3" $member_id &
bash "$SCRIPT4" $member_id &

# Wait for all background jobs to finish
wait

echo "All loop_through_events scripts have completed."