#!/bin/bash
member_id=$1

# Define regions and their input files
declare -A regions=(
  ["JS"]="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d/summary_events_DJF90_present_overstep_24_CAOI_JS_filtered /net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d/summary_events_DJF90_future_overstep_24_CAOI_JS_filtered"
  ["NA"]="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d/summary_events_DJF90_present_overstep_24_CAOI_NA_filtered /net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d/summary_events_DJF90_future_overstep_24_CAOI_NA_filtered"
)

# Base output directory
output_dir="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d"

# Function to process files for a specific region
process_region() {
  local region=$1
  local files=$2
  local output_file="${output_dir}/dates_CAO_max_intens_1990:2100_${region}"

  # Initialize unique event ID counter
  event_id=1

  # Write header to the output file (optional, remove if not needed)
  # echo "# Date-Time Intensity Event-ID Other" > "$output_file"

  # Process each input file for the region
  for file in $files; do
    while read -r line; do
      # Skip empty lines or lines starting with a comment
      [[ -z "$line" || "$line" =~ ^# ]] && continue

      # Extract fields from the line

      date_time=$(echo "$line" | awk '{print $1}')
      intensity=$(echo "$line" | awk '{print $2}')
      other=$(echo "$line" | awk '{print $4}')

      # Write to the output file with a unique event ID
      echo "$date_time $intensity $event_id $other" >> "$output_file"

      # Increment the event ID
      ((event_id++))
    done < "$file"
  done

  # Make the output file executable (optional, remove if not needed)
  chmod +x "$output_file"
}

# Process each region
for region in "${!regions[@]}"; do
  rm -f "${output_dir}/dates_CAO_max_intens_1990:2100_${region}"
  # Create a new output file for the region
  process_region "$region" "${regions[$region]}"
done
