#!/bin/bash

input_file="$1"
output_file="${input_file}.fixed"

echo "Processing: $input_file"

# Much simpler and faster AWK script
awk '
BEGIN {
    traj_line = 0
    th_col = 0
}

# Header line
$1 == "time" {
    for (i = 1; i <= NF; i++) {
        if ($i == "TH") th_col = i
    }
    print $0 "         dTH"
    next
}

# Data lines
NF >= 10 && $1 ~ /^[0-9]/ {
    # Store current line and TH value
    current_line = $0
    current_th = (th_col > 0) ? $th_col : 0
    
    # Read next line to get next TH value
    if ((getline next_line) > 0) {
        if (split(next_line, next_fields) >= th_col && th_col > 0) {
            next_th = next_fields[th_col]
            dth_val = sprintf("%.3f", current_th - next_th)
        } else {
            dth_val = "0.000"
        }
    } else {
        dth_val = "0.000"
    }
    
    # Calculate time
    time_val = sprintf("-%.2f", traj_line * 6.0)
    if (traj_line == 0) time_val = "0.00"
    
    # Output modified line
    sub(/^[0-9.-]+/, time_val, current_line)
    print current_line "   " dth_val
    
    # Handle the next line we read
    if (next_line != "") {
        traj_line++
        if (traj_line >= 29) traj_line = 0
        
        # Process next line
        time_val = sprintf("-%.2f", traj_line * 6.0)
        if (traj_line == 0) time_val = "0.00"
        
        sub(/^[0-9.-]+/, time_val, next_line)
        print next_line "   0.000"  # Last line gets 0.000
        
        traj_line++
        if (traj_line >= 29) traj_line = 0
    }
    next
}

# Non-data lines
{
    traj_line = 0  # Reset counter at trajectory boundaries
    print $0
}
' "$input_file" > "$output_file"